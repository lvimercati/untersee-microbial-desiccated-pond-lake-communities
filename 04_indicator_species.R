# =============================================================================
# INDICATOR SPECIES ANALYSIS — Bacteria, Cyanobacteria, Eukarya
# Full script: setup, shared functions, all three panels
# Fixes applied: Desiccated (red, top) -> Pond (blue, middle) -> Lake (green,
# bottom); bold axis/legend text; bigger bold titles; simplified panel titles
# =============================================================================

library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(Hmisc); packageVersion("Hmisc")
library(ggplot2); packageVersion("ggplot2")
library(plotly); packageVersion("plotly")
library(mctoolsr)
library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(RVAideMemoire)
library(car)
library(cowplot)
library(phyloseq)
library(biomformat)
library(ggrepel)
library(scales)
library(ape)
library(picante)
library(phyloseq)
library(mctoolsr)
library(zCompositions)
library(compositions)
library(lmodel2)
library(reshape2)
library(indicspecies)
theme_set(theme_bw())
show_col(hue_pal()(3))

# -----------------------------------------------------------------------
# Guard against plyr masking dplyr::summarise/group_by (same issue we
# fixed in the FAPROTAX script — applies here too, since plyr is loaded)
# -----------------------------------------------------------------------
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

# ===============================================
# SHARED SETTINGS AND FUNCTIONS
# ===============================================

tuft_names <- c("filament 1", "filament 2", "filament 3",
                "filament 4", "filament 5")

REL_ABUND_THRESHOLD <- 0.0005  # 0.05%
N_TOP <- 15

add_final_column <- function(data) {
  data$Final <- apply(data[, c("Pond", "Desiccated", "Lake")], 1, max)
  return(data)
}

get_best_taxonomy <- function(taxonomy_data) {
  taxonomy_data$best_taxonomy <- apply(taxonomy_data, 1, function(row) {
    for (level in c("taxonomy6","taxonomy5","taxonomy4","taxonomy3","taxonomy2","taxonomy1")) {
      if (!is.na(row[level]) && row[level] != "NA" && row[level] != "") return(row[level])
    }
    return(NA)
  })
  return(taxonomy_data)
}

apply_max_env_filter <- function(otu_table, metadata, threshold = 0.001) {
  rel_abund <- sweep(otu_table, 2, colSums(otu_table), "/")
  envs <- unique(metadata$Environment)
  
  env_means <- sapply(envs, function(e) {
    samples <- intersect(rownames(metadata[metadata$Environment == e, ]), colnames(rel_abund))
    if(length(samples) > 1) {
      return(rowMeans(rel_abund[, samples]))
    } else if(length(samples) == 1) {
      return(rel_abund[, samples])
    } else {
      return(rep(0, nrow(rel_abund)))
    }
  })
  
  max_means <- apply(env_means, 1, max)
  keep_asvs <- names(max_means[max_means >= threshold])
  
  cat("ASVs retained (Max Env Mean >= ", threshold * 100, "%):", length(keep_asvs), "\n")
  cat("ASVs removed:", nrow(otu_table) - length(keep_asvs), "\n")
  
  return(otu_table[keep_asvs, ])
}

run_indicator_analysis <- function(otu_table, metadata, env_threshold = 0.5) {
  grouping_var <- metadata[colnames(otu_table), "Environment"]
  indicator    <- multipatt(t(otu_table), grouping_var, control = how(nperm = 999))
  
  sign_results <- rownames_to_column(indicator$sign, var = "ESV")
  sign_results <- sign_results[!is.na(sign_results$p.value), ]
  
  indicator_values <- as.data.frame(indicator$str)
  indicator_values <- rownames_to_column(indicator_values, var = "ESV")
  indicator_values <- indicator_values[, 1:(ncol(indicator_values) - 4)]
  
  merged <- merge(sign_results, indicator_values, by = "ESV")
  significant_results <- merged[merged$p.value < 0.05, ]
  
  significant_results$Environment <- apply(
    significant_results[, c("Desiccated", "Pond", "Lake")], 1, function(row) {
      if (max(row) >= env_threshold) {
        colnames(significant_results[, c("Desiccated","Pond","Lake")])[which.max(row)]
      } else { NA }
    })
  
  significant_results_final <- add_final_column(significant_results)
  significant_results_final <- significant_results_final[, c(1, 11, 12)]
  return(significant_results_final)
}

# ===============================================
# BACTERIA — ALL
# ===============================================
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
cat("\n=== BACTERIA ===\n")

input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')
cat("Bacteria samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

otu_table <- input_filt_rare$data_loaded
metadata  <- input_filt_rare$map_loaded

sample_match <- all(colnames(otu_table) %in% rownames(metadata))
message("Sample IDs match: ", sample_match)

apply_rel_abund_filter <- function(otu_table, threshold = REL_ABUND_THRESHOLD) {
  rel_abund <- sweep(otu_table, 2, colSums(otu_table), "/")
  mean_rel_abund <- rowMeans(rel_abund)
  keep_asvs <- names(mean_rel_abund[mean_rel_abund >= threshold])
  
  cat("ASVs retained (mean rel. abund >=", threshold * 100, "%):", length(keep_asvs), "\n")
  cat("ASVs removed:", nrow(otu_table) - length(keep_asvs), "\n")
  
  return(otu_table[keep_asvs, ])
}

otu_table_filt <- apply_rel_abund_filter(otu_table)

significant_results_final <- run_indicator_analysis(otu_table_filt, metadata,
                                                    env_threshold = 0.5)

taxonomy <- get_best_taxonomy(rownames_to_column(input_filt_rare$taxonomy_loaded, var = "ESV"))
merged_data <- merge(significant_results_final,
                     taxonomy[, c("ESV", "best_taxonomy")],
                     by = "ESV", all.x = TRUE)

merged_data$Environment <- factor(merged_data$Environment,
                                  levels = c("Desiccated", "Pond", "Lake"))

ordered_data <- merged_data %>%
  filter(!is.na(Environment)) %>%
  group_by(Environment) %>%
  slice_max(order_by = Final, n = N_TOP, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Environment, desc(Final))

ordered_data$ESV_label <- paste0(ordered_data$ESV, " (", ordered_data$best_taxonomy, ")")
ordered_data$ESV_label <- factor(ordered_data$ESV_label,
                                 levels = rev(unique(ordered_data$ESV_label)))

x_min_bac <- round(min(ordered_data$Final) - 0.02, 2)
cat("Bacteria x-axis minimum:", x_min_bac, "\n")

pdf("Indicator_species_Bacteria_notufts_final.pdf", width = 6, height = 12)
ggplot(ordered_data, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_min_bac, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5,
               na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(x_min_bac, 1),
                     oob = scales::squish) +
  labs(title = "Bacteria", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 8, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "right",
        legend.title       = element_text(face = "bold"),
        legend.text        = element_text(face = "bold"))
dev.off()

ordered_data_bacteria <- ordered_data

# ===============================================
# CYANOBACTERIA (REVISED FOR LAKE SENSITIVITY)
# ===============================================
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
cat("\n=== CYANOBACTERIA ===\n")

input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')
cyano_only      <- filter_taxa_from_input(input_filt_rare, taxa_to_keep = "Cyanobacteria")

otu_table <- cyano_only$data_loaded
metadata  <- cyano_only$map_loaded

otu_table_filt <- apply_max_env_filter(otu_table, metadata, threshold = 0.0005)
significant_results_final <- run_indicator_analysis(otu_table_filt, metadata, env_threshold = 0.5)

taxonomy <- get_best_taxonomy(rownames_to_column(cyano_only$taxonomy_loaded, var = "ESV"))
merged_data <- merge(significant_results_final,
                     taxonomy[, c("ESV", "best_taxonomy")],
                     by = "ESV", all.x = TRUE)

merged_data$Environment <- factor(merged_data$Environment,
                                  levels = c("Desiccated", "Pond", "Lake"))

ordered_data <- merged_data %>%
  filter(!is.na(Environment), Final >= 0.5) %>%
  group_by(Environment) %>%
  slice_max(order_by = Final, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Environment, desc(Final))

ordered_data$ESV_label <- paste0(ordered_data$ESV, " (", ordered_data$best_taxonomy, ")")
ordered_data$ESV_label <- factor(ordered_data$ESV_label,
                                 levels = rev(unique(ordered_data$ESV_label)))

ordered_data_cyanobacteria <- ordered_data

print(ordered_data[grep("Tychonema", ordered_data$best_taxonomy), ])

pdf("Indicator_species_cyano_notufts_final.pdf", width = 6, height = 6)
ggplot(ordered_data, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.65, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5,
               na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(0.65, 0.95),
                     oob = scales::squish) +
  labs(title = "Cyanobacteria", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 8, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "right",
        legend.title       = element_text(face = "bold"),
        legend.text        = element_text(face = "bold"))
dev.off()

# ===============================================
# EUKARYA
# ===============================================
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
cat("\n=== EUKARYA ===\n")

input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')
cat("Eukarya samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

otu_table <- input_filt_rare$data_loaded
metadata  <- input_filt_rare$map_loaded

sample_match <- all(colnames(otu_table) %in% rownames(metadata))
message("Sample IDs match: ", sample_match)

otu_table_filt <- apply_rel_abund_filter(otu_table)

significant_results_final <- run_indicator_analysis(otu_table_filt, metadata,
                                                    env_threshold = 0.7)

taxonomy <- get_best_taxonomy(rownames_to_column(input_filt_rare$taxonomy_loaded, var = "ESV"))
merged_data <- merge(significant_results_final,
                     taxonomy[, c("ESV", "best_taxonomy")],
                     by = "ESV", all.x = TRUE)
merged_data$Environment <- factor(merged_data$Environment,
                                  levels = c("Desiccated", "Pond", "Lake"))

ordered_data <- merged_data %>%
  filter(!is.na(Environment), Final >= 0.7) %>%
  filter(!grepl("Eukaryota", best_taxonomy, ignore.case = TRUE)) %>%
  group_by(Environment) %>%
  slice_max(order_by = Final, n = N_TOP, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Environment, desc(Final))

ordered_data$ESV_label <- paste0(ordered_data$ESV, " (", ordered_data$best_taxonomy, ")")
ordered_data$ESV_label <- factor(ordered_data$ESV_label,
                                 levels = rev(unique(ordered_data$ESV_label)))

pdf("Indicator_species_Eukarya_notufts_final.pdf", width = 8, height = 8)
ggplot(ordered_data, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.7, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5,
               na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(0.7, 1),
                     oob = scales::squish) +
  labs(title = "Eukarya", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 8, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "right",
        legend.title       = element_text(face = "bold"),
        legend.text        = element_text(face = "bold"))
dev.off()

ordered_data_eukarya <- ordered_data

cat("\nAll three panels generated with fixes applied:\n")
cat("- Desiccated (red) top, Pond (blue) middle, Lake (green) bottom\n")
cat("- Bold axis/legend text, larger bold titles\n")

# =============================================================================
# COMBINE THE THREE PANELS INTO ONE FIGURE (side by side: Bacteria | 
# Cyanobacteria | Eukarya), matching your reference image layout
# =============================================================================

library(magick)

# -----------------------------------------------------------------------
# STEP 1: Rebuild each panel as a ggplot object (not just saved to PDF),
# WITHOUT individual legends, so we can extract ONE shared legend
# (all three use the identical Desiccated/Pond/Lake colour scheme)
# -----------------------------------------------------------------------
panel_bacteria <- ggplot(ordered_data_bacteria, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = x_min_bac, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(x_min_bac, 1), oob = scales::squish) +
  labs(title = "Bacteria", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 12, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "none")

panel_cyano <- ggplot(ordered_data_cyanobacteria, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.65, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(0.65, 0.95), oob = scales::squish) +
  labs(title = "Cyanobacteria", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 12, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "none")

panel_eukarya <- ggplot(ordered_data_eukarya, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.7, xend = Final - 0.015,
                   y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  scale_x_continuous(limits = c(0.7, 1), oob = scales::squish) +
  labs(title = "Eukarya", x = "Indicator Value", y = NULL) +
  theme_bw() +
  theme(axis.text.y        = element_text(size = 12, face = "bold"),
        axis.text.x        = element_text(face = "bold"),
        axis.title.x       = element_text(face = "bold"),
        plot.title          = element_text(size = 16, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "none")

# -----------------------------------------------------------------------
# STEP 2: Build ONE shared legend (Environment: Desiccated/Pond/Lake)
# -----------------------------------------------------------------------
legend_source <- ggplot(ordered_data_bacteria, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4", "Lake" = "#33A02C")) +
  theme(legend.title = element_text(face = "bold", size = 13),
        legend.text  = element_text(face = "bold", size = 11),
        legend.key.size = unit(0.6, "cm"))

legend_grob <- cowplot::get_legend(legend_source)
legend_plot <- ggdraw() + draw_grob(legend_grob)

# -----------------------------------------------------------------------
# STEP 3: Export panels side-by-side and the legend separately
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

panels_row <- panel_bacteria | panel_cyano | panel_eukarya

ggsave("Figure3_panels_only.png", panels_row,
       width = 21, height = 12, dpi = 300, bg = "white")

ggsave("Figure3_legend.png", legend_plot,
       width = 3, height = 4, dpi = 300, bg = "white")

# -----------------------------------------------------------------------
# STEP 4: Composite — legend placed to the right, vertically centered
# -----------------------------------------------------------------------
dpi <- 300
cm_to_px <- function(cm) round(cm * dpi / 2.54)

panels_img <- image_read("Figure3_panels_only.png")
legend_img <- image_trim(image_read("Figure3_legend.png"))

panels_info <- image_info(panels_img)
legend_info <- image_info(legend_img)

gap_cm   <- 1
x_offset <- panels_info$width + cm_to_px(gap_cm)
y_offset <- (panels_info$height - legend_info$height) / 2

canvas_width  <- panels_info$width + legend_info$width + cm_to_px(gap_cm) + cm_to_px(1)
canvas_height <- panels_info$height

final_figure3 <- image_blank(width = canvas_width, height = canvas_height, color = "white")
final_figure3 <- image_composite(final_figure3, panels_img, offset = "+0+0")
final_figure3 <- image_composite(final_figure3, legend_img,
                                 offset = paste0("+", round(x_offset), "+", round(y_offset)))

# -----------------------------------------------------------------------
# STEP 5: Save — PDF and high-res JPEG
# -----------------------------------------------------------------------
image_write(final_figure3, path = "Figure3_Combined_Final.jpeg", format = "jpeg", quality = 100)
image_write(final_figure3, path = "Figure3_Combined_Final.pdf", format = "pdf", density = dpi)

cat("Saved: Figure3_Combined_Final.jpeg and Figure3_Combined_Final.pdf\n")