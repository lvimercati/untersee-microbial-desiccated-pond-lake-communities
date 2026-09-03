# =============================================================================
# COMBINED 16S + 18S FIGURES — Alpha/Beta Diversity (Figure 5) and Phylum-
# Level Composition (Figure 2)
# Requires objects from BOTH 01_16S_processing_and_analysis.R and
# 02_18S_processing_and_analysis.R to be in the environment (or re-derived
# from the rarefied, tuft-excluded RDS objects for both domains, as this
# script does internally where noted).
# =============================================================================

# =============================================================================
# COMBINED 4-PANEL FIGURE: Bacterial & Eukaryotic Alpha and Beta Diversity
# Row 1: Alpha diversity — 16S (A, left) | 18S (B, right)
# Row 2: Beta diversity  — 16S (C, left) | 18S (D, right)
# Column 1 = Bacteria throughout, Column 2 = Eukarya throughout
# Uses the same env_colours palette, theme_bw(), and convex-hull PCoA style
# as all other figures in this script. Requires objects from BOTH the
# bacterial (this script) and eukaryotic (lakeponddesic_eukarya_notufts.R)
# scripts to already be in the environment (meta_bac, meta_euk, micro.hulls.bac,
# micro.hulls.euk, cld_results for both domains).
# =============================================================================
# =============================================================================
# SELF-CONTAINED SETUP FOR COMBINED 4-PANEL FIGURE
# Rebuilds meta_bac and meta_euk from scratch with ALL required columns:
# richness, significance letters, and PCoA axes — regardless of what has
# or hasn't already run earlier in the session.
# =============================================================================

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

# --- Bacterial richness ---
meta_bac <- read.delim("Drymats_metadata_final.txt")
bac_asv  <- read.csv("ASV_table_rare_notufts.csv", check.names = FALSE, row.names = 1)
bac_asv  <- as.data.frame(t(bac_asv))
bac_asv  <- cbind(Sample = rownames(bac_asv), bac_asv)
rownames(bac_asv) <- NULL
meta_bac$Sample <- as.character(meta_bac$Sample)
bac_asv$Sample  <- as.character(bac_asv$Sample)
meta_bac <- meta_bac[meta_bac$Sample %in% bac_asv$Sample, ]
bac_asv  <- bac_asv[match(meta_bac$Sample, bac_asv$Sample), ]
row.names(bac_asv) <- meta_bac$Sample
bac_asv  <- bac_asv[, -1]
meta_bac$bac_Rich <- specnumber(bac_asv)

# --- Significance letters (from your Dunn test results) ---
cld_results <- data.frame(
  Environment = c("Desiccated", "Lake", "Pond"),
  Letter      = c("a", "b", "ab")
)
meta_bac <- merge(meta_bac, cld_results, by = "Environment", all.x = TRUE)
meta_bac$Environment <- factor(meta_bac$Environment, levels = c("Desiccated", "Pond", "Lake"))

# --- Bacterial PCoA axes (re-derive from the saved UniFrac matrix, avoiding
#     needing to rebuild the full phyloseq/tree object again) ---
varespec.uni_bac <- as.dist(read.csv("16S_UniFrac_Matrix_notufts.csv", row.names = 1, check.names = FALSE))
ordu_bac <- cmdscale(varespec.uni_bac, eig = TRUE, k = 2)
pcoa_df_bac <- as.data.frame(ordu_bac$points)
colnames(pcoa_df_bac) <- c("Axis.1", "Axis.2")
pcoa_df_bac$Sample <- rownames(pcoa_df_bac)
meta_bac <- merge(meta_bac, pcoa_df_bac, by.x = "Sample", by.y = "Sample")

# --- Repeat the same pattern for eukaryotes ---
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
meta_euk <- read.delim("Drymats_metadata_final.txt")
euk_asv  <- read.csv("ASV_table_rare_notufts.csv", check.names = FALSE, row.names = 1)  # adjust filename if different
euk_asv  <- as.data.frame(t(euk_asv))
euk_asv  <- cbind(Sample = rownames(euk_asv), euk_asv)
rownames(euk_asv) <- NULL
meta_euk$Sample <- as.character(meta_euk$Sample)
euk_asv$Sample  <- as.character(euk_asv$Sample)
meta_euk <- meta_euk[meta_euk$Sample %in% euk_asv$Sample, ]
euk_asv  <- euk_asv[match(meta_euk$Sample, euk_asv$Sample), ]
row.names(euk_asv) <- meta_euk$Sample
euk_asv  <- euk_asv[, -1]
meta_euk$euk_Rich <- specnumber(euk_asv)
meta_euk$Environment <- factor(meta_euk$Environment, levels = c("Desiccated", "Pond", "Lake"))

varespec.uni_euk <- as.dist(read.csv("18S_UniFrac_Matrix_notufts.csv", row.names = 1, check.names = FALSE))
ordu_euk <- cmdscale(varespec.uni_euk, eig = TRUE, k = 2)
pcoa_df_euk <- as.data.frame(ordu_euk$points)
colnames(pcoa_df_euk) <- c("Axis.1", "Axis.2")
pcoa_df_euk$Sample <- rownames(pcoa_df_euk)
meta_euk <- merge(meta_euk, pcoa_df_euk, by.x = "Sample", by.y = "Sample")

meta_bac$Letter <- meta_bac$Letter.x
meta_bac$Letter.x <- NULL
meta_bac$Letter.y <- NULL

library(patchwork)

env_colours <- c("Desiccated" = "#E31A1C",
                 "Pond"       = "#1F78B4",
                 "Lake"       = "#33A02C")

type_colours <- c(
  "Untersee Desiccated Mat"       = "#CC0000",
  "Snow Petrel Desiccated Mat"    = "#FF6600",
  "Avalanche pond"                = "#00008B",
  "Eastern lateral moraine pond"  = "#0057FF",
  "Southern pond"                 = "#00BFFF",
  "Western lateral moraine pond"  = "#87CEEB",
  "flat mat"                      = "#005000",
  "pinnacle"                      = "#00A550",
  "cone"                          = "#00FF7F"
)

env_shapes <- c("Desiccated" = 15,
                "Pond"       = 17,
                "Lake"       = 16)

cld_results <- data.frame(
  Environment = c("Desiccated", "Lake", "Pond"),
  Letter      = c("a", "b", "ab")
)
meta_bac <- merge(meta_bac, cld_results, by = "Environment", all.x = TRUE)
meta_bac$Environment <- factor(meta_bac$Environment, levels = c("Desiccated", "Pond", "Lake"))

cld_results_euk <- data.frame(
  Environment = c("Desiccated", "Lake", "Pond"),
  Letter      = c("a", "b", "a")
)
meta_euk <- merge(meta_euk, cld_results_euk, by = "Environment", all.x = TRUE)
meta_euk$Environment <- factor(meta_euk$Environment, levels = c("Desiccated", "Pond", "Lake"))
# -----------------------------------------------------------------------
# Panel A — Bacterial (16S) Alpha Diversity
# -----------------------------------------------------------------------
panel_A <- ggplot(meta_bac, aes(x = Environment, y = bac_Rich, colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 4, alpha = 0.5) +
  scale_colour_manual(values = env_colours) +
  labs(x = NULL, y = "ASV Richness", title = "16S rRNA gene Alpha Diversity") +
  theme_bw() +
  ylim(0, max(meta_bac$bac_Rich, na.rm = TRUE) * 1.15) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 14),
        axis.title.x  = element_text(face = "bold", size = 14, vjust = -2),
        axis.text.x   = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(face = "bold", size = 14, vjust = 4),
        plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_text(data = unique(meta_bac[, c("Environment", "Letter")]),
            aes(x = Environment, y = max(meta_bac$bac_Rich, na.rm = TRUE) * 1.05, label = Letter),
            color = "black", size = 5, fontface = "bold", vjust = 0)

# -----------------------------------------------------------------------
# Panel B — Eukaryotic (18S) Alpha Diversity
# NOTE: adjust y-limits/column names (euk_Rich, Letter) to match whatever
# object names the eukarya script actually produced — mirrored here from
# the bacterial script's structure.
# -----------------------------------------------------------------------
panel_B <- ggplot(meta_euk, aes(x = Environment, y = euk_Rich, colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 4, alpha = 0.5) +
  scale_colour_manual(values = env_colours) +
  labs(x = NULL, y = "ASV Richness", title = "18S rRNA gene Alpha Diversity") +
  theme_bw() +
  ylim(0, max(meta_euk$euk_Rich, na.rm = TRUE) * 1.25) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 14),
        axis.title.x  = element_text(face = "bold", size = 14, vjust = -2),
        axis.text.x   = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(face = "bold", size = 14, vjust = 4),
        plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_text(data = unique(meta_euk[, c("Environment", "Letter")]),
            aes(x = Environment, y = max(meta_euk$euk_Rich, na.rm = TRUE) * 1.15, label = Letter),
            color = "black", size = 5, fontface = "bold", vjust = 0)

# If euk alpha diversity also has significance letters (cld_results equivalent),
# add the same geom_text() layer as panel_A here, matching column names.

# -----------------------------------------------------------------------
# Ensure Type factor levels are set in the desired legend order:
# all Desiccated types, then all Pond types, then all Lake types
# -----------------------------------------------------------------------
type_order <- c(
  "Untersee Desiccated Mat",
  "Snow Petrel Desiccated Mat",
  "Avalanche pond",
  "Eastern lateral moraine pond",
  "Southern pond",
  "Western lateral moraine pond",
  "flat mat",
  "pinnacle",
  "cone"
)

meta_bac$Type <- factor(meta_bac$Type, levels = type_order)
meta_euk$Type <- factor(meta_euk$Type, levels = type_order)

# -----------------------------------------------------------------------
# Convex hulls grouped by Type (not Environment), matching the Type colour
# scheme used for points
# -----------------------------------------------------------------------
find_hulls <- function(df) df[chull(df$Axis.1, df$Axis.2), ]

hulls_bac_type <- meta_bac %>% group_by(Type) %>% do(find_hulls(.))
hulls_euk_type <- meta_euk %>% group_by(Type) %>% do(find_hulls(.))

# -----------------------------------------------------------------------
# Panel C — Bacterial (16S) Beta Diversity (PCoA), Type colour + hulls,
# Environment shape
# -----------------------------------------------------------------------
panel_C <- ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_polygon(data = hulls_bac_type, aes(colour = Type, fill = Type),
               alpha = 0.1, linewidth = 0.25, show.legend = FALSE) +
  geom_point(size = 3.5, aes(colour = Type, shape = Environment), alpha = 0.7) +
  scale_colour_manual(values = type_colours, breaks = type_order) +
  scale_fill_manual(values = type_colours, breaks = type_order) +
  scale_shape_manual(values = env_shapes) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%", title = "16S rRNA gene Beta Diversity") +
  theme_bw() +
  theme(legend.position = "none",   # legend shown once, on panel D
        plot.title    = element_text(size = 14),
        axis.title.x  = element_text(face = "bold", size = 14, vjust = 2),
        axis.title.y  = element_text(face = "bold", size = 14, vjust = 2),
        axis.text     = element_text(size = 11),
        plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# -----------------------------------------------------------------------
# Panel D — Eukaryotic (18S) Beta Diversity (PCoA), Type colour + hulls,
# Environment shape
# -----------------------------------------------------------------------
panel_D <- ggplot(meta_euk, aes(Axis.1, Axis.2)) +
  geom_polygon(data = hulls_euk_type, aes(colour = Type, fill = Type),
               alpha = 0.1, linewidth = 0.25, show.legend = FALSE) +
  geom_point(size = 3.5, aes(colour = Type, shape = Environment), alpha = 0.7) +
  scale_colour_manual(values = type_colours, breaks = type_order) +
  scale_fill_manual(values = type_colours, breaks = type_order) +
  scale_shape_manual(values = env_shapes) +
  labs(x = "PC1: 59.1%", y = "PC2: 46.79%", title = "18S rRNA gene Beta Diversity") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title    = element_text(size = 14),
        axis.title.x  = element_text(face = "bold", size = 14, vjust = 2),
        axis.title.y  = element_text(face = "bold", size = 14, vjust = 2),
        axis.text     = element_text(size = 11),
        legend.title  = element_text(face = "bold", size = 11),
        legend.text   = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        plot.margin   = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  guides(colour = guide_legend(override.aes = list(size = 3), order = 1),
         shape  = guide_legend(override.aes = list(size = 3), order = 2))

# -----------------------------------------------------------------------
# Combine into a single 2x2 figure with panel labels A-D and one shared legend
# -----------------------------------------------------------------------
combined_4panel <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/")

# 1. Widen the discrete x-axis spacing in panels A and B specifically,
#    so the boxplots aren't crushed together
panel_A <- panel_A + scale_x_discrete(expand = expansion(add = 0.6))
panel_B <- panel_B + scale_x_discrete(expand = expansion(add = 0.6))

# 2. Increase overall figure width to accommodate panel D's legend
#    without squeezing the plot areas
pdf("Combined_Alpha_Beta_Diversity_4panel_notufts.pdf", width = 14, height = 10)
print(combined_4panel)
dev.off()

ggsave("Combined_Alpha_Beta_Diversity_4panel_notufts.png", combined_4panel,
       width = 12, height = 10, dpi = 300)

cat("Saved: Combined_Alpha_Beta_Diversity_4panel_notufts.pdf/.png\n")

# =============================================================================
# FIGURE 2 (REVISED): Combined bacterial (A) and eukaryotic (B) phylum-level
# composition, faceted by habitat (Desiccated | Pond | Lake) with habitat
# labels shown as strip titles above each panel — matching the reference
# figure layout. Self-contained: regenerates all data needed from the rarefied,
# tuft-excluded RDS objects for both domains.
# =============================================================================

# =============================================================================
# FIGURE 2 (FINAL): Combined bacterial and eukaryotic phylum-level composition
# Two-step process: (1) generate aligned panels + legends as separate images,
# (2) composite with precise cm-based positioning using magick.
# Outputs: PDF (vector, for manuscript submission) and high-res JPEG.
# =============================================================================

library(mctoolsr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(magick)
library(scales)

# -----------------------------------------------------------------------
# Shared conventions
# -----------------------------------------------------------------------
taxa_palette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                  "#E31A1C", "#FDBF6F", "plum1", "orange", "lightskyblue",
                  "#FFFF99", "#B15928", "grey80")

habitat_order <- c("Desiccated", "Pond", "Lake")

custom_order <- c("Untersee Desiccated Mat 1", "Untersee Desiccated Mat 2",
                  "Untersee Desiccated Mat 3", "Untersee Desiccated Mat 4",
                  "Untersee Desiccated Mat 5", "Untersee Desiccated Mat 6",
                  "Untersee Desiccated Mat 7", "Snow Petrel Desiccated Mat 2",
                  "Snow Petrel Desiccated Mat 3", "Snow Petrel Desiccated Mat 4",
                  "Snow Petrel Desiccated Mat 5", "Snow Petrel Desiccated Mat 6",
                  "Eastern lateral moraine pond 1", "Eastern lateral moraine pond 2",
                  "Eastern lateral moraine pond 3", "Avalanche pond 1",
                  "Avalanche pond 2", "Southern pond 1", "Southern pond 3",
                  "Western lateral moraine pond 1", "Western lateral moraine pond 2",
                  "Western lateral moraine pond 3", "flat mat 1", "flat mat 2",
                  "flat mat 3", "pinnacle 1", "pinnacle 2", "pinnacle 3",
                  "pinnacle 4", "cone 1", "cone 2", "cone 3")

# -----------------------------------------------------------------------
# BACTERIA DATA
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
input_filt_rare_bac <- readRDS("bac_input_filt_rare_notufts.rds")

# Custom phylum-level column: phylum everywhere, except Proteobacteria
# split into its component classes (Alpha-, Beta-, Gammaproteobacteria, etc.)
input_filt_rare_bac$taxonomy_loaded$taxonomy2_custom <- as.character(input_filt_rare_bac$taxonomy_loaded$taxonomy2)
proteo_rows <- input_filt_rare_bac$taxonomy_loaded$taxonomy2 == "Proteobacteria"
input_filt_rare_bac$taxonomy_loaded$taxonomy2_custom[proteo_rows] <-
  as.character(input_filt_rare_bac$taxonomy_loaded$taxonomy3[proteo_rows])

cat("Custom phylum-level categories:\n")
print(sort(table(input_filt_rare_bac$taxonomy_loaded$taxonomy2_custom), decreasing = TRUE))

input_filt_rare_bac$taxonomy_loaded$taxonomy2 <- input_filt_rare_bac$taxonomy_loaded$taxonomy2_custom

phyla_bac <- summarize_taxonomy(input_filt_rare_bac, level = 2, report_higher_tax = FALSE)
input_filt_rare_bac$map_loaded <- cbind(SampleID = rownames(input_filt_rare_bac$map_loaded), input_filt_rare_bac$map_loaded)
pdat_bac <- plot_taxa_bars(phyla_bac, input_filt_rare_bac$map_loaded, "SampleID",
                           num_taxa = 12, data_only = TRUE)
pdat_bac$taxon <- as.factor(pdat_bac$taxon)
l_bac <- levels(pdat_bac$taxon)
pdat_bac$taxon <- factor(pdat_bac$taxon, levels = l_bac[1:13])

env_lookup_bac <- input_filt_rare_bac$map_loaded[, c("SampleID", "Environment")]
pdat_bac <- merge(pdat_bac, env_lookup_bac, by.x = "group_by", by.y = "SampleID")
pdat_bac$Environment <- factor(pdat_bac$Environment, levels = habitat_order)
pdat_bac$group_by <- factor(pdat_bac$group_by,
                            levels = custom_order[custom_order %in% pdat_bac$group_by])

pdat_bac$taxon_label <- as.character(pdat_bac$taxon)
pdat_bac$taxon_label[is.na(pdat_bac$taxon_label) | pdat_bac$taxon_label == "NA"] <- "Unclassified"

# -----------------------------------------------------------------------
# EUKARYA DATA
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
input_filt_rare_euk <- readRDS("euk_input_filt_rare_notufts.rds")

phyla_euk <- summarize_taxonomy(input_filt_rare_euk, level = 2, report_higher_tax = FALSE)
input_filt_rare_euk$map_loaded <- cbind(SampleID = rownames(input_filt_rare_euk$map_loaded),
                                        input_filt_rare_euk$map_loaded)

pdat_euk <- plot_taxa_bars(phyla_euk, input_filt_rare_euk$map_loaded, "SampleID",
                           num_taxa = 12, data_only = TRUE)
pdat_euk$taxon <- as.factor(pdat_euk$taxon)
l_euk <- levels(pdat_euk$taxon)
pdat_euk$taxon <- factor(pdat_euk$taxon, levels = l_euk[1:13])

env_lookup_euk <- input_filt_rare_euk$map_loaded[, c("SampleID", "Environment")]
pdat_euk <- merge(pdat_euk, env_lookup_euk, by.x = "group_by", by.y = "SampleID")
pdat_euk$Environment <- factor(pdat_euk$Environment, levels = habitat_order)
pdat_euk$group_by <- factor(pdat_euk$group_by,
                            levels = custom_order[custom_order %in% pdat_euk$group_by])

pdat_euk$taxon_label <- as.character(pdat_euk$taxon)
pdat_euk$taxon_label[is.na(pdat_euk$taxon_label) | pdat_euk$taxon_label == "NA"] <- "Unclassified"

# -----------------------------------------------------------------------
# BUILD PANEL A — Bacteria (legend title size increased: 9 -> 11)
# -----------------------------------------------------------------------
bacteria_palette <- c(
  "Actinobacteriota"    = "#A6CEE3",
  "Alphaproteobacteria" = "#1F78B4",
  "Armatimonadota"      = "#B2DF8A",
  "Bacteroidota"        = "#6A3D9A",
  "Betaproteobacteria"  = "#FB9A99",
  "Chloroflexi"         = "#E31A1C",
  "Cyanobacteria"       = "#33A02C",   # changed from peach to purple — was too close to Gammaproteobacteria's orange
  "Firmicutes"          = "#8DD3C7",   # changed from plum1 — avoids being close to the new Cyanobacteria purple
  "Gammaproteobacteria" = "orange",    # unchanged — now clearly distinct from Cyanobacteria
  "Gemmatimonadota"     = "#FDBF6F",   # changed from lightskyblue — was too close to Actinobacteriota's light blue
  "Other"               = "#FFFF99",
  "Planctomycetota"     = "#B15928",
  "Verrucomicrobiota"   = "grey80"
)

panel_A <- ggplot(pdat_bac, aes(x = group_by, y = mean_value, fill = taxon_label)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  facet_grid(~ Environment, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative Abundance", fill = "Bacteria") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = bacteria_palette) +
  theme_bw() +
  theme(axis.text.y   = element_text(size = 9, face = "bold"),
        axis.text.x   = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y  = element_text(face = "bold", size = 11),
        strip.text    = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey90"),
        legend.text   = element_text(size = 8, face = "bold"),
        legend.title  = element_text(face = "bold", size = 11),   # increased from 9
        legend.key.size = unit(0.4, "cm"),
        panel.spacing = unit(0.3, "lines"))

# -----------------------------------------------------------------------
# BUILD PANEL B — Eukarya (legend title size increased: 9 -> 11)
# -----------------------------------------------------------------------
eukarya_palette <- c(
  "Cercozoa"     = "#A6CEE3",
  "Chlorophyta"  = "#33A02C",
  "Ciliophora"   = "#B2DF8A",
  "Discoba"      = "#6A3D9A",
  "Discosea"     = "#FB9A99",
  "Fungi"        = "#E31A1C",
  "Gyrista"      = "#1F78B4",     # changed from peach — avoids clashing with Other's orange below
  "Metazoa"      = "plum1",
  "Other"        = "orange",
  "Rhodophyta"   = "#FDBF6F",
  "Streptophyta" = "#FFFF99",
  "Tubulinea"    = "#B15928",
  "Unclassified" = "grey80"
)

panel_B <- ggplot(pdat_euk, aes(x = group_by, y = mean_value, fill = taxon_label)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  facet_grid(~ Environment, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative Abundance", fill = "Eukarya") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = eukarya_palette) +
  theme_bw() +
  theme(axis.text.y   = element_text(size = 9, face = "bold"),
        axis.text.x   = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.title.y  = element_text(face = "bold", size = 11),
        strip.text    = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey90"),
        legend.text   = element_text(size = 8, face = "bold"),
        legend.title  = element_text(face = "bold", size = 11),   # increased from 9
        legend.key.size = unit(0.4, "cm"),
        panel.spacing = unit(0.3, "lines"))

# -----------------------------------------------------------------------
# STEP 1: Export panels (no legends) and legends separately
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

panel_A_noleg <- panel_A + theme(legend.position = "none")
panel_B_noleg <- panel_B + theme(legend.position = "none")

panels_combined <- panel_A_noleg / panel_B_noleg

ggsave("Figure2_panels_only.png", panels_combined,
       width = 13, height = 11, dpi = 300, bg = "white")

legend_A <- cowplot::get_legend(panel_A)
legend_B <- cowplot::get_legend(panel_B)

legend_A_plot <- ggdraw() + draw_grob(legend_A)
legend_B_plot <- ggdraw() + draw_grob(legend_B)

ggsave("Figure2_legend_Bacteria.png", legend_A_plot,
       width = 4, height = 5, dpi = 300, bg = "white")

ggsave("Figure2_legend_Eukarya.png", legend_B_plot,
       width = 4, height = 5, dpi = 300, bg = "white")

# -----------------------------------------------------------------------
# STEP 2: Composite — precise cm-based horizontal gap, vertical centering
# on each row's bar-plot area (confirmed working parameters)
# -----------------------------------------------------------------------
dpi <- 300
cm_to_px <- function(cm) round(cm * dpi / 2.54)

panels_img <- image_read("Figure2_panels_only.png")
legend_bac <- image_trim(image_read("Figure2_legend_Bacteria.png"))
legend_euk <- image_trim(image_read("Figure2_legend_Eukarya.png"))

panels_info     <- image_info(panels_img)
legend_bac_info <- image_info(legend_bac)
legend_euk_info <- image_info(legend_euk)

row_height_px <- panels_info$height / 2

plot_center_fraction_A <- 0.38
plot_center_fraction_B <- 0.38

x_offset <- panels_info$width + cm_to_px(1)

y_offset_A <- (row_height_px * plot_center_fraction_A) - (legend_bac_info$height / 2)
y_offset_B <- row_height_px + (row_height_px * plot_center_fraction_B) - (legend_euk_info$height / 2)

canvas_width  <- panels_info$width + max(legend_bac_info$width, legend_euk_info$width) + cm_to_px(2)
canvas_height <- panels_info$height

final_figure <- image_blank(width = canvas_width, height = canvas_height, color = "white")
final_figure <- image_composite(final_figure, panels_img, offset = "+0+0")
final_figure <- image_composite(final_figure, legend_bac,
                                offset = paste0("+", round(x_offset), "+", round(y_offset_A)))
final_figure <- image_composite(final_figure, legend_euk,
                                offset = paste0("+", round(x_offset), "+", round(y_offset_B)))

# -----------------------------------------------------------------------
# STEP 3: Save final outputs — high-res JPEG and PDF
# -----------------------------------------------------------------------
image_write(final_figure, path = "Figure2_Combined_Final.jpeg",
            format = "jpeg", quality = 100)

pdf_canvas_in_width  <- canvas_width  / dpi
pdf_canvas_in_height <- canvas_height / dpi

image_write(final_figure, path = "Figure2_Combined_Final.pdf", format = "pdf",
            density = dpi)

cat("Saved: Figure2_Combined_Final.jpeg (high-res) and Figure2_Combined_Final.pdf\n")
cat("Final canvas size:", round(pdf_canvas_in_width, 2), "x", round(pdf_canvas_in_height, 2), "inches\n")
