library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(Hmisc); packageVersion("Hmisc")
library(ggplot2); packageVersion("ggplot2")
library(plotly); packageVersion("plotly")
library(mctoolsr)
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/")

tax_table_fp = '~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/seqtab_wTax_mctoolsr.txt'
map_fp       = '~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/Drymats_metadata_final.txt'
input        = load_taxa_table(tax_table_fp, map_fp)

# Filter out Bacteria, chloroplast and mitochondria
input_filt <- filter_taxa_from_input(input, taxa_to_remove = c("Bacteria", "chloroplast", "mitochondria"))
# Filter out NAs at Kingdom level
input_filt <- filter_taxa_from_input(input_filt, at_spec_level = 1, taxa_to_remove = "NA")
# Filter out singletons
input_filt <- filter_taxa_from_input(input_filt, filter_thresh = 0.027)
sort(rowSums(input_filt$data_loaded))
sort(colSums(input_filt$data_loaded))
# Rarefy at 3220
set.seed(42)  # same seed as bacterial script — keep consistent across both
input_filt_rare <- single_rarefy(input = input_filt, depth = 3220)
sort(colSums(input_filt_rare$data_loaded))

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
write.csv(input_filt_rare$data_loaded, "ASV_table_rare.csv")
library(biomformat)
b <- make_biom(input_filt_rare$data_loaded)
write_biom(b, "ASV_table_rare.biom")
saveRDS(input_filt_rare, "euk_input_filt_rare.rds") # only do this once!
input_filt_rare <- readRDS("euk_input_filt_rare.rds")

################################################################################
# REMOVE TUFT SAMPLES
# Filament tufts (filament 1-5) are structurally atypical — loose cyanobacterial
# surface layers with no morphological equivalent in pond/desiccated dataset.
# NOTE: tufts are not eukaryotic outliers in 18S (unlike 16S where they are
# near-monoculture Cyanobacteriia), but are removed for consistency with the
# 16S analysis and because they lack structural equivalents in the other habitats.
# Removed following Reviewer 1 and Reviewer 2 recommendations.
################################################################################

# Verify Type levels before filtering
print(table(input_filt_rare$map_loaded$Type))

# Remove filament tuft samples
input_filt_rare <- filter_data(input_filt_rare,
                                'Type',
                                filter_vals = 'filament')

# Verify removal — should show 32 samples, no filament type
print(table(input_filt_rare$map_loaded$Type))
print(paste("Samples remaining:", ncol(input_filt_rare$data_loaded)))
sort(colSums(input_filt_rare$data_loaded))

# Save tuft-removed object
write.csv(input_filt_rare$data_loaded, "ASV_table_rare_notufts.csv")
library(biomformat)
b_nt <- make_biom(input_filt_rare$data_loaded)
write_biom(b_nt, "ASV_table_rare_notufts.biom")
saveRDS(input_filt_rare, "euk_input_filt_rare_notufts.rds")
message("Saved: euk_input_filt_rare_notufts.rds (filament tufts removed)")

################################################################################

# Libraries
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
library(RColorBrewer)
library(indicspecies)
theme_set(theme_bw())
show_col(hue_pal()(3))

# Custom sample order — filament samples removed
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

#PHYLA
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
taxonomy_loaded <- input_filt_rare$taxonomy_loaded
taxonomy_loaded <- taxonomy_loaded %>%
  mutate(taxonomy2 = case_when(taxonomy2 == "NA" ~ "Other", TRUE ~ taxonomy2))
input_filt_rare$taxonomy_loaded <- taxonomy_loaded
phyla <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = FALSE)
input_filt_rare$map_loaded <- cbind(SampleID = rownames(input_filt_rare$map_loaded), input_filt_rare$map_loaded)
pdat <- plot_taxa_bars(phyla, input_filt_rare$map_loaded, "SampleID", num_taxa = 12, data_only = TRUE)
pdat$taxon    <- as.factor(pdat$taxon)
l             <- levels(pdat$taxon)
pdat$taxon    <- factor(pdat$taxon, levels = c(l[1:13]))
pdat$group_by <- factor(pdat$group_by, levels = custom_order)

pdf("Phyla_notufts.pdf", width = 13, height = 6)
ggplot(data = pdat, aes(x = group_by, y = mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black") +
  labs(x = "Samples", y = "Relative Abundance", fill = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                               "#E31A1C","#FDBF6F", "plum1", "orange", "lightskyblue",
                               "#FFFF99", "#B15928","grey80")) +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 8, angle = 50, hjust = 1),
        axis.title   = element_text(face = "bold", size = 14))
dev.off()

#ORDER
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
taxonomy_loaded <- input_filt_rare$taxonomy_loaded
taxonomy_loaded <- taxonomy_loaded %>%
  mutate(taxonomy3 = case_when(taxonomy3 == "NA" ~ "Other", TRUE ~ taxonomy3))
input_filt_rare$taxonomy_loaded <- taxonomy_loaded
order_tax <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = FALSE)
input_filt_rare$map_loaded <- cbind(SampleID = rownames(input_filt_rare$map_loaded), input_filt_rare$map_loaded)
pdat <- plot_taxa_bars(order_tax, input_filt_rare$map_loaded, "SampleID", num_taxa = 12, data_only = TRUE)
pdat$taxon    <- as.factor(pdat$taxon)
l             <- levels(pdat$taxon)
pdat$taxon    <- factor(pdat$taxon, levels = c(l[1:13]))
pdat$group_by <- factor(pdat$group_by, levels = custom_order)

pdf("Order_notufts.pdf", width = 13, height = 6)
ggplot(data = pdat, aes(x = group_by, y = mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black") +
  labs(x = "Samples", y = "Relative Abundance", fill = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                               "#E31A1C","#FDBF6F", "plum1", "orange", "lightskyblue",
                               "#FFFF99", "#B15928","grey80")) +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 8, angle = 50, hjust = 1),
        axis.title   = element_text(face = "bold", size = 14))
dev.off()

#FAMILY
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
taxonomy_loaded <- input_filt_rare$taxonomy_loaded
taxonomy_loaded <- taxonomy_loaded %>%
  mutate(taxonomy4 = case_when(taxonomy4 == "NA" ~ "Other", TRUE ~ taxonomy4))
input_filt_rare$taxonomy_loaded <- taxonomy_loaded
family <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = FALSE)
input_filt_rare$map_loaded <- cbind(SampleID = rownames(input_filt_rare$map_loaded), input_filt_rare$map_loaded)
pdat <- plot_taxa_bars(family, input_filt_rare$map_loaded, "SampleID", num_taxa = 12, data_only = TRUE)
pdat$taxon    <- as.factor(pdat$taxon)
l             <- levels(pdat$taxon)
pdat$taxon    <- factor(pdat$taxon, levels = c(l[1:13]))
pdat$group_by <- factor(pdat$group_by, levels = custom_order)

pdf("Family_notufts.pdf", width = 13, height = 6)
ggplot(data = pdat, aes(x = group_by, y = mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black") +
  labs(x = "Samples", y = "Relative Abundance", fill = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                               "#E31A1C","#FDBF6F", "plum1", "orange", "lightskyblue",
                               "#FFFF99", "#B15928","grey80")) +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 8, angle = 50, hjust = 1),
        axis.title   = element_text(face = "bold", size = 14))
dev.off()

#GENUS
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
taxonomy_loaded <- input_filt_rare$taxonomy_loaded
taxonomy_loaded <- taxonomy_loaded %>%
  mutate(taxonomy6 = case_when(taxonomy6 == "NA" ~ "Other", TRUE ~ taxonomy6))
input_filt_rare$taxonomy_loaded <- taxonomy_loaded
genus <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = FALSE)
input_filt_rare$map_loaded <- cbind(SampleID = rownames(input_filt_rare$map_loaded), input_filt_rare$map_loaded)
pdat <- plot_taxa_bars(genus, input_filt_rare$map_loaded, "SampleID", num_taxa = 12, data_only = TRUE)
pdat$taxon    <- as.factor(pdat$taxon)
l             <- levels(pdat$taxon)
pdat$taxon    <- factor(pdat$taxon, levels = c(l[1:13]))
pdat$group_by <- factor(pdat$group_by, levels = custom_order)

pdf("Genus_notufts.pdf", width = 13, height = 6)
ggplot(data = pdat, aes(x = group_by, y = mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = "black") +
  labs(x = "Samples", y = "Relative Abundance", fill = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                               "#E31A1C","#FDBF6F", "plum1", "orange", "lightskyblue",
                               "#FFFF99", "#B15928","grey80")) +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 8, angle = 50, hjust = 1),
        axis.title   = element_text(face = "bold", size = 14))
dev.off()

# ===============================================
# SUMMARY STATISTICS FOR MANUSCRIPT TEXT — 18S
# Updated percentages after tuft removal
# ===============================================

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")

# Relative abundance matrix
rel_abund <- sweep(input_filt_rare$data_loaded, 2,
                   colSums(input_filt_rare$data_loaded), "/")

tax  <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded

# Function to get mean ± SD by environment for a given taxon — returns a data frame
get_env_stats_18s <- function(taxon_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == taxon_name]
  if (length(asvs) == 0) return(NULL)
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  agg <- data.frame(
    Taxon       = taxon_name,
    Tax_Level   = tax_level,
    Environment = names(tapply(rel_env, meta$Environment, mean)),
    Mean_Pct    = round(as.numeric(tapply(rel_env, meta$Environment, mean)), 1),
    SD_Pct      = round(as.numeric(tapply(rel_env, meta$Environment, sd)), 1)
  )
  agg$Mean_SD <- paste0(agg$Mean_Pct, " ± ", agg$SD_Pct, "%")
  agg
}

# Automatically get EVERY taxon at each level
get_all_level_stats_18s <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  do.call(rbind, lapply(taxa_at_level, get_env_stats_18s, tax_level = tax_level))
}

phylum_results_18s <- get_all_level_stats_18s("taxonomy2")
order_results_18s  <- get_all_level_stats_18s("taxonomy3")
family_results_18s <- get_all_level_stats_18s("taxonomy4")
genus_results_18s  <- get_all_level_stats_18s("taxonomy6")

all_results_18s <- rbind(phylum_results_18s, order_results_18s, family_results_18s, genus_results_18s)
write.csv(all_results_18s, "manuscript_taxa_summary_notufts_18S_ALL.csv", row.names = FALSE)

# ===============================================
# EUKARYOTIC TAXON-LEVEL STATISTICAL TESTS
# ===============================================
library(FSA)

# Reload notufts data
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')

rel_abund <- sweep(input_filt_rare$data_loaded, 2,
                   colSums(input_filt_rare$data_loaded), "/")
tax  <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded

test_taxon_18s <- function(taxon_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == taxon_name]
  if (length(asvs) == 0) return(NULL)
  
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  df      <- data.frame(abundance = rel_env, Environment = meta$Environment)
  kw      <- kruskal.test(abundance ~ Environment, data = df)
  
  out <- data.frame(Taxon = taxon_name, Tax_Level = tax_level,
                    KW_p = round(kw$p.value, 4), Comparison = NA, Dunn_p_adj = NA)
  
  if (kw$p.value < 0.05) {
    dunn <- dunnTest(abundance ~ Environment, data = df, method = "bh")
    pairwise <- data.frame(Taxon = taxon_name, Tax_Level = tax_level,
                           KW_p = round(kw$p.value, 4),
                           Comparison = dunn$res$Comparison,
                           Dunn_p_adj = round(dunn$res$P.adj, 4))
    out <- rbind(out, pairwise)
  }
  out
}

get_all_level_tests_18s <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  do.call(rbind, lapply(taxa_at_level, test_taxon_18s, tax_level = tax_level))
}


### Alpha Diversity — 18S
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
meta_euk <- read.delim("Drymats_metadata_final.txt")
euk_asv_raw <- read.csv("ASV_table_rare_notufts.csv", check.names = F, row.names = 1)
euk_asv <- as.data.frame(t(euk_asv_raw))
meta_euk <- meta_euk[meta_euk$Sample %in% rownames(euk_asv), ]
euk_asv  <- euk_asv[match(meta_euk$Sample, rownames(euk_asv)), ]
identical(rownames(euk_asv), meta_euk$Sample)
print(table(meta_euk$Environment))
row.names(euk_asv) <- meta_euk$Sample
meta_euk$Type      <- as.factor(meta_euk$Type)
meta_euk$euk_Rich  <- specnumber(euk_asv)

leveneTest(euk_Rich ~ Environment, data = meta_euk)
m1 <- aov(euk_Rich ~ Environment, data = meta_euk)
summary(m1)
library(FSA)
dunn_results <- dunnTest(euk_Rich ~ Environment, data = meta_euk, method = "bh")
print(dunn_results)

cld_results <- data.frame(
  Environment = c("Lake", "Pond", "Desiccated"),
  Letter      = c("a", "b", "b")
)
meta_euk <- merge(meta_euk, cld_results, by = "Environment", all.x = TRUE)
aggregate(meta_euk$euk_Rich, by = list(meta_euk$Environment), FUN = mean)

# Set factor levels: Desiccated, Pond, Lake left to right
meta_euk$Environment <- factor(meta_euk$Environment,
                               levels = c("Desiccated", "Pond", "Lake"))

# Colour palette — consistent with 16S
env_colours <- c("Desiccated" = "#E31A1C",
                 "Pond"       = "#1F78B4",
                 "Lake"       = "#33A02C")

pdf("Euk_Alpha_Diversity_notufts.pdf", width = 8, height = 6)
ggplot(meta_euk, aes(x = Environment, y = euk_Rich, colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 4, alpha = 0.5) +
  scale_colour_manual(values = env_colours) +
  labs(x = "Environment", y = "ASV Richness",
       title = "18S rRNA gene Alpha Diversity (tufts excluded)") +
  theme_bw() +
  ylim(15, 180) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 18),
        axis.title.x  = element_text(face = "bold", size = 18, vjust = -2),
        axis.text.x   = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 18, vjust = 4),
        plot.margin   = unit(c(0.7, 0.7, 0.7, 0.7), "cm")) +
  geom_text(data = unique(meta_euk[, c("Environment", "Letter")]),
            aes(x = Environment, y = 180, label = Letter),
            color = "black", size = 6, fontface = "bold")
dev.off()

# ===============================================
# TYPE-LEVEL ALPHA DIVERSITY — 18S
# ===============================================

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

meta_euk$Type <- factor(meta_euk$Type,
                        levels = c(
                          "Untersee Desiccated Mat",
                          "Snow Petrel Desiccated Mat",
                          "Avalanche pond",
                          "Eastern lateral moraine pond",
                          "Southern pond",
                          "Western lateral moraine pond",
                          "flat mat",
                          "pinnacle",
                          "cone"))

kw_type <- kruskal.test(euk_Rich ~ Type, data = meta_euk)
cat("Kruskal-Wallis across sample types: p =", round(kw_type$p.value, 4), "\n")
dunn_type <- dunnTest(euk_Rich ~ Type, data = meta_euk, method = "bh")
print(dunn_type$res[dunn_type$res$P.adj < 0.05, c("Comparison", "P.adj")])

pdf("Euk_Alpha_Diversity_Type_notufts.pdf", width = 10, height = 6)
ggplot(meta_euk, aes(x = Type, y = euk_Rich, colour = Type, shape = Environment)) +
  geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.7) +
  geom_point(size = 4, alpha = 0.7) +
  scale_colour_manual(values = type_colours) +
  scale_shape_manual(values = c("Desiccated" = 15,
                                "Pond"       = 17,
                                "Lake"       = 16)) +
  labs(x = "Sample Type", y = "ASV Richness",
       title = "18S rRNA gene Alpha Diversity by Sample Type (tufts excluded)") +
  theme_bw() +
  ylim(15, 180) +
  theme(legend.position    = "right",
        plot.title         = element_text(size = 14),
        axis.title.x       = element_text(face = "bold", size = 14, vjust = -2),
        axis.text.x        = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y        = element_text(size = 12),
        axis.title.y       = element_text(face = "bold", size = 14, vjust = 4),
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.45, "cm"),
        plot.margin        = unit(c(0.7, 0.7, 0.7, 0.7), "cm")) +
  geom_vline(xintercept = c(2.5, 6.5), linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  annotate("text", x = 1.5, y = 175, label = "Desiccated",
           size = 4, fontface = "bold", colour = "#E31A1C") +
  annotate("text", x = 4.5, y = 175, label = "Pond",
           size = 4, fontface = "bold", colour = "#1F78B4") +
  annotate("text", x = 8,   y = 175, label = "Lake",
           size = 4, fontface = "bold", colour = "#33A02C")
dev.off()

### Rarefaction Curves (18S)
# Generated to assess sequencing depth adequacy at rarefaction depth of 3220
# Supports methods text on 18S rarefaction depth

library(vegan)

# Use the filtered but NOT rarefied data for rarefaction curves
input_filt_norare <- filter_data(
  filter_taxa_from_input(
    filter_taxa_from_input(
      filter_taxa_from_input(input, taxa_to_remove = c("Bacteria", "chloroplast", "mitochondria")),
      at_spec_level = 1, taxa_to_remove = "NA"),
    filter_thresh = 0.027),
  'Type', filter_vals = 'filament')

# Transpose so samples are rows
asv_mat <- t(input_filt_norare$data_loaded)

# Generate rarefaction curve data
rarecurve_data <- rarecurve(asv_mat,
                             step    = 100,
                             sample  = 3220,
                             label   = FALSE,
                             tidy    = TRUE)

# Add metadata
meta_rare <- input_filt_norare$map_loaded
rarecurve_data$Environment <- meta_rare[rarecurve_data$Site, "Environment"]
rarecurve_data$Type        <- meta_rare[rarecurve_data$Site, "Type"]
rarecurve_data$Environment <- factor(rarecurve_data$Environment,
                                      levels = c("Lake", "Pond", "Desiccated"))

# Plot 1: Coloured by Environment
pdf("18S_Rarefaction_Curves_Environment_notufts.pdf", width = 8, height = 6)
ggplot(rarecurve_data, aes(x = Sample, y = Species,
                            group = Site, colour = Environment)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  geom_vline(xintercept = 3220, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  annotate("text", x = 3220, y = max(rarecurve_data$Species) * 0.95,
           label = "Rarefaction\ndepth (3220)",
           hjust = -0.05, size = 3.5, colour = "black") +
  labs(x      = "Number of sequences",
       y      = "ASV Richness",
       title  = "18S rRNA gene Rarefaction Curves (tufts excluded)",
       colour = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 14),
        axis.title      = element_text(face = "bold", size = 13),
        axis.text       = element_text(size = 11),
        legend.text     = element_text(size = 11))
dev.off()

type_colours <- c(
  "flat mat"                      = "#005000",
  "pinnacle"                      = "#00A550",
  "cone"                          = "#00FF7F",
  "Avalanche pond"                = "#00008B",
  "Eastern lateral moraine pond"  = "#0057FF",
  "Southern pond"                 = "#00BFFF",
  "Western lateral moraine pond"  = "#87CEEB",
  "Untersee Desiccated Mat"       = "#CC0000",
  "Snow Petrel Desiccated Mat"    = "#FF6600"
)

# Plot 2: Coloured by Type (finer resolution)
pdf("18S_Rarefaction_Curves_Type_notufts.pdf", width = 9, height = 6)
ggplot(rarecurve_data, aes(x = Sample, y = Species,
                            group = Site, colour = Type)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  geom_vline(xintercept = 3220, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  annotate("text", x = 3220, y = max(rarecurve_data$Species) * 0.95,
           label = "Rarefaction\ndepth (3220)",
           hjust = -0.05, size = 3.5, colour = "black") +
  scale_colour_manual(values = type_colours) +
  labs(x      = "Number of sequences",
       y      = "ASV Richness",
       title  = "18S rRNA gene Rarefaction Curves by sample type (tufts excluded)",
       colour = "Sample Type") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 13),
        axis.title      = element_text(face = "bold", size = 13),
        axis.text       = element_text(size = 11),
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
dev.off()

### Beta Diversity
install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
tree      <- read.tree("repset_aln.tre")
tax_table <- tax_table(as.matrix(input_filt_rare$taxonomy_loaded))
otu_table <- otu_table(input_filt_rare$data_loaded, taxa_are_rows = TRUE)
sam_table <- sample_data(input_filt_rare$map_loaded)
tot       <- phyloseq(otu_table, tax_table, sam_table, tree)
varespec.uni <- UniFrac(tot, weighted = TRUE)
write.csv(as.matrix(varespec.uni), "18S_UniFrac_Matrix_notufts.csv")

# Reload metadata aligned to distance matrix
meta_euk <- read.delim("Drymats_metadata_final.txt")
meta_euk <- meta_euk[meta_euk$Sample %in% labels(varespec.uni), ]
meta_euk <- meta_euk[match(labels(varespec.uni), meta_euk$Sample), ]
stopifnot(all(meta_euk$Sample == labels(varespec.uni)))

# PERMANOVA
m <- adonis2(varespec.uni ~ meta_euk$Environment, permutations = 999)
print(m) #0.001 ***
dispersion <- betadisper(varespec.uni, meta_euk$Environment)
anova(dispersion) #0.0008639 ***
permutest(dispersion, permutations = 999) #0.003 **

# Pairwise PERMANOVA
pairwise.adonis2 <- function(dist_matrix, groups, p.adjust.m = "bonferroni", permutations = 999) {
  co      <- combn(unique(as.character(groups)), 2)
  pairs   <- c(); F.Model <- c(); R2 <- c(); p.value <- c()
  for (i in 1:ncol(co)) {
    group1        <- co[1, i]; group2 <- co[2, i]
    subset_dist   <- as.dist(as.matrix(dist_matrix)[groups %in% c(group1, group2),
                                                     groups %in% c(group1, group2)])
    subset_groups <- groups[groups %in% c(group1, group2)]
    ad      <- adonis2(subset_dist ~ subset_groups, permutations = permutations)
    pairs   <- c(pairs,   paste(group1, "vs", group2))
    F.Model <- c(F.Model, ad$F[1])
    R2      <- c(R2,      ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  data.frame(Comparison = pairs, F.Model = F.Model, R2 = R2,
             p.value = p.value, p.adjusted = p.adjusted)
}
results <- pairwise.adonis2(varespec.uni, meta_euk$Environment)
print(results)

#          Comparison  F.Model         R2 p.value p.adjusted
#Desiccated vs Pond 1.957955 0.08916836   0.100      0.300
#Desiccated vs Lake 8.495080 0.29812446   0.001      0.003
#      Pond vs Lake 6.515445 0.26576900   0.001      0.003

# PERMANOVA with Year as covariate (sequential, Type I)
meta_euk$Year <- as.factor(ifelse(meta_euk$Environment == "Lake", 2011, 2019))
set.seed(123)
m_year_env <- adonis2(varespec.uni ~ Year + Environment,
                      data = meta_euk, permutations = 999, by = "terms")
print(m_year_env)
#Year         1  0.63324 0.20180 7.9561  0.001 ***
#Environment  1  0.19654 0.06263 2.4694  0.045 *  
set.seed(42)
m_marginal <- adonis2(varespec.uni ~ Year + Environment,
                      data = meta_euk, permutations = 999, by = "margin")
print(m_marginal)
#Environment  1  0.19654 0.06263 2.4694  0.045 *

# PERMANOVA at each taxonomic level — mirrors bacterial script
run_level_permanova_18s <- function(tax_level) {
  taxa_names <- tax[[tax_level]]
  agg <- aggregate(rel_abund, by = list(Taxon = taxa_names), FUN = sum)
  rownames(agg) <- agg$Taxon
  agg$Taxon <- NULL
  agg_t <- t(agg)
  
  bray <- vegdist(agg_t, method = "bray")
  ad <- adonis2(bray ~ meta_euk$Environment, permutations = 999)
  cat("\n", tax_level, "— PERMANOVA p =", ad$`Pr(>F)`[1],
      "R2 =", round(ad$R2[1], 4), "\n")
  return(ad)
}

phylum_permanova_18s <- run_level_permanova_18s("taxonomy2")
order_permanova_18s  <- run_level_permanova_18s("taxonomy3")
family_permanova_18s <- run_level_permanova_18s("taxonomy4")
genus_permanova_18s  <- run_level_permanova_18s("taxonomy6")

# PCoA
ordu <- ordinate(tot, "PCoA", "unifrac", weighted = TRUE)
df   <- as.data.frame(as.matrix(ordu$vectors))
df$sample <- row.names(df)
meta_euk$Axis.1 <- df$Axis.1
meta_euk$Axis.2 <- df$Axis.2
eig2 <- ordu$values$Eigenvalues
pct1 <- round(eig2[1] / sum(eig2) * 100, 2)
pct2 <- round(eig2[2] / sum(eig2) * 100, 2)

find_hull  <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
find_hulls <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
micro.hulls.euk <- meta_euk %>% group_by(Environment) %>% do(find_hulls(.))

pdf("Euk_Beta_Diversity_notufts.pdf", width = 7.78, height = 5)
ggplot(meta_euk, aes(Axis.1, Axis.2)) +
  geom_point(size = 4, aes(colour = Environment), alpha = 0.5) +
  labs(x = "PC1: 59.1%", y = "PC2: 46.79%",
       title = "18S rRNA gene Beta Diversity (tufts excluded)") +
  geom_polygon(data = micro.hulls.euk, aes(colour = Environment, fill = Environment),
               alpha = 0.1, linewidth = 0.25, show.legend = F) +
  theme_bw() +
  theme(legend.position = "right",
        plot.title    = element_text(size = 16),
        axis.title.x  = element_text(face = "bold", size = 16, vjust = 2),
        axis.text.x   = element_text(size = 14),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 16),
        plot.margin   = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# SIMPER
rel_euk  <- decostand(t(as.data.frame(as.matrix(tot@otu_table))), "total")
meta_euk_sim <- meta_euk[match(rownames(rel_euk), meta_euk$Sample), ]
stopifnot(all(rownames(rel_euk) == meta_euk_sim$Sample))
sim   <- with(meta_euk_sim, simper(rel_euk, Environment))
s_euk <- summary(sim)
head(s_euk$Desiccated_Pond, n = 10)
head(s_euk$Desiccated_Lake, n = 10)
head(s_euk$Pond_Lake,       n = 10)

# Top ESVs for vectors — update after checking SIMPER output
top_esvs_euk <- unique(c(
  rownames(head(s_euk$Desiccated_Pond, 5)),
  rownames(head(s_euk$Desiccated_Lake, 5)),
  rownames(head(s_euk$Pond_Lake,       5))
))
input_filt_rare$taxonomy_loaded[top_esvs_euk, ]

# SIMPER vectors — update ESV IDs and shortnames from output above
tot_transposed <- t(tot@otu_table)
w    <- wascores(x = ordu$vectors, w = tot_transposed)
wdf  <- as.data.frame(w)
wdf$Species <- rownames(w)

# Update these ESV IDs and shortnames based on your new SIMPER output
sub_euk <- subset(wdf, Species %in% c("ESV_9", "ESV_1", "ESV_12",
                                       "ESV_22", "ESV_20", "ESV_2", "ESV_4"))
sub_euk$shortnames <- c("Oxytrichidae", "Adineta vaga", "Adineta vaga", "Glaciozyma antarctica", "Sporobolomyces inositophilus", "Chlorococcum microstigmatum", "Leptophryidae")

micro.hulls.euk <- ddply(meta_euk, "Environment", find_hull)
meta_euk$Axis.1 <- df$Axis.1
meta_euk$Axis.2 <- df$Axis.2
micro.hulls.euk$Environment <- factor(micro.hulls.euk$Environment)
meta_euk$Environment        <- factor(meta_euk$Environment)
max_sample_distance  <- max(abs(meta_euk$Axis.1), abs(meta_euk$Axis.2))
max_species_distance <- max(abs(sub_euk$Axis.1), abs(sub_euk$Axis.2))
scaling_factor       <- 0.5 * (max_sample_distance / max_species_distance)
sub_euk$Axis.1       <- sub_euk$Axis.1 * scaling_factor
sub_euk$Axis.2       <- sub_euk$Axis.2 * scaling_factor
label_offset         <- 1.1; perp_offset <- 0.05
sub_euk$LabelX       <- sub_euk$Axis.1 * label_offset - perp_offset * sub_euk$Axis.2
sub_euk$LabelY       <- sub_euk$Axis.2 * label_offset + perp_offset * sub_euk$Axis.1

pdf("BetaDiversity_simper_scaled_notufts.pdf", width = 7.78, height = 5)
ggplot(meta_euk, aes(Axis.1, Axis.2)) +
  geom_point(size = 4, aes(colour = Environment), alpha = 0.5) +
  geom_polygon(data = micro.hulls.euk, aes(colour = Environment, fill = Environment),
               alpha = 0.1, linewidth = 0.25, show.legend = FALSE) +
  geom_segment(data = sub_euk, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey", linewidth = 0.6, inherit.aes = FALSE) +
  geom_text_repel(data = sub_euk, aes(x = LabelX, y = LabelY, label = shortnames),
                  size = 4, colour = "black", nudge_x = 0.02, nudge_y = 0.02,
                  segment.size = 0, segment.color = NA, inherit.aes = FALSE) +
  labs(x = "PC1: 59.1%", y = "PC2: 46.79%",
       title = "18S rRNA gene Beta Diversity (tufts excluded)") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title    = element_text(size = 16),
        axis.title.x  = element_text(face = "bold", size = 16, vjust = 2),
        axis.text.x   = element_text(size = 14),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 16),
        plot.margin   = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

### Type-level PCoA (addresses Reviewer 1 — finer resolution)
meta_euk$Type <- factor(meta_euk$Type,
                         levels = c("flat mat", "pinnacle", "cone",
                                    "Avalanche pond", "Eastern lateral moraine pond",
                                    "Southern pond", "Western lateral moraine pond",
                                    "Untersee Desiccated Mat", "Snow Petrel Desiccated Mat"))
meta_euk$Environment <- factor(meta_euk$Environment, levels = c("Lake", "Pond", "Desiccated"))

micro.hulls.type <- meta_euk %>% group_by(Type) %>% do(find_hulls(.))
micro.hulls.env  <- meta_euk %>% group_by(Environment) %>% do(find_hulls(.))

type_colours <- c(
  "flat mat"                      = "#005000",
  "pinnacle"                      = "#00A550",
  "cone"                          = "#00FF7F",
  "Avalanche pond"                = "#00008B",
  "Eastern lateral moraine pond"  = "#0057FF",
  "Southern pond"                 = "#00BFFF",
  "Western lateral moraine pond"  = "#87CEEB",
  "Untersee Desiccated Mat"       = "#CC0000",
  "Snow Petrel Desiccated Mat"    = "#FF6600"
)

# Figure 1: Type coloured, no environment hulls
pdf("Euk_BetaDiversity_Type_notufts.pdf", width = 9, height = 6)
ggplot(meta_euk, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.type,
               aes(colour = Type, fill = Type),
               alpha = 0.08, linewidth = 0.25, show.legend = FALSE) +
  geom_point(aes(colour = Type, shape = Environment), size = 4, alpha = 0.9) +
  scale_colour_manual(values = type_colours) +
  scale_fill_manual(values   = type_colours) +
  scale_shape_manual(values  = c("Lake" = 16, "Pond" = 17, "Desiccated" = 15)) +
  labs(x = "PC1: 59.1%", y = "PC2: 46.79%",
       title = "18S rRNA gene - Weighted UniFrac by sample type (tufts excluded)",
       colour = "Sample Type", shape = "Habitat") +
  theme_bw() +
  theme(legend.position  = "right",
        plot.title       = element_text(size = 13),
        axis.title.x     = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x      = element_text(size = 12),
        axis.text.y      = element_text(size = 12),
        axis.title.y     = element_text(face = "bold", size = 14),
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.45, "cm"),
        plot.margin      = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# Figure 2: Type colours with dashed environment hulls overlaid
pdf("Euk_BetaDiversity_Type_EnvHulls_notufts.pdf", width = 9, height = 6)
ggplot(meta_euk, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.env,
               aes(group = Environment, colour = Environment),
               fill = NA, linewidth = 0.6, linetype = "dashed", show.legend = FALSE) +
  geom_polygon(data = micro.hulls.type,
               aes(colour = Type, fill = Type),
               alpha = 0.1, linewidth = 0.2, show.legend = FALSE) +
  geom_point(aes(colour = Type, shape = Environment), size = 4, alpha = 0.9) +
  scale_colour_manual(values = type_colours) +
  scale_fill_manual(values   = type_colours) +
  scale_shape_manual(values  = c("Lake" = 16, "Pond" = 17, "Desiccated" = 15)) +
  labs(x = "PC1: 59.1%", y = "PC2: 46.79%",
       title = "18S rRNA gene - Sample types within habitat boundaries (tufts excluded)",
       colour = "Sample Type", shape = "Habitat") +
  theme_bw() +
  theme(legend.position  = "right",
        plot.title       = element_text(size = 12),
        axis.title.x     = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x      = element_text(size = 12),
        axis.text.y      = element_text(size = 12),
        axis.title.y     = element_text(face = "bold", size = 14),
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.45, "cm"),
        plot.margin      = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

### Network analysis 18S
library(ggplot2); library(igraph); library(scales)
library(mctoolsr); library(Hmisc); library(dplyr)
library(RColorBrewer)
set.seed(500)

input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
top18s    <- as.data.frame(head(sort(rowSums(input_filt_rare$data_loaded), decreasing = TRUE), n = 100))
write.csv(top18s, "top18s_notufts.csv")
top18s    <- read.csv("top18s_notufts.csv", row.names = 1)
top200euk <- filter_taxa_from_input(input_filt_rare, taxa_IDs_to_keep = rownames(top18s))
cor.cutoff <- 0.8

# Taxonomy colours — 18S has more phyla, use colorRampPalette
top200euk$taxonomy_loaded$taxonomy2 <- as.character(top200euk$taxonomy_loaded$taxonomy2)
top200euk$taxonomy_loaded$taxonomy2[top200euk$taxonomy_loaded$taxonomy2 == "NA"] <- "Unassigned"
top200euk$taxonomy_loaded$taxonomy2 <- as.factor(top200euk$taxonomy_loaded$taxonomy2)
l     <- levels(top200euk$taxonomy_loaded$taxonomy2)
colrs <- c(colorRampPalette(brewer.pal(12, "Paired"))(length(l) - 1), "grey60")

build_network_euk <- function(data, cor.cutoff) {
  cor.mat <- rcorr(t(data), type = "spearman")
  diag(cor.mat$r) <- 0
  cor.mat$r[is.na(cor.mat$r)] <- 0
  net <- graph_from_adjacency_matrix(cor.mat$r, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[abs(weight) < cor.cutoff])
  list(net = net, cor = cor.mat)
}

plot_network_euk <- function(net, cor.mat, deg, taxonomy, colrs, filename) {
  V(net)$phylum <- factor(taxonomy$taxonomy2[match(V(net)$name, rownames(taxonomy))],
                           levels = levels(taxonomy$taxonomy2))
  abs_w     <- abs(E(net)$weight)
  net_abs   <- net; E(net_abs)$weight <- abs_w
  comm      <- cluster_fast_greedy(net_abs)
  layout_cl <- layout_with_fr(net_abs)
  pdf(filename, width = 9, height = 4.8)
  par(mar = c(0, 0, 0, 0)); par(xpd = TRUE)
  plot(net,
       vertex.color       = colrs[as.numeric(V(net)$phylum)],
       vertex.size        = deg * 0.3,
       vertex.shape       = "circle",
       vertex.frame.color = "black",
       vertex.label       = NA,
       edge.color         = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved        = 0.2,
       edge.width         = abs_w * 0.5,
       layout             = layout_cl * 10)
  legend(x = -1.7, y = 0.8, levels(taxonomy$taxonomy2),
         pch = 21, col = NA, pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)
  dev.off()
}

#### 1. Lake ####
top200euk_lake <- filter_data(top200euk, 'Environment', keep_vals = 'Lake')
lake_data      <- top200euk$data_loaded[, rownames(top200euk_lake$map_loaded)]
write.csv(lake_data, "lake_data_notufts.csv")
lake_res  <- build_network_euk(lake_data, cor.cutoff)
net_lake  <- lake_res$net
deg_lake  <- degree(net_lake, mode = "all")
cat("Lake 18S network (no tufts): Edges =", length(E(net_lake)),
    "| Mean degree =", round(mean(deg_lake), 2),
    "| Transitivity =", round(transitivity(net_lake), 3),
    "| Positive =", sum(E(net_lake)$weight > 0),
    "| Negative =", sum(E(net_lake)$weight < 0), "\n")
plot_network_euk(net_lake, lake_res$cor, deg_lake,
                 top200euk$taxonomy_loaded, colrs,
                 "Lake_circle18S_notufts.pdf")

#### 2. Pond ####
top200euk_pond <- filter_data(top200euk, 'Environment', keep_vals = 'Pond')
pond_data      <- top200euk$data_loaded[, rownames(top200euk_pond$map_loaded)]
write.csv(pond_data, "pond_data_notufts.csv")
pond_res  <- build_network_euk(pond_data, cor.cutoff)
net_pond  <- pond_res$net
deg_pond  <- degree(net_pond, mode = "all")
cat("Pond 18S network: Edges =", length(E(net_pond)),
    "| Mean degree =", round(mean(deg_pond), 2),
    "| Transitivity =", round(transitivity(net_pond), 3),
    "| Positive =", sum(E(net_pond)$weight > 0),
    "| Negative =", sum(E(net_pond)$weight < 0), "\n")
plot_network_euk(net_pond, pond_res$cor, deg_pond,
                 top200euk$taxonomy_loaded, colrs,
                 "Pond_circle18S_notufts.pdf")

#### 3. Desiccated ####
top200euk_des <- filter_data(top200euk, 'Environment', keep_vals = 'Desiccated')
des_data      <- top200euk$data_loaded[, rownames(top200euk_des$map_loaded)]
write.csv(des_data, "des_data_notufts.csv")
des_res   <- build_network_euk(des_data, cor.cutoff)
net_des   <- des_res$net
deg_des   <- degree(net_des, mode = "all")
cat("Desiccated 18S network: Edges =", length(E(net_des)),
    "| Mean degree =", round(mean(deg_des), 2),
    "| Transitivity =", round(transitivity(net_des), 3),
    "| Positive =", sum(E(net_des)$weight > 0),
    "| Negative =", sum(E(net_des)$weight < 0), "\n")
plot_network_euk(net_des, des_res$cor, deg_des,
                 top200euk$taxonomy_loaded, colrs,
                 "Desiccated_circle18S_notufts.pdf")

### Indicator Species — All Eukaryota
input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
otu_table_ind <- input_filt_rare$data_loaded
metadata_ind  <- input_filt_rare$map_loaded
grouping_var  <- metadata_ind[colnames(otu_table_ind), "Environment"]
indicator     <- multipatt(t(otu_table_ind), grouping_var, control = how(nperm = 999))
summary(indicator)
sign_results     <- rownames_to_column(indicator$sign, var = "ESV")
sign_results     <- sign_results[!is.na(sign_results$p.value), ]
indicator_values <- rownames_to_column(as.data.frame(indicator$str), var = "ESV")
indicator_values <- indicator_values[, 1:(ncol(indicator_values) - 4)]
merged_ind       <- merge(sign_results, indicator_values, by = "ESV")
significant_results <- merged_ind[merged_ind$p.value < 0.05, ]
significant_results$Environment <- apply(
  significant_results[, c("Desiccated", "Pond", "Lake")], 1,
  function(row) if (max(row) >= 0.7) colnames(significant_results[, c("Desiccated","Pond","Lake")])[which.max(row)] else NA)
significant_results$Final <- apply(significant_results[, c("Pond","Desiccated","Lake")], 1, max)
significant_results_final <- significant_results[, c(1, which(colnames(significant_results) %in% c("Environment","Final")))]
taxonomy_ind <- rownames_to_column(input_filt_rare$taxonomy_loaded, var = "ESV")
get_best_taxonomy <- function(taxonomy_data) {
  taxonomy_data$best_taxonomy <- apply(taxonomy_data, 1, function(row) {
    for (level in c("taxonomy6","taxonomy5","taxonomy4","taxonomy3","taxonomy2","taxonomy1")) {
      if (!is.na(row[level]) && row[level] != "NA" && row[level] != "") return(row[level])
    }; return(NA)})
  taxonomy_data}
taxonomy_ind  <- get_best_taxonomy(taxonomy_ind)
merged_data   <- merge(significant_results_final, taxonomy_ind[, c("ESV","best_taxonomy")], by = "ESV", all.x = TRUE)
merged_data$Environment <- factor(merged_data$Environment, levels = c("Lake","Pond","Desiccated"))
ordered_data  <- merged_data[order(merged_data$Environment, merged_data$Final), ]
ordered_data$ESV_label <- factor(paste0(ordered_data$ESV," (",ordered_data$best_taxonomy,")"),
                                  levels = rev(unique(paste0(ordered_data$ESV," (",ordered_data$best_taxonomy,")"))))
ordered_data  <- ordered_data[ordered_data$Final >= 0.7, ]
# Remove unassigned Eukaryota entries
ordered_data  <- ordered_data %>%
  filter(!grepl("Eukaryota", best_taxonomy, ignore.case = TRUE))

pdf("Indicator_species_Eukarya_notufts.pdf", width = 8, height = 8)
ggplot(ordered_data, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.7, xend = Final - 0.015, y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("Lake" = "darkblue","Pond" = "darkgreen","Desiccated" = "brown")) +
  scale_x_continuous(limits = c(0.7, 1)) +
  labs(title = "Indicator Species by Environment — 18S (tufts excluded)",
       x = "Indicator Value", y = "ESV") +
  theme_bw() +
  theme(axis.text.y          = element_text(size = 8, face = "bold"),
        panel.grid.major.y   = element_line(color = "gray90"),
        panel.grid.minor     = element_blank(),
        legend.position      = "right")
dev.off()

### SourceTracker visualisation (18S)
# Re-run SourceTracker2 externally using ASV_table_rare_notufts.biom first:
# sourcetracker2 gibbs -i ASV_table_rare_notufts.biom \
#   -m Drymats_metadata_final.txt -o sourcetracker_notufts/ --jobs 4

### SourceTracker visualisation (18S)
# SourceTracker2 was run on NON-rarefied data (ASV_table_NONrare.biom) — correct
# for this algorithm which handles library size variation internally.
# Three separate analyses: each habitat designated as sink in turn.
# Results interpreted as COMPOSITIONAL OVERLAP not directional dispersal.
# Lake-involving comparisons carry additional cross-study caveat (Greco et al. 2020).

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_18S/")
meta_source <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_18S/Drymats_metadata_final.txt")

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/Sourcetracker_final/")
library(dplyr)

# ── Read mixing proportions from each sink subfolder ──────────────────────────
desic_sink <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_18S/sourcetracker_results_desiccated/mixing_proportions.txt",
                         check.names = FALSE)
lake_sink  <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_18S/sourcetracker_results_lake/mixing_proportions.txt",
                         check.names = FALSE)
pond_sink  <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_18S/sourcetracker_results_pond/mixing_proportions.txt",
                         check.names = FALSE)

# Add Environment from metadata to each sink
id_col <- colnames(desic_sink)[1]
desic_sink <- merge(desic_sink, meta_source, by.x = id_col, by.y = "SampleID", all.x = TRUE)
lake_sink  <- merge(lake_sink,  meta_source, by.x = id_col, by.y = "SampleID", all.x = TRUE)
pond_sink  <- merge(pond_sink,  meta_source, by.x = id_col, by.y = "SampleID", all.x = TRUE)

# Verify column names — source columns should be Lake, Pond, Desiccated, Unknown
print(head(desic_sink))
print(head(lake_sink))
print(head(pond_sink))

# ── Compute average source contributions per sink ─────────────────────────────
desic_avg <- colMeans(desic_sink[, c("Lake", "Pond", "Unknown")], na.rm = TRUE)
pond_avg  <- colMeans(pond_sink[,  c("Lake", "Desiccated", "Unknown")], na.rm = TRUE)
lake_avg  <- colMeans(lake_sink[,  c("Pond", "Desiccated", "Unknown")], na.rm = TRUE)

# ── Visualisation: three-panel pie chart ──────────────────────────────────────
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_18S/")

pdf("SourceTracker_all_sinks_18S_notufts.pdf", width = 14, height = 5)
par(mfrow = c(1, 3))

# Panel 1: Desiccated as sink (primary)
pie(c(desic_avg["Lake"], desic_avg["Pond"], desic_avg["Unknown"]),
    labels = c(paste0("Lake\n",    round(desic_avg["Lake"]    * 100, 1), "%"),
               paste0("Pond\n",    round(desic_avg["Pond"]    * 100, 1), "%"),
               paste0("Unknown\n", round(desic_avg["Unknown"] * 100, 1), "%")),
    col    = c("#33A02C", "#1F78B4", "grey80"),
    main   = "Desiccated as sink\n(primary analysis)",
    cex.main = 0.95)

# Panel 2: Pond as sink
pie(c(pond_avg["Lake"], pond_avg["Desiccated"], pond_avg["Unknown"]),
    labels = c(paste0("Lake\n",        round(pond_avg["Lake"]        * 100, 1), "%"),
               paste0("Desiccated\n",  round(pond_avg["Desiccated"]  * 100, 1), "%"),
               paste0("Unknown\n",     round(pond_avg["Unknown"]     * 100, 1), "%")),
    col    = c("#33A02C", "#E31A1C", "grey80"),
    main   = "Pond as sink",
    cex.main = 0.95)

# Panel 3: Lake as sink (exploratory)
pie(c(lake_avg["Pond"], lake_avg["Desiccated"], lake_avg["Unknown"]),
    labels = c(paste0("Pond\n",        round(lake_avg["Pond"]        * 100, 1), "%"),
               paste0("Desiccated*\n", round(lake_avg["Desiccated"]  * 100, 1), "%"),
               paste0("Unknown\n",     round(lake_avg["Unknown"]     * 100, 1), "%")),
    col    = c("#1F78B4", "#E31A1C", "grey80"),
    main   = "Lake as sink\n(exploratory — see caveats)",
    cex.main = 0.95)

mtext("* Desiccated source reflects shared ancestry rather than dispersal (circular)",
      side = 1, line = -1, outer = TRUE, cex = 0.75, adj = 0.98)

par(mfrow = c(1, 1))
dev.off()

# ── Print summary table for manuscript ────────────────────────────────────────
cat("\n=== SourceTracker: Average compositional overlap (18S, tufts excluded) ===\n")
cat("\nDesiccated as sink (primary analysis):\n")
print(round(desic_avg * 100, 1))
cat("\nPond as sink:\n")
print(round(pond_avg * 100, 1))
cat("\nLake as sink (exploratory):\n")
print(round(lake_avg * 100, 1))

# =============================================================================
# EUKARYOTIC (18S) TAXON-LEVEL STATISTICAL TESTS — ALL LEVELS
# Standalone, self-contained script: loads data fresh, defines all needed
# functions, and runs Kruskal-Wallis + Dunn tests at every taxonomic level
# (division/phylum, class, order, family, genus) — including family
# (taxonomy5), which was defined but never actually called in the original
# script, and so was missing from prior statistical output.
# =============================================================================

library(dplyr)
library(mctoolsr)
library(FSA)

# Guard against plyr masking dplyr::summarise/group_by, which can silently
# break grouped operations if plyr is loaded earlier in the same session
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

# -----------------------------------------------------------------------
# 1. LOAD DATA — rarefied, tuft-excluded eukaryotic dataset
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")

input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')

cat("Samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

rel_abund <- sweep(input_filt_rare$data_loaded, 2,
                   colSums(input_filt_rare$data_loaded), "/")
tax  <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded

# -----------------------------------------------------------------------
# 2. FUNCTIONS — per-taxon Kruskal-Wallis + Dunn test, and the wrapper
#    that runs it across every taxon at a given taxonomic level
# -----------------------------------------------------------------------
test_taxon_18s <- function(taxon_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == taxon_name]
  if (length(asvs) == 0) return(NULL)
  
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  df      <- data.frame(abundance = rel_env, Environment = meta$Environment)
  kw      <- kruskal.test(abundance ~ Environment, data = df)
  
  out <- data.frame(Taxon = taxon_name, Tax_Level = tax_level,
                    KW_p = round(kw$p.value, 4), Comparison = NA, Dunn_p_adj = NA)
  
  if (kw$p.value < 0.05) {
    dunn <- dunnTest(abundance ~ Environment, data = df, method = "bh")
    pairwise <- data.frame(Taxon = taxon_name, Tax_Level = tax_level,
                           KW_p = round(kw$p.value, 4),
                           Comparison = dunn$res$Comparison,
                           Dunn_p_adj = round(dunn$res$P.adj, 4))
    out <- rbind(out, pairwise)
  }
  out
}

get_all_level_tests_18s <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  do.call(rbind, lapply(taxa_at_level, test_taxon_18s, tax_level = tax_level))
}

# -----------------------------------------------------------------------
# 3. RUN TESTS AT EVERY TAXONOMIC LEVEL
#    (taxonomy1 = Supergroup is excluded, consistent with the rest of the
#    manuscript's analyses, which begin reporting at taxonomy2/Division)
# -----------------------------------------------------------------------
division_tests_18s <- get_all_level_tests_18s("taxonomy2")
class_tests_18s    <- get_all_level_tests_18s("taxonomy3")
order_tests_18s    <- get_all_level_tests_18s("taxonomy4")
family_tests_18s   <- get_all_level_tests_18s("taxonomy5")   # previously missing
genus_tests_18s    <- get_all_level_tests_18s("taxonomy6")

all_tests_18s <- rbind(division_tests_18s, class_tests_18s, order_tests_18s,
                       family_tests_18s, genus_tests_18s)

# -----------------------------------------------------------------------
# 4. SAVE
# -----------------------------------------------------------------------
write.csv(all_tests_18s, "manuscript_stats_ALL_taxa_notufts_18S_COMPLETE.csv", row.names = FALSE)

cat("\nSaved: manuscript_stats_ALL_taxa_notufts_18S_COMPLETE.csv\n")
cat("\nTaxa tested per level:\n")
for (lvl in c("taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5", "taxonomy6")) {
  n <- length(unique(all_tests_18s$Taxon[all_tests_18s$Tax_Level == lvl]))
  cat(" ", lvl, ":", n, "taxa\n")
}

# -----------------------------------------------------------------------
# 5. QUICK CHECK — any significant taxa specifically at family level
#    (taxonomy5), the level that was previously never tested
# -----------------------------------------------------------------------
family_sig <- family_tests_18s[is.na(family_tests_18s$Comparison), ]
family_sig <- family_sig[!is.na(family_sig$KW_p) & family_sig$KW_p < 0.05, ]

cat("\n=== Significant taxa at FAMILY level (KW p < 0.05) ===\n")
if (nrow(family_sig) > 0) {
  print(family_sig[order(family_sig$KW_p), c("Taxon", "KW_p")])
  
  cat("\nPairwise comparisons for these significant family-level taxa:\n")
  for (t in family_sig$Taxon) {
    cat("\n--", t, "--\n")
    print(all_tests_18s[all_tests_18s$Taxon == t & all_tests_18s$Tax_Level == "taxonomy5" &
                          !is.na(all_tests_18s$Comparison),
                        c("Comparison", "Dunn_p_adj")])
  }
} else {
  cat("None found.\n")
}
