library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(Hmisc); packageVersion("Hmisc")
library(ggplot2); packageVersion("ggplot2")
library(plotly); packageVersion("plotly")
library(mctoolsr)
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/")

#First: Change taxonomy of Burkholderiales from Gammaproteobacteria to Betaproteobacteria
data <- read.table("seqtab_wTax_mctoolsr_temp.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
data <- data %>%
  mutate(taxonomy = ifelse(
    grepl("Burkholderiales", taxonomy),
    gsub("Gammaproteobacteria", "Betaproteobacteria", taxonomy),
    taxonomy
  ))
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
write.table(data, "mod_seqtab_wTax_mctoolsr_temp.txt", sep = "\t", row.names = FALSE, quote = FALSE)

tax_table_fp = '~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/mod_seqtab_wTax_mctoolsr_temp.txt'
map_fp       = '~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/Drymats_metadata_final.txt'
input        = load_taxa_table(tax_table_fp, map_fp)

# Filter out chloroplast and mitochondria
input_filt <- filter_taxa_from_input(input, taxa_to_remove = c("Chloroplast","Mitochondria", "Eukaryota"))
# Filter out NAs at Kingdom level
input_filt <- filter_taxa_from_input(input_filt, at_spec_level = 1, taxa_to_remove = "NA")
# Filter out singletons
input_filt <- filter_taxa_from_input(input_filt, filter_thresh = 0.027)
sort(rowSums(input_filt$data_loaded))
sort(colSums(input_filt$data_loaded))
# Rarefy
set.seed(42)  # pick any integer — just document it and never change it
input_filt_rare <- single_rarefy(input = input_filt, depth = 20804)
sort(colSums(input_filt_rare$data_loaded))
write.csv(input_filt_rare$data_loaded, "ASV_table_rare.csv")
library(biomformat)
b <- make_biom(input_filt_rare$data_loaded)
write_biom(b, "ASV_table_rare.biom")
saveRDS(input_filt_rare, "bac_input_filt_rare.rds") # only do this once!
input_filt_rare <- readRDS("bac_input_filt_rare.rds")

################################################################################
# REMOVE TUFT SAMPLES
# Filament tufts (filament 1-5) are structurally atypical — loose cyanobacterial
# surface layers with no morphological equivalent in pond/desiccated dataset,
# near-monoculture Cyanobacteriia in 16S. Removed following Reviewer 1 and
# Reviewer 2 recommendations.
################################################################################

# Verify sample names before filtering
print(table(input_filt_rare$map_loaded$Type))

# Remove filament tufts — note lowercase "filament" to match metadata
input_filt_rare <- filter_data(input_filt_rare,
                                'Type',
                                filter_vals = 'filament')

# Verify removal — should show 32 samples, no filament type
print(table(input_filt_rare$map_loaded$Type))
print(paste("Samples remaining:", ncol(input_filt_rare$data_loaded)))
sort(colSums(input_filt_rare$data_loaded))

# Save tuft-removed object as new RDS
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
write.csv(input_filt_rare$data_loaded, "ASV_table_rare_notufts.csv")
saveRDS(input_filt_rare, "bac_input_filt_rare_notufts.rds")
message("Saved: bac_input_filt_rare_notufts.rds (filament tufts removed)")

################################################################################
# NORMALIZATION STRATEGY (hybrid approach)
# ─────────────────────────────────────────────────────────────────────────────
# Alpha diversity (ASV richness): rarefied counts at 20804 reads/sample
#   Rarefaction is appropriate here because richness is sensitive to library size
#
# Beta diversity (weighted UniFrac, PERMANOVA, PCoA): rarefied counts
#   UniFrac is a compositional metric; rarefied counts reduce library size effects
#
# Taxonomy bar plots: relative abundance via summarize_taxonomy() — NOT rarefied
#   Relative abundance is standard for compositional visualization
#
# SIMPER: relative abundance via decostand("total") — NOT rarefied
#   Relative abundance removes library size bias from pairwise dissimilarity
#
# Co-occurrence networks: rarefied counts
#   Correlation-based networks use count data; rarefied counts are standard
#
# Indicator species: rarefied counts
#   IndVal requires count data; rarefied counts ensure comparability
#
# For 16S the library size range is narrow (all samples ≥ 20804 after rarefaction)
# so rarefaction is well-justified. See rarefaction curves below.
################################################################################

### Rarefaction Curves (16S)
# Generated from filtered but non-rarefied data to assess depth adequacy
# Uses tuft-removed filtered data — must be run after the filtering section above

library(vegan)

# Rebuild filtered non-rarefied object with tufts removed for rarefaction curves
input_filt_norare <- filter_data(input_filt, 'Type', filter_vals = 'filament')

# Transpose so samples are rows
asv_mat_bac <- t(input_filt_norare$data_loaded)

# Generate rarefaction curve data
rarecurve_data_bac <- rarecurve(asv_mat_bac,
                                 step    = 500,
                                 sample  = 20804,
                                 label   = FALSE,
                                 tidy    = TRUE)

# Add metadata
meta_rare_bac <- input_filt_norare$map_loaded
rarecurve_data_bac$Environment <- meta_rare_bac[rarecurve_data_bac$Site, "Environment"]
rarecurve_data_bac$Type        <- meta_rare_bac[rarecurve_data_bac$Site, "Type"]
rarecurve_data_bac$Environment <- factor(rarecurve_data_bac$Environment,
                                          levels = c("Lake", "Pond", "Desiccated"))

# Type colour palette — defined early so available throughout script
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

# Plot 1: Coloured by Environment
pdf("16S_Rarefaction_Curves_Environment_notufts.pdf", width = 8, height = 6)
ggplot(rarecurve_data_bac, aes(x = Sample, y = Species,
                                group = Site, colour = Environment)) +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  geom_vline(xintercept = 20804, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  annotate("text", x = 20804, y = max(rarecurve_data_bac$Species) * 0.95,
           label = "Rarefaction\ndepth (20,804)",
           hjust = -0.05, size = 3.5, colour = "black") +
  labs(x      = "Number of sequences",
       y      = "ASV Richness",
       title  = "16S rRNA gene Rarefaction Curves (tufts excluded)",
       colour = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 14),
        axis.title      = element_text(face = "bold", size = 13),
        axis.text       = element_text(size = 11),
        legend.text     = element_text(size = 11))
dev.off()

# Plot 2: Coloured by Type
pdf("16S_Rarefaction_Curves_Type_notufts.pdf", width = 9, height = 6)
ggplot(rarecurve_data_bac, aes(x = Sample, y = Species,
                                group = Site, colour = Type)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  geom_vline(xintercept = 20804, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  annotate("text", x = 20804, y = max(rarecurve_data_bac$Species) * 0.95,
           label = "Rarefaction\ndepth (20,804)",
           hjust = -0.05, size = 3.5, colour = "black") +
  scale_colour_manual(values = type_colours) +
  labs(x      = "Number of sequences",
       y      = "ASV Richness",
       title  = "16S rRNA gene Rarefaction Curves by sample type (tufts excluded)",
       colour = "Sample Type") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 13),
        axis.title      = element_text(face = "bold", size = 13),
        axis.text       = element_text(size = 11),
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
dev.off()

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

#Stackedbar plot, Single bar per samples
#PHYLA
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")

# Custom phylum-level column: phylum everywhere, except Proteobacteria
# split into its component classes (Alpha-, Beta-, Gammaproteobacteria, etc.)
input_filt_rare$taxonomy_loaded$taxonomy2_custom <- as.character(input_filt_rare$taxonomy_loaded$taxonomy2)
proteo_rows <- input_filt_rare$taxonomy_loaded$taxonomy2 == "Proteobacteria"
input_filt_rare$taxonomy_loaded$taxonomy2_custom[proteo_rows] <-
  as.character(input_filt_rare$taxonomy_loaded$taxonomy3[proteo_rows])

cat("Custom phylum-level categories:\n")
print(sort(table(input_filt_rare$taxonomy_loaded$taxonomy2_custom), decreasing = TRUE))

input_filt_rare$taxonomy_loaded$taxonomy2 <- input_filt_rare$taxonomy_loaded$taxonomy2_custom

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
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
order_tax <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = FALSE)
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
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
taxonomy_loaded <- input_filt_rare$taxonomy_loaded
taxonomy_loaded <- taxonomy_loaded %>%
  mutate(taxonomy5 = case_when(
    taxonomy5 == "NA"             ~ "Unknown",
    taxonomy5 == "Unknown_Family" ~ "Unknown",
    TRUE                          ~ taxonomy5))
input_filt_rare$taxonomy_loaded <- taxonomy_loaded
family <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = FALSE)
input_filt_rare$map_loaded <- cbind(SampleID = rownames(input_filt_rare$map_loaded), input_filt_rare$map_loaded)
pdat <- plot_taxa_bars(family, input_filt_rare$map_loaded, "SampleID", num_taxa = 12, data_only = TRUE)
pdat$taxon    <- as.factor(pdat$taxon)
l             <- levels(pdat$taxon)
pdat$taxon    <- factor(pdat$taxon, levels = c(l[1:10], l[12], l[11], l[13]))
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
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
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
# SUMMARY STATISTICS FOR MANUSCRIPT TEXT
# Updated percentages after tuft removal
# ===============================================

input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")

# Relative abundance matrix (samples as columns)
rel_abund <- sweep(input_filt_rare$data_loaded, 2,
                   colSums(input_filt_rare$data_loaded), "/")

tax  <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded

# Function to get mean ± SD by environment for a given taxonomic class
get_env_stats <- function(class_name, tax_level = "taxonomy3") {
  asvs    <- rownames(tax)[tax[[tax_level]] == class_name]
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  result  <- tapply(rel_env, meta$Environment, function(x)
    paste0(round(mean(x), 1), " ± ", round(sd(x), 1), "%"))
  cat("\n", class_name, ":\n")
  print(result)
}

get_env_stats("Cyanobacteriia")
get_env_stats("Verrucomicrobiae")
get_env_stats("Bacteroidia")
get_env_stats("Actinobacteria")
get_env_stats("Planctomycetes")
get_env_stats("Gammaproteobacteria")
get_env_stats("Alphaproteobacteria")
get_env_stats("Chloroflexia")
get_env_stats("Betaproteobacteria")

# Function to get mean ± SD by environment for a given taxon
get_env_stats <- function(class_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == class_name]
  if (length(asvs) == 0) return(NULL)
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  agg <- data.frame(
    Taxon       = class_name,
    Tax_Level   = tax_level,
    Environment = names(tapply(rel_env, meta$Environment, mean)),
    Mean_Pct    = round(as.numeric(tapply(rel_env, meta$Environment, mean)), 1),
    SD_Pct      = round(as.numeric(tapply(rel_env, meta$Environment, sd)), 1)
  )
  agg$Mean_SD <- paste0(agg$Mean_Pct, " ± ", agg$SD_Pct, "%")
  agg
}

# Automatically get EVERY taxon at each level, rather than a manual list
get_all_level_stats <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "")]
  do.call(rbind, lapply(taxa_at_level, get_env_stats, tax_level = tax_level))
}

phylum_results <- get_all_level_stats("taxonomy2")  # phylum
class_results <- get_all_level_stats("taxonomy3")  # class
order_results  <- get_all_level_stats("taxonomy4")  # order
family_results <- get_all_level_stats("taxonomy5")  # family
genus_results  <- get_all_level_stats("taxonomy6")  # genus

all_results <- rbind(phylum_results, order_results, family_results, genus_results)
write.csv(all_results, "manuscript_taxa_summary_notufts_ALL.csv", row.names = FALSE)

# ===============================================
# MULTI-LEVEL STATISTICAL TESTS
# ===============================================
library(FSA)

test_taxon <- function(taxon_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == taxon_name]
  if(length(asvs) == 0) {
    cat("\n", taxon_name, ": NOT FOUND at", tax_level, "\n")
    return(NULL)
  }
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  df      <- data.frame(abundance   = rel_env,
                        Environment = meta$Environment)
  kw      <- kruskal.test(abundance ~ Environment, data = df)
  cat("\n", taxon_name, "(", tax_level, ") — Kruskal-Wallis p =",
      round(kw$p.value, 4), "\n")
  if(kw$p.value < 0.05) {
    dunn <- dunnTest(abundance ~ Environment, data = df, method = "bh")
    print(dunn$res[, c("Comparison", "P.adj")])
  }
}

# PHYLUM / CLASS level (taxonomy2/3)
cat("\n===== PHYLUM/CLASS LEVEL =====\n")
test_taxon("Cyanobacteriia",      "taxonomy3")
test_taxon("Bacteroidia",         "taxonomy3")
test_taxon("Chloroflexia",        "taxonomy3")
test_taxon("Actinobacteria",      "taxonomy3")
test_taxon("Gammaproteobacteria", "taxonomy3")
test_taxon("Verrucomicrobiae",    "taxonomy3")
test_taxon("Planctomycetes",      "taxonomy3")
test_taxon("Alphaproteobacteria", "taxonomy3")
test_taxon("Chlamydiae",          "taxonomy2")

# ORDER level (taxonomy4)
cat("\n===== ORDER LEVEL =====\n")
test_taxon("Cyanobacteriales",            "taxonomy4")
test_taxon("Oxyphotobacteria_Incertae_Sedis", "taxonomy4")
test_taxon("RD011",                       "taxonomy4")
test_taxon("Chloroflexales",              "taxonomy4")
test_taxon("Micrococcales",               "taxonomy4")
test_taxon("Flavobacteriales",            "taxonomy4")
test_taxon("Leptolyngbyales",             "taxonomy4")

# FAMILY level (taxonomy5)
cat("\n===== FAMILY LEVEL =====\n")
test_taxon("Phormidiaceae",       "taxonomy5")
test_taxon("Chloroflexaceae",     "taxonomy5")
test_taxon("Microbacteriaceae",   "taxonomy5")
test_taxon("Flavobacteriaceae",   "taxonomy5")
test_taxon("Leptolyngbyaceae",    "taxonomy5")
test_taxon("Caulobacteraceae",    "taxonomy5")

# GENUS level (taxonomy6)
cat("\n===== GENUS LEVEL =====\n")
test_taxon("Tychonema_CCAP_1459-11B", "taxonomy6")
test_taxon("Leptolyngbya_ANT.L67.1",  "taxonomy6")
test_taxon("Chloronema",              "taxonomy6")
test_taxon("Flavobacterium",          "taxonomy6")
test_taxon("Brevundimonas",           "taxonomy6")
test_taxon("Cryobacterium",           "taxonomy6")

# ===============================================
# COMPREHENSIVE STATISTICAL TEST — ALL TAXA, ALL LEVELS
# ===============================================
library(FSA)

# Modified test_taxon() that RETURNS a data frame row instead of just printing
test_taxon_df <- function(taxon_name, tax_level) {
  asvs <- rownames(tax)[tax[[tax_level]] == taxon_name]
  if (length(asvs) == 0) return(NULL)
  
  rel_env <- colSums(rel_abund[asvs, , drop = FALSE]) * 100
  df      <- data.frame(abundance = rel_env, Environment = meta$Environment)
  kw      <- kruskal.test(abundance ~ Environment, data = df)
  
  # Base row with overall KW result
  out <- data.frame(
    Taxon      = taxon_name,
    Tax_Level  = tax_level,
    KW_p       = round(kw$p.value, 4),
    Comparison = NA,
    Dunn_p_adj = NA
  )
  
  # If significant, add pairwise Dunn results as additional rows
  if (kw$p.value < 0.05) {
    dunn <- dunnTest(abundance ~ Environment, data = df, method = "bh")
    pairwise <- data.frame(
      Taxon      = taxon_name,
      Tax_Level  = tax_level,
      KW_p       = round(kw$p.value, 4),
      Comparison = dunn$res$Comparison,
      Dunn_p_adj = round(dunn$res$P.adj, 4)
    )
    out <- rbind(out, pairwise)
  }
  out
}

# Automatically pull every unique taxon at each level, excluding placeholders
get_all_level_tests <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  results_list  <- lapply(taxa_at_level, test_taxon_df, tax_level = tax_level)
  do.call(rbind, results_list)
}

# Run across all four taxonomic levels
phylum_tests <- get_all_level_tests("taxonomy2")   # phylum
class_tests  <- get_all_level_tests("taxonomy3")   # class
order_tests  <- get_all_level_tests("taxonomy4")   # order
family_tests <- get_all_level_tests("taxonomy5")   # family
genus_tests  <- get_all_level_tests("taxonomy6")   # genus

all_tests <- rbind(phylum_tests, order_tests, family_tests, genus_tests)
write.csv(all_tests, "manuscript_stats_ALL_taxa_notufts.csv", row.names = FALSE)

cat("\nTotal taxa tested:", length(unique(all_tests$Taxon)), "\n")
cat("Significant at phylum/class level (KW p<0.05):",
    sum(phylum_tests$KW_p < 0.05 & is.na(phylum_tests$Comparison)), "\n")


### Alpha Diversity
# Uses rarefied counts (20804 reads/sample) — appropriate for richness metrics
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
meta_bac <- read.delim("Drymats_metadata_final.txt")
bac_asv <- read.csv("ASV_table_rare_notufts.csv", check.names = F, row.names = 1)
bac_asv <- as.data.frame(t(bac_asv))
bac_asv <- cbind(Sample = rownames(bac_asv), bac_asv)
rownames(bac_asv) <- NULL
meta_bac$Sample <- as.character(meta_bac$Sample)
bac_asv$Sample  <- as.character(bac_asv$Sample)
# Keep only samples present after tuft removal
meta_bac <- meta_bac[meta_bac$Sample %in% bac_asv$Sample, ]
bac_asv  <- bac_asv[match(meta_bac$Sample, bac_asv$Sample), ]
identical(bac_asv$Sample, meta_bac$Sample)
row.names(bac_asv) <- meta_bac$Sample
bac_asv  <- bac_asv[,-1]
meta_bac$Type <- as.factor(meta_bac$Type)
meta_bac$bac_Rich <- specnumber(bac_asv)

# Environment — non-homogeneous variance, use Kruskal-Wallis + Dunn
leveneTest(bac_Rich ~ Environment, data = meta_bac) #0.1355
m1 <- kruskal.test(bac_Rich ~ Environment, data = meta_bac)
print(m1) #4.798e-05

### Alpha Diversity
# Uses rarefied counts (20804 reads/sample) — appropriate for richness metrics
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
meta_bac <- read.delim("Drymats_metadata_final.txt")
bac_asv <- read.csv("ASV_table_rare_notufts.csv", check.names = F, row.names = 1)
bac_asv <- as.data.frame(t(bac_asv))
bac_asv <- cbind(Sample = rownames(bac_asv), bac_asv)
rownames(bac_asv) <- NULL
meta_bac$Sample <- as.character(meta_bac$Sample)
bac_asv$Sample  <- as.character(bac_asv$Sample)
meta_bac <- meta_bac[meta_bac$Sample %in% bac_asv$Sample, ]
bac_asv  <- bac_asv[match(meta_bac$Sample, bac_asv$Sample), ]
identical(bac_asv$Sample, meta_bac$Sample)
row.names(bac_asv) <- meta_bac$Sample
bac_asv  <- bac_asv[,-1]
meta_bac$Type    <- as.factor(meta_bac$Type)
meta_bac$bac_Rich <- specnumber(bac_asv)

leveneTest(bac_Rich ~ Environment, data = meta_bac)
m1 <- kruskal.test(bac_Rich ~ Environment, data = meta_bac)
print(m1)
library(FSA)
library(multcompView)
dunn_results <- dunnTest(bac_Rich ~ Environment, data = meta_bac, method = "bh")
p_values     <- dunn_results$res[, c("Comparison", "P.adj")]
colnames(p_values) <- c("Comparison", "p.value")

cld_results <- data.frame(
  Environment = c("Desiccated", "Lake", "Pond"),
  Letter      = c("a", "b", "ab")
)
meta_bac <- merge(meta_bac, cld_results, by = "Environment", all.x = TRUE)

# Set factor levels for correct left-to-right order
meta_bac$Environment <- factor(meta_bac$Environment,
                               levels = c("Desiccated", "Pond", "Lake"))

# Colour palette
env_colours <- c("Desiccated" = "#E31A1C",
                 "Pond"       = "#1F78B4",
                 "Lake"       = "#33A02C")

pdf("Bac_Alpha_Diversity_notufts.pdf", width = 8, height = 6)
ggplot(meta_bac, aes(x = Environment, y = bac_Rich, colour = Environment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 4, alpha = 0.5) +
  scale_colour_manual(values = env_colours) +
  labs(x = "Environment", y = "ASV Richness",
       title = "16S rRNA gene Alpha Diversity") +
  theme_bw() +
  ylim(100, 500) +
  theme(legend.position = "none",
        plot.title    = element_text(size = 18),
        axis.title.x  = element_text(face = "bold", size = 18, vjust = -2),
        axis.text.x   = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 18, vjust = 4),
        plot.margin   = unit(c(0.7, 0.7, 0.7, 0.7), "cm")) +
  geom_text(data = unique(meta_bac[, c("Environment", "Letter")]),
            aes(x = Environment, y = 500, label = Letter),
            color = "black", size = 6, fontface = "bold")
dev.off()

# ===============================================
# TYPE-LEVEL ALPHA DIVERSITY
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

meta_bac$Type <- factor(meta_bac$Type,
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

kw_type <- kruskal.test(bac_Rich ~ Type, data = meta_bac)
cat("Kruskal-Wallis across sample types: p =", round(kw_type$p.value, 4), "\n")
dunn_type <- dunnTest(bac_Rich ~ Type, data = meta_bac, method = "bh")
print(dunn_type$res[dunn_type$res$P.adj < 0.05, c("Comparison", "P.adj")])

pdf("Bac_Alpha_Diversity_Type_notufts.pdf", width = 10, height = 6)
ggplot(meta_bac, aes(x = Type, y = bac_Rich, colour = Type, shape = Environment)) +
  geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.7) +
  geom_point(size = 4, alpha = 0.7) +
  scale_colour_manual(values = type_colours) +
  scale_shape_manual(values = c("Desiccated" = 15,
                                "Pond"       = 17,
                                "Lake"       = 16)) +
  labs(x = "Sample Type", y = "ASV Richness",
       title = "16S rRNA gene Alpha Diversity by Sample Type") +
  theme_bw() +
  ylim(100, 500) +
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
  annotate("text", x = 1.5, y = 490, label = "Desiccated",
           size = 4, fontface = "bold", colour = "#E31A1C") +
  annotate("text", x = 4.5, y = 490, label = "Pond",
           size = 4, fontface = "bold", colour = "#1F78B4") +
  annotate("text", x = 8,   y = 490, label = "Lake",
           size = 4, fontface = "bold", colour = "#33A02C")
dev.off()

### Beta Diversity
# Uses rarefied counts — weighted UniFrac is compositional; rarefied counts
# reduce library size effects on phylogenetic distance calculation
install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
tree      <- read.tree("repset_aln.tre")
tax_table <- tax_table(as.matrix(input_filt_rare$taxonomy_loaded))
otu_table <- otu_table(input_filt_rare$data_loaded, taxa_are_rows = TRUE)
sam_table <- sample_data(input_filt_rare$map_loaded)
tot       <- phyloseq(otu_table, tax_table, sam_table, tree)
varespec.uni <- UniFrac(tot, weighted = TRUE)
write.csv(as.matrix(varespec.uni), "16S_UniFrac_Matrix_notufts.csv")

# Reload metadata aligned to distance matrix
meta_bac <- read.delim("Drymats_metadata_final.txt")
meta_bac <- meta_bac[meta_bac$Sample %in% labels(varespec.uni), ]
meta_bac <- meta_bac[match(labels(varespec.uni), meta_bac$Sample), ]
stopifnot(all(meta_bac$Sample == labels(varespec.uni)))

# PERMANOVA
m <- adonis2(varespec.uni ~ meta_bac$Environment, permutations = 999) #0.001 ***
print(m)
dispersion <- betadisper(varespec.uni, meta_bac$Environment)
anova(dispersion) #0.2665
permutest(dispersion, permutations = 999) #0.273

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
results <- pairwise.adonis2(varespec.uni, meta_bac$Environment)
print(results)
#          Comparison  F.Model        R2 p.value p.adjusted
#Desiccated vs Pond  2.76777 0.1215653   0.030      0.090
#Desiccated vs Lake 12.05764 0.3761238   0.001      0.003
#      Pond vs Lake 11.59660 0.3918221   0.001      0.003

# PERMANOVA with Year as covariate (sequential, Type I)
meta_bac$Year <- as.factor(ifelse(meta_bac$Environment == "Lake", 2011, 2019))
set.seed(123)
m_year_env <- adonis2(varespec.uni ~ Year + Environment,
                      data = meta_bac, permutations = 999, by = "terms")
print(m_year_env)
#Df SumOfSqs      R2       F Pr(>F)    
#Year         1  0.53053 0.30012 13.7419  0.001 ***
#Environment  1  0.11761 0.06653  3.0464  0.010 ** 
set.seed(123)
m_marginal <- adonis2(varespec.uni ~ Year + Environment,
                      data = meta_bac, permutations = 999, by = "margin")
print(m_marginal) #0.01 **

# PERMANOVA at each taxonomic level, using relative-abundance Bray-Curtis
# (run once per level, using the taxonomy-collapsed abundance tables)
library(vegan)

run_level_permanova <- function(tax_level) {
  # Aggregate ASV relative abundances by named taxon at this level
  taxa_names <- tax[[tax_level]]
  agg <- aggregate(rel_abund, by = list(Taxon = taxa_names), FUN = sum)
  rownames(agg) <- agg$Taxon
  agg$Taxon <- NULL
  agg_t <- t(agg)  # samples as rows
  
  bray <- vegdist(agg_t, method = "bray")
  ad <- adonis2(bray ~ meta$Environment, permutations = 999)
  cat("\n", tax_level, "— PERMANOVA p =", ad$`Pr(>F)`[1], "\n")
  return(ad)
}

phylum_permanova <- run_level_permanova("taxonomy3")
order_permanova  <- run_level_permanova("taxonomy4")
family_permanova <- run_level_permanova("taxonomy5")
genus_permanova  <- run_level_permanova("taxonomy6")

# Print full results, not just p — R2 and F should differ across levels
# even when p is bottomed out at the permutation floor
print(phylum_permanova)
print(order_permanova)
print(family_permanova)
print(genus_permanova)

# Quick side-by-side comparison of R2 and F across levels
compare_permanova <- data.frame(
  Level = c("Phylum", "Order", "Family", "Genus"),
  R2    = c(phylum_permanova$R2[1], order_permanova$R2[1],
            family_permanova$R2[1], genus_permanova$R2[1]),
  F     = c(phylum_permanova$F[1], order_permanova$F[1],
            family_permanova$F[1], genus_permanova$F[1]),
  p     = c(phylum_permanova$`Pr(>F)`[1], order_permanova$`Pr(>F)`[1],
            family_permanova$`Pr(>F)`[1], genus_permanova$`Pr(>F)`[1])
)
print(compare_permanova)

# PCoA
ordu <- ordinate(tot, "PCoA", "unifrac", weighted = TRUE)
df   <- as.data.frame(as.matrix(ordu$vectors))
df$sample <- row.names(df)
meta_bac$Axis.1 <- df$Axis.1
meta_bac$Axis.2 <- df$Axis.2
eig2 <- ordu$values$Eigenvalues
eig2 / sum(eig2)

find_hull  <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
find_hulls <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
micro.hulls.bac <- meta_bac %>% group_by(Environment) %>% do(find_hulls(.))

pdf("Bac_Beta_Diversity_notufts.pdf", width = 7.78, height = 5)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_point(size = 4, aes(colour = Environment), alpha = 0.5) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%", title = "16S rRNA gene Beta Diversity (tufts excluded)") +
  geom_polygon(data = micro.hulls.bac, aes(colour = Environment, fill = Environment),
               alpha = 0.1, size = 0.25, show.legend = F) +
  theme_bw() +
  theme(legend.position = "right",
        plot.title    = element_text(size = 16),
        axis.title.x  = element_text(face = "bold", size = 16, vjust = 2),
        axis.text.x   = element_text(size = 14),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 16),
        plot.margin   = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# SIMPER — uses relative abundance (decostand "total"), NOT rarefied counts
# Relative abundance removes library size bias from pairwise dissimilarity
bac_asv  <- read.csv("ASV_table_rare_notufts.csv", check.names = F)
meta_bac_sim <- read.delim("Drymats_metadata_final.txt") %>%
  filter(Sample %in% colnames(bac_asv)[-1]) %>%
  column_to_rownames(var = "Sample")
# Use phyloseq otu_table which has correct ESV IDs as row names
rel_bac <- decostand(
  t(as.data.frame(as.matrix(tot@otu_table))),
  "total"
)

# Make sure meta_bac_sim row order matches
meta_bac_sim <- meta_bac_sim[match(rownames(rel_bac), rownames(meta_bac_sim)), ]
stopifnot(all(rownames(rel_bac) == rownames(meta_bac_sim)))

sim   <- with(meta_bac_sim, simper(rel_bac, Environment))
s_bac <- summary(sim)

# Check ESV IDs now
rownames(head(s_bac$Desiccated_Lake, 5))
# Should now return ESV_1, ESV_2 etc.
#"ESV_1"  "ESV_2"  "ESV_3"  "ESV_31" "ESV_5" 
# Top contributors
top_esvs <- unique(c(
  rownames(head(s_bac$Desiccated_Pond, 5)),
  rownames(head(s_bac$Desiccated_Lake, 5)),
  rownames(head(s_bac$Pond_Lake,       5))
))
top_esvs

# Pull taxonomy
tax_top <- input_filt_rare$taxonomy_loaded[top_esvs, ]
print(tax_top)

# SIMPER vectors on PCoA — update ESV numbers from new SIMPER output
tot_transposed <- t(tot@otu_table)
w    <- wascores(x = ordu$vectors, w = tot_transposed)
wdf  <- as.data.frame(w)
wdf$Species <- rownames(w)
# Update ESV IDs below based on your new SIMPER output
sub_bac <- subset(wdf, Species %in% c("ESV_1","ESV_2","ESV_10","ESV_3","ESV_31","ESV_14"))
sub_bac$shortnames <- c("Tychonema", "Leptolyngbya", "Chloronema",
                        "Rhodanobacter", "Microbacteriaceae", "RD011")
micro.hulls.bac <- ddply(meta_bac, "Environment", find_hull)
meta_bac$Axis.1 <- df$Axis.1
meta_bac$Axis.2 <- df$Axis.2
micro.hulls.bac$Environment <- factor(micro.hulls.bac$Environment)
meta_bac$Environment        <- factor(meta_bac$Environment)
max_sample_distance  <- max(abs(meta_bac$Axis.1), abs(meta_bac$Axis.2))
max_species_distance <- max(abs(sub_bac$Axis.1), abs(sub_bac$Axis.2))
scaling_factor       <- 0.5 * (max_sample_distance / max_species_distance)
sub_bac$Axis.1       <- sub_bac$Axis.1 * scaling_factor
sub_bac$Axis.2       <- sub_bac$Axis.2 * scaling_factor
label_offset         <- 1.1; perp_offset <- 0.05
sub_bac$LabelX       <- sub_bac$Axis.1 * label_offset - perp_offset * sub_bac$Axis.2
sub_bac$LabelY       <- sub_bac$Axis.2 * label_offset + perp_offset * sub_bac$Axis.1

pdf("BetaDiversity_simper_scaled_notufts.pdf", width = 7.78, height = 5)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_point(size = 4, aes(colour = Environment), alpha = 0.5) +
  geom_polygon(data = micro.hulls.bac, aes(colour = Environment, fill = Environment),
               alpha = 0.1, size = 0.25, show.legend = FALSE) +
  geom_segment(data = sub_bac, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "darkgrey", size = 0.6) +
  geom_text_repel(data = sub_bac, aes(x = LabelX, y = LabelY, label = shortnames),
                  size = 4, colour = "black", nudge_x = 0.02, nudge_y = 0.02,
                  segment.size = 0, segment.color = NA) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%", title = "16S rRNA gene Beta Diversity (tufts excluded)") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title    = element_text(size = 16),
        axis.title.x  = element_text(face = "bold", size = 16, vjust = 2),
        axis.text.x   = element_text(size = 14),
        axis.text.y   = element_text(size = 14),
        axis.title.y  = element_text(face = "bold", size = 16),
        plot.margin   = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# ── Type-level hulls ──────────────────────────────────────────────────────────
micro.hulls.type <- meta_bac %>% group_by(Type) %>% do(find_hulls(.))

# ── Plot 1: Type-level PCoA without SIMPER vectors ───────────────────────────
pdf("BetaDiversity_Type_notufts.pdf", width = 9, height = 6)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.type,
               aes(colour = Type, fill = Type),
               alpha = 0.1, linewidth = 0.25, show.legend = FALSE) +
  geom_point(size = 4, aes(colour = Type, shape = Environment), alpha = 0.8) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%",
       title = "16S rRNA gene Beta Diversity by sample type (tufts excluded)",
       colour = "Sample Type",
       shape  = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title   = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        plot.margin  = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# ── Plot 2: Type-level PCoA with SIMPER vectors ───────────────────────────────
pdf("BetaDiversity_Type_simper_notufts.pdf", width = 9, height = 6)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.type,
               aes(colour = Type, fill = Type),
               alpha = 0.1, linewidth = 0.25, show.legend = FALSE) +
  geom_point(size = 4, aes(colour = Type, shape = Environment), alpha = 0.8) +
  geom_segment(data = sub_bac,
               aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "darkgrey", linewidth = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = sub_bac,
                  aes(x = LabelX, y = LabelY, label = shortnames),
                  size = 4, colour = "black",
                  nudge_x = 0.02, nudge_y = 0.02,
                  segment.size = 0, segment.color = NA,
                  inherit.aes = FALSE) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%",
       title = "16S rRNA gene Beta Diversity by sample type (tufts excluded)",
       colour = "Sample Type",
       shape  = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title   = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title.y = element_text(face = "bold", size = 14),
        plot.margin  = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# Ensure Type and Environment are factors with meaningful order
meta_bac$Type <- factor(meta_bac$Type,
                        levels = c("flat mat", "pinnacle", "cone",
                                   "Avalanche pond", "Eastern lateral moraine pond",
                                   "Southern pond", "Western lateral moraine pond",
                                   "Untersee Desiccated Mat", "Snow Petrel Desiccated Mat"))

meta_bac$Environment <- factor(meta_bac$Environment,
                               levels = c("Lake", "Pond", "Desiccated"))

# Convex hulls by Type and Environment
micro.hulls.type <- meta_bac %>% group_by(Type) %>% do(find_hulls(.))
micro.hulls.env  <- meta_bac %>% group_by(Environment) %>% do(find_hulls(.))

# Colour palette for Type — defined at top of script, available here

# Figure 1: Type coloured, no environment hulls
pdf("BetaDiversity_Type_notufts_final.pdf", width = 9, height = 6)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.type,
               aes(colour = Type, fill = Type),
               alpha = 0.08, linewidth = 0.25, show.legend = FALSE) +
  geom_point(aes(colour = Type, shape = Environment), size = 4, alpha = 0.9) +
  scale_colour_manual(values = type_colours) +
  scale_fill_manual(values   = type_colours) +
  scale_shape_manual(values  = c("Lake" = 16, "Pond" = 17, "Desiccated" = 15)) +
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%",
       title = "16S rRNA gene - Weighted UniFrac by sample type (tufts excluded)",
       colour = "Sample Type", shape = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 13),
        axis.title.x    = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        axis.title.y    = element_text(face = "bold", size = 14),
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"),
        plot.margin     = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

# Figure 2: Type colours with dashed environment hulls overlaid
pdf("BetaDiversity_Type_EnvHulls_notufts.pdf", width = 9, height = 6)
ggplot(meta_bac, aes(Axis.1, Axis.2)) +
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
  labs(x = "PC1: 35.15%", y = "PC2: 26.58%",
       title = "16S rRNA gene - Sample types within habitat boundaries (tufts excluded)",
       colour = "Sample Type", shape = "Habitat") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 12),
        axis.title.x    = element_text(face = "bold", size = 14, vjust = 2),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        axis.title.y    = element_text(face = "bold", size = 14),
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"),
        plot.margin     = unit(c(0, 0.1, 0, 0.1), "cm"))
dev.off()

### Network analysis 16S
library(ggplot2); library(igraph); library(scales)
library(mctoolsr); library(Hmisc); library(dplyr)
library(RColorBrewer)
set.seed(500)

input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
top16s    <- as.data.frame(head(sort(rowSums(input_filt_rare$data_loaded), decreasing = TRUE), n = 200))
write.csv(top16s, "top16s_notufts.csv")
top16s    <- read.csv("top16s_notufts.csv", row.names = 1)
top200bac <- filter_taxa_from_input(input_filt_rare, taxa_IDs_to_keep = rownames(top16s))
cor.cutoff <- 0.8

# Set up taxonomy colours — shared across all three environment networks
top200bac$taxonomy_loaded$taxonomy2 <- as.character(top200bac$taxonomy_loaded$taxonomy2)
top200bac$taxonomy_loaded$taxonomy2[top200bac$taxonomy_loaded$taxonomy2 == "NA"] <- "Unassigned"
top200bac$taxonomy_loaded$taxonomy2 <- as.factor(top200bac$taxonomy_loaded$taxonomy2)
l <- levels(top200bac$taxonomy_loaded$taxonomy2)
top200bac$taxonomy_loaded$taxonomy2 <- factor(top200bac$taxonomy_loaded$taxonomy2,
                                               levels = c(l[1:11], l[13], l[12]))
colrs <- c(brewer.pal(12, "Paired"), "grey60")

build_network <- function(data, cor.cutoff) {
  cor.mat <- rcorr(t(data), type = "spearman")
  diag(cor.mat$r) <- 0
  cor.mat$r[is.na(cor.mat$r)] <- 0
  net <- graph_from_adjacency_matrix(cor.mat$r, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[abs(weight) < cor.cutoff])
  list(net = net, cor = cor.mat)
}

plot_network <- function(net, cor.mat, deg, taxonomy, colrs, filename, title) {
  V(net)$phylum <- factor(taxonomy$taxonomy2[match(V(net)$name, rownames(taxonomy))],
                           levels = levels(taxonomy$taxonomy2))
  V(net)$color  <- colrs[as.numeric(V(net)$phylum)]
  abs_w  <- abs(E(net)$weight)
  net_abs <- net; E(net_abs)$weight <- abs_w
  comm    <- cluster_fast_greedy(net_abs)
  layout_cl <- layout_with_fr(net_abs)
  pdf(filename, width = 9, height = 4.8)
  par(mar = c(0, 0, 0, 0)); par(xpd = TRUE)
  plot(net, vertex.color = V(net)$color, vertex.size = deg * 0.3,
       vertex.shape = "circle", vertex.frame.color = "black", vertex.label = NA,
       edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved = 0.2, edge.width = abs_w * 0.5, layout = layout_cl * 10)
  legend(x = -1.7, y = 0.8, levels(taxonomy$taxonomy2),
         pch = 21, col = NA, pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)
  dev.off()
}

#### 1. Lake (no tufts — filament samples already removed) ####
top200bac_lake <- filter_data(top200bac, 'Environment', keep_vals = 'Lake')
lake_data      <- top200bac$data_loaded[, rownames(top200bac_lake$map_loaded)]
write.csv(lake_data, "lake_data_notufts.csv")
lake_res  <- build_network(lake_data, cor.cutoff)
net_lake  <- lake_res$net
deg_lake  <- degree(net_lake, mode = "all")
cat("Lake network (no tufts): Edges =", length(E(net_lake)),
    "| Mean degree =", round(mean(deg_lake), 2),
    "| Transitivity =", round(transitivity(net_lake), 3),
    "| Positive =", sum(E(net_lake)$weight > 0),
    "| Negative =", sum(E(net_lake)$weight < 0), "\n")
plot_network(net_lake, lake_res$cor, deg_lake,
             top200bac$taxonomy_loaded, colrs,
             "Lake_circle16S_notufts.pdf", "Lake")

#### 2. Pond ####
top200bac_pond <- filter_data(top200bac, 'Environment', keep_vals = 'Pond')
pond_data      <- top200bac$data_loaded[, rownames(top200bac_pond$map_loaded)]
write.csv(pond_data, "pond_data_notufts.csv")
pond_res  <- build_network(pond_data, cor.cutoff)
net_pond  <- pond_res$net
deg_pond  <- degree(net_pond, mode = "all")
cat("Pond network: Edges =", length(E(net_pond)),
    "| Mean degree =", round(mean(deg_pond), 2),
    "| Transitivity =", round(transitivity(net_pond), 3),
    "| Positive =", sum(E(net_pond)$weight > 0),
    "| Negative =", sum(E(net_pond)$weight < 0), "\n")
plot_network(net_pond, pond_res$cor, deg_pond,
             top200bac$taxonomy_loaded, colrs,
             "Pond_circle16S_notufts.pdf", "Pond")

#### 3. Desiccated ####
top200bac_des <- filter_data(top200bac, 'Environment', keep_vals = 'Desiccated')
des_data      <- top200bac$data_loaded[, rownames(top200bac_des$map_loaded)]
write.csv(des_data, "des_data_notufts.csv")
des_res   <- build_network(des_data, cor.cutoff)
net_des   <- des_res$net
deg_des   <- degree(net_des, mode = "all")
cat("Desiccated network: Edges =", length(E(net_des)),
    "| Mean degree =", round(mean(deg_des), 2),
    "| Transitivity =", round(transitivity(net_des), 3),
    "| Positive =", sum(E(net_des)$weight > 0),
    "| Negative =", sum(E(net_des)$weight < 0), "\n")
plot_network(net_des, des_res$cor, deg_des,
             top200bac$taxonomy_loaded, colrs,
             "Desiccated_circle16S_notufts.pdf", "Desiccated")

### Indicator Species — All Bacteria
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
library(indicspecies)
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
  function(row) if (max(row) >= 0.5) colnames(significant_results[, c("Desiccated","Pond","Lake")])[which.max(row)] else NA)
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
ordered_data  <- ordered_data[ordered_data$Final >= 0.875, ]

pdf("Indicator_species_Bacteria_notufts.pdf", width = 6, height = 12)
ggplot(ordered_data, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.85, xend = Final - 0.015, y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("Lake" = "darkblue","Pond" = "darkgreen","Desiccated" = "brown")) +
  scale_x_continuous(limits = c(0.85, 1)) +
  labs(title = "Indicator Species by Environment (tufts excluded)", x = "Indicator Value", y = "ESV") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "right")
dev.off()

### Indicator Species — Cyanobacteria only
cyano_only    <- filter_taxa_from_input(input_filt_rare, taxa_to_keep = "Cyanobacteria")
otu_table_cy  <- cyano_only$data_loaded
metadata_cy   <- cyano_only$map_loaded
grouping_cy   <- metadata_cy[colnames(otu_table_cy), "Environment"]
indicator_cy  <- multipatt(t(otu_table_cy), grouping_cy, control = how(nperm = 999))
summary(indicator_cy)
sign_cy       <- rownames_to_column(indicator_cy$sign, var = "ESV")
sign_cy       <- sign_cy[!is.na(sign_cy$p.value), ]
indval_cy     <- rownames_to_column(as.data.frame(indicator_cy$str), var = "ESV")
indval_cy     <- indval_cy[, 1:(ncol(indval_cy) - 4)]
merged_cy     <- merge(sign_cy, indval_cy, by = "ESV")
sig_cy        <- merged_cy[merged_cy$p.value < 0.05, ]
sig_cy$Environment <- apply(sig_cy[, c("Desiccated","Pond","Lake")], 1,
  function(row) if (max(row) >= 0.5) colnames(sig_cy[, c("Desiccated","Pond","Lake")])[which.max(row)] else NA)
sig_cy$Final  <- apply(sig_cy[, c("Pond","Desiccated","Lake")], 1, max)
sig_cy_final  <- sig_cy[, c(1, which(colnames(sig_cy) %in% c("Environment","Final")))]
tax_cy        <- get_best_taxonomy(rownames_to_column(cyano_only$taxonomy_loaded, var = "ESV"))
merged_cy2    <- merge(sig_cy_final, tax_cy[, c("ESV","best_taxonomy")], by = "ESV", all.x = TRUE)
merged_cy2$Environment <- factor(merged_cy2$Environment, levels = c("Lake","Pond","Desiccated"))
ord_cy        <- merged_cy2[order(merged_cy2$Environment, merged_cy2$Final), ]
ord_cy$ESV_label <- factor(paste0(ord_cy$ESV," (",ord_cy$best_taxonomy,")"),
                            levels = rev(unique(paste0(ord_cy$ESV," (",ord_cy$best_taxonomy,")"))))
ord_cy        <- ord_cy[ord_cy$Final >= 0.7, ]

pdf("Indicator_species_cyano_notufts.pdf", width = 6, height = 6)
ggplot(ord_cy, aes(x = Final, y = ESV_label, color = Environment)) +
  geom_point(size = 6) +
  geom_segment(aes(x = 0.65, xend = Final - 0.015, y = ESV_label, yend = ESV_label),
               color = "grey", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("Lake" = "darkblue","Pond" = "darkgreen","Desiccated" = "brown")) +
  scale_x_continuous(limits = c(0.65, 0.95)) +
  labs(title = "Cyanobacteria Indicator Species (tufts excluded)", x = "Indicator Value", y = "ESV") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor   = element_blank(),
        legend.position    = "right")
dev.off()

### SourceTracker visualisation
# SourceTracker2 was run on NON-rarefied data (ASV_table_NONrare.biom) — correct
# for this algorithm which handles library size variation internally.
# Three separate analyses: each habitat designated as sink in turn.
# Results interpreted as COMPOSITIONAL OVERLAP not directional dispersal.
# Lake-involving comparisons carry additional cross-study caveat (Greco et al. 2020).

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_16S/")
meta_source <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_16S/Drymats_metadata_final.txt")

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/Sourcetracker_final/")
library(dplyr)

# ── Read mixing proportions from each sink subfolder ──────────────────────────
desic_sink <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_16S/sourcetracker_results_desiccated/mixing_proportions.txt", check.names = FALSE)
lake_sink  <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_16S/sourcetracker_results_lake/mixing_proportions.txt",check.names = FALSE)
pond_sink  <- read.delim("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/lake_pond_desiccated/Sourcetracker_final/Sourcetracker_16S/sourcetracker_results_pond/mixing_proportions.txt",check.names = FALSE)

# Add Environment from metadata to each sink
id_col <- colnames(desic_sink)[1]  # typically "SampleID" or "#SampleID"
desic_sink <- merge(desic_sink, meta_source, by.x = id_col, by.y = "Sample", all.x = TRUE)
lake_sink  <- merge(lake_sink,  meta_source, by.x = id_col, by.y = "Sample", all.x = TRUE)
pond_sink  <- merge(pond_sink,  meta_source, by.x = id_col, by.y = "Sample", all.x = TRUE)

# Verify column names — source columns should be Lake, Pond, Desiccated, Unknown
print(head(desic_sink))
print(head(lake_sink))
print(head(pond_sink))

# ── Compute average source contributions per sink ─────────────────────────────
# Desiccated as sink — primary ecologically justified analysis
# Lake and Pond are both putative sources; pond contribution is within-study (2019)
# Lake contribution carries cross-study caveat
desic_avg <- colMeans(desic_sink[, c("Lake", "Pond", "Unknown")], na.rm = TRUE)

# Pond as sink
pond_avg  <- colMeans(pond_sink[, c("Lake", "Desiccated", "Unknown")], na.rm = TRUE)

# Lake as sink — note: circular for Desiccated source (desiccated originated from lake)
# and cross-study for all comparisons; present as exploratory only
lake_avg  <- colMeans(lake_sink[, c("Pond", "Desiccated", "Unknown")], na.rm = TRUE)

# ── Visualisation: three-panel pie chart ──────────────────────────────────────
# Order: Desiccated as sink (primary), Pond as sink, Lake as sink (exploratory)

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/no_tufts_16S/")

pdf("SourceTracker_all_sinks_notufts.pdf", width = 14, height = 5)
par(mfrow = c(1, 3))

# Panel 1: Desiccated as sink (primary — within-study pond contribution)
pie(c(desic_avg["Lake"], desic_avg["Pond"], desic_avg["Unknown"]),
    labels = c(paste0("Lake\n", round(desic_avg["Lake"] * 100, 1), "%"),
               paste0("Pond\n",  round(desic_avg["Pond"]  * 100, 1), "%"),
               paste0("Unknown\n",round(desic_avg["Unknown"] * 100, 1), "%")),
    col    = c("#33A02C", "#1F78B4", "grey80"),
    main   = "Desiccated as sink\n(primary analysis)",
    cex.main = 0.95)

# Panel 2: Pond as sink
pie(c(pond_avg["Lake"], pond_avg["Desiccated"], pond_avg["Unknown"]),
    labels = c(paste0("Lake\n",       round(pond_avg["Lake"]       * 100, 1), "%"),
               paste0("Desiccated\n", round(pond_avg["Desiccated"] * 100, 1), "%"),
               paste0("Unknown\n",    round(pond_avg["Unknown"]    * 100, 1), "%")),
    col    = c("#33A02C", "#E31A1C", "grey80"),
    main   = "Pond as sink",
    cex.main = 0.95)

# Panel 3: Lake as sink (exploratory — cross-study caveat + circular Desiccated source)
pie(c(lake_avg["Pond"], lake_avg["Desiccated"], lake_avg["Unknown"]),
    labels = c(paste0("Pond\n",       round(lake_avg["Pond"]       * 100, 1), "%"),
               paste0("Desiccated*\n",round(lake_avg["Desiccated"] * 100, 1), "%"),
               paste0("Unknown\n",    round(lake_avg["Unknown"]    * 100, 1), "%")),
    col    = c("#1F78B4", "#E31A1C", "grey80"),
    main   = "Lake as sink\n(exploratory — see caveats)",
    cex.main = 0.95)

mtext("* Desiccated source reflects shared ancestry rather than dispersal (circular)",
      side = 1, line = -1, outer = TRUE, cex = 0.75, adj = 0.98)

par(mfrow = c(1, 1))
dev.off()

# ── Print summary table for manuscript ────────────────────────────────────────
cat("\n=== SourceTracker: Average compositional overlap (16S, tufts excluded) ===\n")
cat("\nDesiccated as sink (primary analysis):\n")
print(round(desic_avg * 100, 1))
cat("\nPond as sink:\n")
print(round(pond_avg * 100, 1))
cat("\nLake as sink (exploratory):\n")
print(round(lake_avg * 100, 1))

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

# =============================================================================
# CORRECTED PHYLUM-LEVEL STATISTICS — true phylum (taxonomy2), with
# Proteobacteria split into its component classes, to match the corrected
# Figure 2. Replaces all earlier taxonomy3-based "phylum" statistics.
# =============================================================================

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")

# --- Build the custom phylum column (same logic as the Figure 2 fix) ---
input_filt_rare$taxonomy_loaded$taxonomy2_custom <- as.character(input_filt_rare$taxonomy_loaded$taxonomy2)
proteo_rows <- input_filt_rare$taxonomy_loaded$taxonomy2 == "Proteobacteria"
input_filt_rare$taxonomy_loaded$taxonomy2_custom[proteo_rows] <-
  as.character(input_filt_rare$taxonomy_loaded$taxonomy3[proteo_rows])

cat("Custom phylum-level categories:\n")
print(sort(table(input_filt_rare$taxonomy_loaded$taxonomy2_custom), decreasing = TRUE))

# Overwrite taxonomy2 so existing functions (which reference "taxonomy2" by
# name) automatically pick up the custom, Proteobacteria-split categories
input_filt_rare$taxonomy_loaded$taxonomy2 <- input_filt_rare$taxonomy_loaded$taxonomy2_custom

# --- Rebuild rel_abund from this object, since it must reflect the same
#     taxonomy table get_env_stats()/test_taxon() will reference ---
rel_abund <- sweep(input_filt_rare$data_loaded, 2, colSums(input_filt_rare$data_loaded), "/")
tax <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded  # confirm this matches whatever object name get_env_stats()/test_taxon() actually use internally — check against their body if this errors

# =============================================================================
# Re-run summary statistics (mean ± SD) at TRUE phylum level (taxonomy2)
# =============================================================================
phylum_results_corrected <- get_all_level_stats("taxonomy2")
write.csv(phylum_results_corrected, "manuscript_taxa_summary_notufts_PHYLUM_CORRECTED.csv", row.names = FALSE)

# =============================================================================
# Re-run Kruskal-Wallis + Dunn tests at TRUE phylum level (taxonomy2)
# =============================================================================
phylum_tests_corrected <- get_all_level_tests("taxonomy2")
write.csv(phylum_tests_corrected, "manuscript_stats_ALL_taxa_notufts_PHYLUM_CORRECTED.csv", row.names = FALSE)

cat("\nDone. New files:\n")
cat(" - manuscript_taxa_summary_notufts_PHYLUM_CORRECTED.csv\n")
cat(" - manuscript_stats_ALL_taxa_notufts_PHYLUM_CORRECTED.csv\n")
cat("\nThese now report TRUE phylum-level statistics (Cyanobacteria, Bacteroidota,\n")
cat("Actinobacteriota, Planctomycetota, Verrucomicrobiota, etc.), with\n")
cat("Proteobacteria split into Alphaproteobacteria / Betaproteobacteria /\n")
cat("Gammaproteobacteria as separate categories — consistent with Figure 2.\n")