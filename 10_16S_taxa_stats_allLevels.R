# =============================================================================
# BACTERIAL (16S) TAXA SUMMARY AND STATISTICS — ALL TAXONOMIC LEVELS
# Standalone, self-contained. Generates BOTH abundance summaries (mean ± SD
# per habitat) AND Kruskal-Wallis + Dunn statistical tests at every level:
# Phylum (taxonomy2, with Proteobacteria split into its component classes,
# matching the corrected Figure 2), Class (taxonomy3), Order (taxonomy4),
# Family (taxonomy5), and Genus (taxonomy6).
# =============================================================================

library(mctoolsr)
library(dplyr)
library(FSA)

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

# -----------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')

cat("Samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

# -----------------------------------------------------------------------
# 2. BUILD CORRECTED PHYLUM COLUMN — true phylum (taxonomy2), with
#    Proteobacteria split into its component classes (Alpha-, Beta-,
#    Gammaproteobacteria), matching the corrected Figure 2 and the
#    phylum-level statistics fix applied earlier in this project.
# -----------------------------------------------------------------------
input_filt_rare$taxonomy_loaded$taxonomy2_custom <- as.character(input_filt_rare$taxonomy_loaded$taxonomy2)
proteo_rows <- input_filt_rare$taxonomy_loaded$taxonomy2 == "Proteobacteria"
input_filt_rare$taxonomy_loaded$taxonomy2_custom[proteo_rows] <-
  as.character(input_filt_rare$taxonomy_loaded$taxonomy3[proteo_rows])
input_filt_rare$taxonomy_loaded$taxonomy2 <- input_filt_rare$taxonomy_loaded$taxonomy2_custom

cat("\nCorrected phylum-level categories (Proteobacteria split):\n")
print(sort(table(input_filt_rare$taxonomy_loaded$taxonomy2), decreasing = TRUE))

# -----------------------------------------------------------------------
# 3. SHARED OBJECTS
# -----------------------------------------------------------------------
rel_abund <- sweep(input_filt_rare$data_loaded, 2,
                   colSums(input_filt_rare$data_loaded), "/")
tax  <- input_filt_rare$taxonomy_loaded
meta <- input_filt_rare$map_loaded

# -----------------------------------------------------------------------
# 4. FUNCTIONS — abundance summary (mean ± SD per habitat)
# -----------------------------------------------------------------------
get_env_stats <- function(taxon_name, tax_level) {
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

get_all_level_stats <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  do.call(rbind, lapply(taxa_at_level, get_env_stats, tax_level = tax_level))
}

# -----------------------------------------------------------------------
# 5. FUNCTIONS — Kruskal-Wallis + Dunn pairwise tests
# -----------------------------------------------------------------------
test_taxon_df <- function(taxon_name, tax_level) {
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

get_all_level_tests <- function(tax_level) {
  taxa_at_level <- unique(na.omit(tax[[tax_level]]))
  taxa_at_level <- taxa_at_level[!taxa_at_level %in% c("NA", "Unknown", "Unclassified", "", "Other")]
  do.call(rbind, lapply(taxa_at_level, test_taxon_df, tax_level = tax_level))
}

# -----------------------------------------------------------------------
# 6. RUN — SUMMARY (mean ± SD) AT ALL FIVE LEVELS
# -----------------------------------------------------------------------
cat("\n=== Generating abundance summaries ===\n")
phylum_summary <- get_all_level_stats("taxonomy2")   # corrected phylum (Proteobacteria split)
class_summary  <- get_all_level_stats("taxonomy3")   # class
order_summary  <- get_all_level_stats("taxonomy4")   # order
family_summary <- get_all_level_stats("taxonomy5")   # family
genus_summary  <- get_all_level_stats("taxonomy6")   # genus

all_summary <- rbind(phylum_summary, class_summary, order_summary, family_summary, genus_summary)
write.csv(all_summary, "manuscript_taxa_summary_notufts_16S_COMPLETE.csv", row.names = FALSE)
cat("Saved: manuscript_taxa_summary_notufts_16S_COMPLETE.csv\n")

# -----------------------------------------------------------------------
# 7. RUN — STATISTICAL TESTS AT ALL FIVE LEVELS
# -----------------------------------------------------------------------
cat("\n=== Running statistical tests ===\n")
phylum_tests <- get_all_level_tests("taxonomy2")
class_tests  <- get_all_level_tests("taxonomy3")
order_tests  <- get_all_level_tests("taxonomy4")
family_tests <- get_all_level_tests("taxonomy5")
genus_tests  <- get_all_level_tests("taxonomy6")

all_tests <- rbind(phylum_tests, class_tests, order_tests, family_tests, genus_tests)
write.csv(all_tests, "manuscript_stats_ALL_taxa_notufts_16S_COMPLETE.csv", row.names = FALSE)
cat("Saved: manuscript_stats_ALL_taxa_notufts_16S_COMPLETE.csv\n")

# -----------------------------------------------------------------------
# 8. SUMMARY OF TAXA TESTED PER LEVEL
# -----------------------------------------------------------------------
cat("\nTaxa tested per level:\n")
for (lvl in c("taxonomy2", "taxonomy3", "taxonomy4", "taxonomy5", "taxonomy6")) {
  n <- length(unique(all_tests$Taxon[all_tests$Tax_Level == lvl]))
  cat(" ", lvl, ":", n, "taxa\n")
}

# -----------------------------------------------------------------------
# 9. FLAG SIGNIFICANT TAXA AT THE PREVIOUSLY UN-REPORTED LEVELS
#    (Class = taxonomy3, Family = taxonomy5) — for direct review before
#    incorporating into the Results paragraph
# -----------------------------------------------------------------------
for (lvl_name in c("Class (taxonomy3)", "Family (taxonomy5)")) {
  lvl <- ifelse(lvl_name == "Class (taxonomy3)", "taxonomy3", "taxonomy5")
  cat("\n=== Significant taxa at", lvl_name, "(KW p < 0.05) ===\n")
  sig <- all_tests[all_tests$Tax_Level == lvl & is.na(all_tests$Comparison), ]
  sig <- sig[!is.na(sig$KW_p) & sig$KW_p < 0.05, ]
  print(sig[order(sig$KW_p), c("Taxon", "KW_p")])
}

cat("\nDone. Review the Class and Family sections above, and cross-check\n")
cat("against manuscript_taxa_summary_notufts_16S_COMPLETE.csv for abundance\n")
cat("direction (mean %) before writing any new Results sentences — do not\n")
cat("rely on significance (p-value) alone, as very low-abundance taxa can\n")
cat("appear statistically significant despite a biologically negligible or\n")
cat("even reversed effect direction.\n")
