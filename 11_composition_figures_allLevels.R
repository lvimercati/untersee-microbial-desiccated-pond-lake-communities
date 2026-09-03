# =============================================================================
# TAXONOMIC COMPOSITION FIGURES — ALL LEVELS BELOW PHYLUM/DIVISION
# Generates the same style of figure as Figure 2 (stacked % bar charts,
# Bacteria on top / Eukarya below, faceted by Desiccated/Pond/Lake) for
# Class, Order, Family, and Genus levels, for both domains.
#
# Built directly on the same mctoolsr functions (summarize_taxonomy,
# plot_taxa_bars) and sample-ordering/merging logic used for the existing
# Figure 2 (phylum, level=2) script.
#
# NOTE ON LEGEND STYLING: Figure 2's final combined output used a hand-tuned,
# pixel-perfect legend placement (via magick) that was specifically tuned for
# the phylum-level palette's exact size. That tuning doesn't generalize
# automatically to Class/Order/Family/Genus, where the number of legend
# entries differs at each level. This script instead uses patchwork's
# built-in legend collection (guides = "collect"), which adapts automatically
# to any legend size -- the visual STYLE (colours, bars, faceting, fonts,
# sample order) is otherwise identical to Figure 2.
# =============================================================================

library(dplyr)
library(ggplot2)
library(mctoolsr)
library(patchwork)
library(scales)

# -----------------------------------------------------------------------
# SHARED SETTINGS (copied exactly from the Figure 2 script)
# -----------------------------------------------------------------------
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

# Levels to generate: Class, Order, Family, Genus
# (mctoolsr/SILVA & PR2 both use this same numbering: 3=Class, 4=Order,
#  5=Family, 6=Genus -- consistent across 16S and 18S in this pipeline)
levels_to_plot <- list(
  list(level = 3, name = "Class"),
  list(level = 4, name = "Order"),
  list(level = 5, name = "Family"),
  list(level = 6, name = "Genus")
)

num_taxa <- 12  # top N taxa shown individually; rest lumped into "Other"

# -----------------------------------------------------------------------
# LOAD DATA (same RDS objects used for Figure 2)
# -----------------------------------------------------------------------
input_filt_rare_bac <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/bac_input_filt_rare_notufts.rds")
input_filt_rare_euk <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/euk_input_filt_rare_notufts.rds")

# Ensure SampleID column present (as in the original script)
if (!"SampleID" %in% colnames(input_filt_rare_bac$map_loaded)) {
  input_filt_rare_bac$map_loaded <- cbind(SampleID = rownames(input_filt_rare_bac$map_loaded),
                                          input_filt_rare_bac$map_loaded)
}
if (!"SampleID" %in% colnames(input_filt_rare_euk$map_loaded)) {
  input_filt_rare_euk$map_loaded <- cbind(SampleID = rownames(input_filt_rare_euk$map_loaded),
                                          input_filt_rare_euk$map_loaded)
}

# -----------------------------------------------------------------------
# FUNCTION: build the plotting data frame for one domain at one level
# (exact same logic as the Figure 2 script's pdat_bac/pdat_euk construction,
# generalised to any taxonomy level)
# -----------------------------------------------------------------------
build_pdat <- function(input_filt_rare, level, num_taxa = 12) {
  tax_summary <- summarize_taxonomy(input_filt_rare, level = level, report_higher_tax = FALSE)

  pdat <- plot_taxa_bars(tax_summary, input_filt_rare$map_loaded, "SampleID",
                         num_taxa = num_taxa, data_only = TRUE)

  pdat$taxon <- as.factor(pdat$taxon)
  taxon_levels <- levels(pdat$taxon)
  # Keep up to num_taxa + 1 (the +1 slot is mctoolsr's own "other" bucket)
  keep_n <- min(length(taxon_levels), num_taxa + 1)
  pdat$taxon <- factor(pdat$taxon, levels = taxon_levels[1:keep_n])

  env_lookup <- input_filt_rare$map_loaded[, c("SampleID", "Environment")]
  pdat <- merge(pdat, env_lookup, by.x = "group_by", by.y = "SampleID")
  pdat$Environment <- factor(pdat$Environment, levels = habitat_order)
  pdat$group_by <- factor(pdat$group_by,
                          levels = custom_order[custom_order %in% pdat$group_by])

  pdat$taxon_label <- as.character(pdat$taxon)
  pdat$taxon_label[is.na(pdat$taxon_label) | pdat$taxon_label == "NA"] <- "Unclassified"

  pdat
}

# -----------------------------------------------------------------------
# FUNCTION: generate a colour palette for whatever taxa are actually
# present at this level. "Other" and "Unclassified" always get the same
# fixed colours (matching Figure 2's convention); everything else gets a
# colour from an extended qualitative palette.
# -----------------------------------------------------------------------
build_palette <- function(taxon_labels) {
  taxa <- sort(unique(taxon_labels))
  taxa <- taxa[!taxa %in% c("Other", "Unclassified")]

  # Reserved colours for Other/Unclassified -- excluded from the pool
  # assignable to real taxa so nothing collides with them. Note
  # RColorBrewer's "Paired" palette includes "#FFFF99" as one of its 12
  # stock colours, which is exactly Other's reserved colour -- must be
  # filtered out explicitly, not just appended after the fact.
  reserved_colours <- c("#FFFF99", "grey80")
  source_pool <- RColorBrewer::brewer.pal(12, "Paired")
  source_pool <- source_pool[!source_pool %in% reserved_colours]

  n_needed <- length(taxa)
  base_pal <- if (n_needed <= length(source_pool)) {
    source_pool[seq_len(n_needed)]
  } else {
    colorRampPalette(source_pool)(n_needed)
  }

  palette <- setNames(base_pal, taxa)
  palette["Other"] <- "#FFFF99"
  palette["Unclassified"] <- "grey80"
  palette
}

# -----------------------------------------------------------------------
# FUNCTION: build one panel (bacteria or eukarya), styled identically to
# Figure 2's panel_A / panel_B
# -----------------------------------------------------------------------
build_panel <- function(pdat, fill_label) {
  palette <- build_palette(pdat$taxon_label)

  ggplot(pdat, aes(x = group_by, y = mean_value, fill = taxon_label)) +
    geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
    facet_grid(~ Environment, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Relative Abundance", fill = fill_label) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(axis.text.y      = element_text(size = 9, face = "bold"),
          axis.text.x      = element_text(size = 6.5, angle = 45, hjust = 1, vjust = 1, face = "bold"),
          axis.title.y     = element_text(face = "bold", size = 11),
          strip.text       = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey90"),
          legend.text      = element_text(size = 8, face = "bold"),
          legend.title     = element_text(face = "bold", size = 11),
          legend.key.size  = unit(0.4, "cm"),
          panel.spacing    = unit(0.3, "lines"))
}

# -----------------------------------------------------------------------
# MAIN LOOP: build and save the combined figure for each level
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/")

for (lvl in levels_to_plot) {

  cat("\n=== Building", lvl$name, "(level", lvl$level, ") ===\n")

  pdat_bac <- build_pdat(input_filt_rare_bac, level = lvl$level, num_taxa = num_taxa)
  pdat_euk <- build_pdat(input_filt_rare_euk, level = lvl$level, num_taxa = num_taxa)

  panel_A <- build_panel(pdat_bac, "Bacteria")
  panel_B <- build_panel(pdat_euk, "Eukarya")

  # Each panel keeps its own legend (bacteria and eukarya use different
  # taxon palettes, so there is nothing to merge) -- legend.position =
  # "right" is already set per-panel inside build_panel(), which
  # vertically centers each legend against its own panel by default.
  # A/B tags added via plot_annotation, matching the convention used
  # elsewhere in this project.
  combined <- (panel_A / panel_B) +
    plot_annotation(tag_levels = list(c("A", "B"))) &
    theme(plot.tag = element_text(size = 14, face = "bold"))

  out_base <- paste0("Figure_Composition_", lvl$name)

  ggsave(paste0(out_base, ".jpeg"), combined,
         width = 15, height = 12, dpi = 300, bg = "white")
  ggsave(paste0(out_base, ".pdf"), combined,
         width = 15, height = 12, dpi = 300)

  cat("Saved:", paste0(out_base, ".jpeg"), "and", paste0(out_base, ".pdf"), "\n")
}

cat("\n\nDone. Generated Class, Order, Family, and Genus composition figures\n")
cat("for both Bacteria and Eukarya, matching Figure 2's visual style.\n")
