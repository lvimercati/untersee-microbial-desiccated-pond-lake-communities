# =============================================================================
# MASTER SCRIPT: Sources both Mantel_16S.R and Mantel_18S.R (both corrected
# so "Mat" displays as "Desiccated" throughout), in ISOLATED environments to
# avoid variable name collisions, then builds the final combined Figure 7
# (16S left/A, 18S right/B, shared legend).
# =============================================================================

library(ggplot2)
library(cowplot)
library(vegan)
library(geosphere)
library(dplyr)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
set.seed(123)

# -----------------------------------------------------------------------
# STEP 1: Source each corrected script in its OWN environment — fully
# isolates their identically-named variables (plot_data_combined,
# combined_mantel, coord_filtered, etc.)
# -----------------------------------------------------------------------
cat("Sourcing 16S Mantel script...\n")
env_16s <- new.env()
source("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/Mantel/Mantel_16S.R",
       local = env_16s)

cat("\nSourcing 18S Mantel script...\n")
env_18s <- new.env()
source("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/Mantel/Mantel_18S.R",
       local = env_18s)

cat("\nBoth scripts sourced successfully.\n")

# -----------------------------------------------------------------------
# STEP 2: Extract needed objects, and print the FULL results tables so
# every stated value in the manuscript can be directly verified against
# these console outputs before finalizing any Results text.
# -----------------------------------------------------------------------
plot_data_16s     <- env_16s$plot_data_combined
mantel_result_16s <- env_16s$combined_mantel

plot_data_18s     <- env_18s$plot_data_combined
mantel_result_18s <- env_18s$combined_mantel

if (is.null(plot_data_16s) || is.null(mantel_result_16s)) {
  stop("Failed to extract 16S plot data — check that Mantel_16S.R defines plot_data_combined and combined_mantel.")
}
if (is.null(plot_data_18s) || is.null(mantel_result_18s)) {
  stop("Failed to extract 18S plot data — check that Mantel_18S.R defines plot_data_combined and combined_mantel.")
}

cat("\n========== FULL 16S MANTEL RESULTS TABLE ==========\n")
print(env_16s$mantel_results)

cat("\n========== FULL 18S MANTEL RESULTS TABLE ==========\n")
print(env_18s$mantel_results)

cat("\n16S: n =", nrow(plot_data_16s), "pairwise comparisons\n")
cat("18S: n =", nrow(plot_data_18s), "pairwise comparisons\n")

# -----------------------------------------------------------------------
# STEP 3: Build both panels with consistent, updated styling
# (labels now read "Desiccated-Desiccated", "Desiccated-Pond", "Pond-Pond")
# -----------------------------------------------------------------------
pair_colors <- c("Desiccated-Desiccated" = "#A23B72", "Desiccated-Pond" = "#F18F01", "Pond-Pond" = "#2E86AB")

p2_16s <- ggplot(plot_data_16s, aes(x = Geographic_Distance_km, y = Community_Dissimilarity,
                                     color = Habitat_Pair)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = pair_colors, na.value = "grey50") +
  theme_minimal() +
  labs(title = "A) Bacterial (16S rRNA gene) communities",
       subtitle = paste0("Mantel r = ", round(mantel_result_16s$result$Mantel_r, 3),
                         ", p = ", mantel_result_16s$result$P_value),
       x = "Geographic Distance (km)",
       y = "Community Dissimilarity (Bray-Curtis)",
       color = "Habitat Pair") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, face = "bold"),
    axis.title    = element_text(face = "bold"),
    legend.title  = element_text(face = "bold"),
    legend.text   = element_text(face = "bold"),
    legend.position = "bottom"
  )

p2_18s <- ggplot(plot_data_18s, aes(x = Geographic_Distance_km, y = Community_Dissimilarity,
                                     color = Habitat_Pair)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = pair_colors, na.value = "grey50") +
  theme_minimal() +
  labs(title = "B) Eukaryotic (18S rRNA gene) communities",
       subtitle = paste0("Mantel r = ", round(mantel_result_18s$result$Mantel_r, 3),
                         ", p = ", mantel_result_18s$result$P_value),
       x = "Geographic Distance (km)",
       y = "Community Dissimilarity (Bray-Curtis)",
       color = "Habitat Pair") +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, face = "bold"),
    axis.title    = element_text(face = "bold"),
    legend.title  = element_text(face = "bold"),
    legend.text   = element_text(face = "bold"),
    legend.position = "bottom"
  )

# -----------------------------------------------------------------------
# STEP 4: Combine side by side with ONE shared legend
# -----------------------------------------------------------------------
combined_figure7 <- plot_grid(
  p2_16s + theme(legend.position = "none"),
  p2_18s + theme(legend.position = "none"),
  ncol = 2
)

shared_legend <- cowplot::get_legend(p2_16s + theme(legend.position = "bottom"))

figure7_final <- plot_grid(combined_figure7, shared_legend,
                           ncol = 1, rel_heights = c(1, 0.08))

# -----------------------------------------------------------------------
# STEP 5: Save — PDF and high-res JPEG
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Bacteria/Mantel/")

ggsave("Figure7_Combined_DistanceDecay_Desiccated.pdf", plot = figure7_final,
       width = 16, height = 7, dpi = 300)

ggsave("Figure7_Combined_DistanceDecay_Desiccated.jpeg", plot = figure7_final,
       width = 16, height = 7, dpi = 300, quality = 100)

cat("\nSaved: Figure7_Combined_DistanceDecay_Desiccated.pdf and .jpeg\n")

cat("\n=== IMPORTANT: verify these printed values against the manuscript text ===\n")
cat("16S combined: r =", round(mantel_result_16s$result$Mantel_r, 3),
    ", p =", mantel_result_16s$result$P_value, "\n")
cat("18S combined: r =", round(mantel_result_18s$result$Mantel_r, 3),
    ", p =", mantel_result_18s$result$P_value, "\n")
