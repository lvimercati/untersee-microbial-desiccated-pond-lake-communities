##18S
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/Eukarya/Mantel/")
# Antarctic Microbial Community Mantel Test Analysis
# Comprehensive analysis with separate and combined habitat testing

# Antarctic 18S rRNA Microbial Community Mantel Test Analysis
# Analysis for samples present in both coordinate and ASV tables

# Load required libraries
library(vegan)
library(geosphere)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# Set up results storage
mantel_results <- data.frame()

# Read 18S ASV data first to see what samples we have
asv_data <- read.csv("18SASV.csv", row.names = 1, check.names = FALSE)
print(paste("18S ASV data dimensions:", nrow(asv_data), "ASVs x", ncol(asv_data), "samples"))
print("Available samples in ASV data:")
print(colnames(asv_data))

# Create coordinate data from updated Mantel_samples.csv
coord_data <- data.frame(
  SampleID = c(
    "Untersee Desiccated Mat 1",
    "Untersee Desiccated Mat 2",
    "Untersee Desiccated Mat 3",
    "Untersee Desiccated Mat 4",
    "Untersee Desiccated Mat 5",
    "Untersee Desiccated Mat 6",
    "Untersee Desiccated Mat 7",
    "Snow Petrel Desiccated Mat 2",
    "Snow Petrel Desiccated Mat 3",
    "Snow Petrel Desiccated Mat 4",
    "Snow Petrel Desiccated Mat 5",
    "Snow Petrel Desiccated Mat 6",
    "Eastern lateral moraine pond 1",
    "Eastern lateral moraine pond 2",
    "Eastern lateral moraine pond 3",
    "Avalanche pond 1",
    "Avalanche pond 2",
    "Western lateral moraine pond 1",
    "Western lateral moraine pond 2",
    "Western lateral moraine pond 3",
    "Southern pond 1",
    "Southern pond 3"
  ),
  Latitude = c(
    -71.339250, -71.353611, -71.354750, -71.354750, -71.354722, -71.356444, -71.354861,
    -71.319639, -71.319472, -71.319389, -71.319917, -71.319833,
    -71.306778, -71.306778, -71.306778, -71.361781, -71.361781,
    -71.312300, -71.312300, -71.312300, -71.367470, -71.367470
  ),
  Longitude = c(
    13.457861, 13.404972, 13.403417, 13.403417, 13.403056, 13.445389, 13.455861,
    13.455194, 13.455194, 13.455111, 13.455167, 13.455167,
    13.572833, 13.572833, 13.572833, 13.443259, 13.443259,
    13.449400, 13.449400, 13.449400, 13.402610, 13.402610
  ),
  Habitat = c(
    "Desiccated", "Desiccated", "Desiccated", "Desiccated", "Desiccated", "Desiccated", "Desiccated",
    "Desiccated", "Desiccated", "Desiccated", "Desiccated", "Desiccated",
    "Pond", "Pond", "Pond", "Pond", "Pond",
    "Pond", "Pond", "Pond", "Pond", "Pond"
  ),
  stringsAsFactors = FALSE
)

# Get matching samples between both datasets
coord_samples <- coord_data$SampleID
asv_samples <- colnames(asv_data)
matching_samples <- intersect(coord_samples, asv_samples)

print("\nSample matching results:")
print(paste("Samples in coordinates:", length(coord_samples)))
print(paste("Samples in ASV data:", length(asv_samples)))
print(paste("Matching samples:", length(matching_samples)))

cat("\nMatching samples:\n")
print(matching_samples)

cat("\nSamples in coordinates but NOT in ASV data:\n")
missing_in_asv <- setdiff(coord_samples, asv_samples)
print(missing_in_asv)

cat("\nSamples in ASV data but NOT in coordinates (will be excluded):\n")
missing_in_coords <- setdiff(asv_samples, coord_samples)
print(missing_in_coords)

if(length(matching_samples) < 4) {
  stop("Not enough matching samples for Mantel test (need at least 4)")
}

# Filter data to matching samples
coord_filtered <- coord_data[coord_data$SampleID %in% matching_samples, ]
asv_filtered <- asv_data[, matching_samples]

# Ensure same order
coord_filtered <- coord_filtered[match(matching_samples, coord_filtered$SampleID), ]
asv_filtered <- asv_filtered[, matching_samples]

print(paste("\nFinal dataset:", nrow(coord_filtered), "samples"))
print("Final habitat distribution:")
table(coord_filtered$Habitat)

# Prepare community data
asv_t <- t(asv_filtered)  # Transpose: samples as rows, ASVs as columns
asv_t <- asv_t[, colSums(asv_t) > 0]  # Remove zero-sum ASVs
print(paste("18S ASVs retained:", ncol(asv_t)))

# Calculate distance matrices
coords_matrix <- as.matrix(coord_filtered[, c("Longitude", "Latitude")])
geo_dist <- distm(coords_matrix, fun = distHaversine) / 1000  # km
rownames(geo_dist) <- colnames(geo_dist) <- coord_filtered$SampleID
geo_dist_obj <- as.dist(geo_dist)

# Community distance matrix
comm_dist <- vegdist(asv_t, method = "bray")

# Function to perform and store Mantel test results
perform_mantel <- function(geo_dist, comm_dist, analysis_name, sample_ids, habitat_info = NULL) {
  mantel_result <- mantel(geo_dist, comm_dist, method = "pearson", permutations = 999)
  
  result_row <- data.frame(
    Analysis = analysis_name,
    N_samples = length(sample_ids),
    Mantel_r = round(mantel_result$statistic, 4),
    P_value = round(mantel_result$signif, 4),
    Significant = mantel_result$signif < 0.05,
    Habitat_info = ifelse(is.null(habitat_info), "All", habitat_info),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== ", analysis_name, " ===\n")
  cat("Samples:", length(sample_ids), "\n")
  cat("Mantel r:", round(mantel_result$statistic, 4), "\n")
  cat("P-value:", round(mantel_result$signif, 4), "\n")
  cat("Significant:", ifelse(mantel_result$signif < 0.05, "YES", "NO"), "\n")
  
  return(list(result = result_row, mantel_obj = mantel_result))
}

# 1. COMBINED ANALYSIS (All samples)
cat("\n", rep("=", 50), "\n")
cat("COMBINED ANALYSIS (ALL SAMPLES)\n")
cat(rep("=", 50), "\n")

combined_mantel <- perform_mantel(geo_dist_obj, comm_dist, "Combined (All samples)", 
                                  coord_filtered$SampleID, "All habitats")
mantel_results <- rbind(mantel_results, combined_mantel$result)

# 2. SEPARATE ANALYSES BY HABITAT
cat("\n", rep("=", 50), "\n")
cat("SEPARATE HABITAT ANALYSES\n")
cat(rep("=", 50), "\n")

# Mats only
mat_samples <- coord_filtered$SampleID[coord_filtered$Habitat == "Desiccated"]
if(length(mat_samples) >= 4) {
  mat_indices <- which(coord_filtered$SampleID %in% mat_samples)
  geo_dist_mat <- as.dist(geo_dist[mat_indices, mat_indices])
  comm_dist_mat <- vegdist(asv_t[mat_samples, ], method = "bray")
  
  mat_mantel <- perform_mantel(geo_dist_mat, comm_dist_mat, "Mat only", 
                               mat_samples, "Desiccated")
  mantel_results <- rbind(mantel_results, mat_mantel$result)
} else {
  cat("\nNot enough mat samples for analysis (", length(mat_samples), "< 4)\n")
}

# Ponds only
pond_samples <- coord_filtered$SampleID[coord_filtered$Habitat == "Pond"]
if(length(pond_samples) >= 4) {
  pond_indices <- which(coord_filtered$SampleID %in% pond_samples)
  geo_dist_pond <- as.dist(geo_dist[pond_indices, pond_indices])
  comm_dist_pond <- vegdist(asv_t[pond_samples, ], method = "bray")
  
  pond_mantel <- perform_mantel(geo_dist_pond, comm_dist_pond, "Pond only", 
                                pond_samples, "Pond")
  mantel_results <- rbind(mantel_results, pond_mantel$result)
} else {
  cat("\nNot enough pond samples for analysis (", length(pond_samples), "< 4)\n")
}

# Print results summary
cat("\n", rep("=", 50), "\n")
cat("RESULTS SUMMARY\n")
cat(rep("=", 50), "\n")
print(mantel_results)

# VISUALIZATIONS
cat("\nGenerating visualizations...\n")

# Create color palette
habitat_colors <- c("Pond" = "#2E86AB", "Desiccated" = "#A23B72")

# Initialize plot list
plot_list <- list()

# 1. Sample locations map
p1 <- ggplot(coord_filtered, aes(x = Longitude, y = Latitude, color = Habitat)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = SampleID), vjust = -0.7, size = 2.5, angle = 45) +
  scale_color_manual(values = habitat_colors) +
  theme_minimal() +
  labs(title = "Sample Locations by Habitat Type",
       subtitle = paste("Total samples:", nrow(coord_filtered)),
       x = "Longitude (°E)", y = "Latitude (°S)") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

plot_list[[1]] <- p1

# 2. Distance-decay plot - Combined analysis with separate habitat pair lines
geo_vec <- as.vector(geo_dist_obj)
comm_vec <- as.vector(comm_dist)

# Create habitat pairs
n_samples <- nrow(coord_filtered)
sample_indices <- which(lower.tri(matrix(1, n_samples, n_samples)), arr.ind = TRUE)

habitat_pairs <- character(length(geo_vec))
for(i in 1:length(geo_vec)) {
  row_idx <- sample_indices[i, 1]
  col_idx <- sample_indices[i, 2]
  hab1 <- coord_filtered$Habitat[row_idx]
  hab2 <- coord_filtered$Habitat[col_idx]
  habitat_pairs[i] <- paste(sort(c(hab1, hab2)), collapse = "-")
}

plot_data_combined <- data.frame(
  Geographic_Distance_km = geo_vec,
  Community_Dissimilarity = comm_vec,
  Habitat_Pair = factor(habitat_pairs),
  stringsAsFactors = FALSE
)

# Create color palette for habitat pairs
pair_colors <- c("Desiccated-Desiccated" = "#A23B72", "Desiccated-Pond" = "#F4A259", "Pond-Pond" = "#2E86AB")

p2 <- ggplot(plot_data_combined, aes(x = Geographic_Distance_km, y = Community_Dissimilarity, 
                                     color = Habitat_Pair)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = pair_colors, na.value = "grey50") +
  theme_minimal() +
  labs(title = "Distance-Decay Relationship: Combined Analysis (18S rRNA)",
       subtitle = paste("Mantel r =", round(combined_mantel$result$Mantel_r, 2), 
                        ", p =", combined_mantel$result$P_value),
       x = "Geographic Distance (km)",
       y = "Community Dissimilarity (Bray-Curtis)",
       color = "Habitat Pair") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

plot_list[[2]] <- p2

# 3. Ponds only distance-decay plot
if(exists("pond_mantel")) {
  geo_vec_pond <- as.vector(geo_dist_pond)
  comm_vec_pond <- as.vector(comm_dist_pond)
  
  plot_data_pond <- data.frame(
    Geographic_Distance_km = geo_vec_pond,
    Community_Dissimilarity = comm_vec_pond
  )
  
  p3 <- ggplot(plot_data_pond, aes(x = Geographic_Distance_km, y = Community_Dissimilarity)) +
    geom_point(color = habitat_colors["Pond"], alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = habitat_colors["Pond"]) +
    theme_minimal() +
    labs(title = "Distance-Decay: Pond Only (18S rRNA)",
         subtitle = paste("Mantel r =", round(pond_mantel$result$Mantel_r, 3), 
                          ", p =", pond_mantel$result$P_value,
                          ", n =", pond_mantel$result$N_samples),
         x = "Geographic Distance (km)",
         y = "Community Dissimilarity") +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  plot_list[[3]] <- p3
}

# 4. Mats only distance-decay plot
if(exists("mat_mantel")) {
  geo_vec_mat <- as.vector(geo_dist_mat)
  comm_vec_mat <- as.vector(comm_dist_mat)
  
  plot_data_mat <- data.frame(
    Geographic_Distance_km = geo_vec_mat,
    Community_Dissimilarity = comm_vec_mat
  )
  
  p4 <- ggplot(plot_data_mat, aes(x = Geographic_Distance_km, y = Community_Dissimilarity)) +
    geom_point(color = habitat_colors["Desiccated"], alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = habitat_colors["Desiccated"]) +
    theme_minimal() +
    labs(title = "Distance-Decay: Mat Only (18S rRNA)",
         subtitle = paste("Mantel r =", round(mat_mantel$result$Mantel_r, 2), 
                          ", p =", mat_mantel$result$P_value,
                          ", n =", mat_mantel$result$N_samples),
         x = "Geographic Distance (km)",
         y = "Community Dissimilarity") +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  plot_list[[4]] <- p4
}

# 5. NMDS Ordination
cat("Running NMDS ordination...\n")
nmds <- metaMDS(asv_t, distance = "bray", k = 2, try = 20, trace = FALSE)

if(nmds$converged) {
  cat("NMDS converged successfully. Stress =", round(nmds$stress, 3), "\n")
  
  # Extract scores
  if("sites" %in% names(nmds)) {
    nmds_scores <- nmds$points
  } else {
    nmds_scores <- scores(nmds, display = "sites")
  }
  
  if(is.null(dim(nmds_scores)) || ncol(nmds_scores) != 2) {
    nmds_scores <- scores(nmds)
    if(is.null(dim(nmds_scores))) {
      nmds_scores <- nmds$points
    }
  }
  
  nmds_df <- data.frame(
    SampleID = rownames(nmds_scores),
    NMDS1 = nmds_scores[, 1],
    NMDS2 = nmds_scores[, 2],
    Habitat = coord_filtered$Habitat[match(rownames(nmds_scores), coord_filtered$SampleID)]
  )
  
  p5 <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Habitat)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = SampleID), vjust = -0.7, size = 2, angle = 45) +
    scale_color_manual(values = habitat_colors) +
    stat_ellipse(level = 0.68, linetype = 2, linewidth = 1) +
    theme_minimal() +
    labs(title = "NMDS Ordination: 18S rRNA Microbial Communities",
         subtitle = paste("Stress =", round(nmds$stress, 3), "| Based on Bray-Curtis dissimilarity"),
         x = "NMDS1", y = "NMDS2") +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
  
  plot_list[[length(plot_list) + 1]] <- p5
} else {
  cat("Warning: NMDS did not converge. Skipping NMDS plot.\n")
}

# 6. Hierarchical clustering with habitat colors
hc <- hclust(comm_dist, method = "average")
plot(hc, main = "Hierarchical Clustering: 18S rRNA Microbial Communities", 
     sub = paste("Bray-Curtis dissimilarity |", ncol(asv_t), "ASVs |", nrow(asv_t), "samples"),
     xlab = "", cex = 0.7)

# Add habitat color bar
habitat_colors_vec <- habitat_colors[coord_filtered$Habitat[match(hc$labels, coord_filtered$SampleID)]]
names(habitat_colors_vec) <- hc$labels
legend("topright", legend = names(habitat_colors), fill = habitat_colors, 
       title = "Habitat", cex = 0.8)

# Display all plots
if(length(plot_list) == 3) {
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 2, nrow = 2)
} else if(length(plot_list) >= 4) {
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], ncol = 2, nrow = 2)
  if(length(plot_list) >= 5) plot_list[[5]]  # NMDS plot separately
} else {
  grid.arrange(plot_list[[1]], plot_list[[2]], ncol = 2)
}

# SAVE ALL PLOTS TO FILES
cat("\nSaving plots to files...\n")

# Save individual plots
ggsave("18S_01_sample_locations_final.png", plot_list[[1]], width = 10, height = 8, dpi = 300)
ggsave("18S_02_distance_decay_combined_final.png", plot_list[[2]], width = 10, height = 8, dpi = 300)

if(length(plot_list) >= 3) {
  ggsave("18S_03_distance_decay_ponds_final.png", plot_list[[3]], width = 10, height = 8, dpi = 300)
}

if(length(plot_list) >= 4) {
  ggsave("18S_04_distance_decay_mats_final.png", plot_list[[4]], width = 10, height = 8, dpi = 300)
}

if(length(plot_list) >= 5) {
  ggsave("18S_05_NMDS_ordination_final.png", plot_list[[5]], width = 10, height = 8, dpi = 300)
}

# Save combined plot panel - matching the 16S layout
if(length(plot_list) >= 4) {
  combined_plot <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], 
                                ncol = 2, nrow = 2)
  ggsave("18S_combined_analysis_panel_final.png", combined_plot, width = 16, height = 12, dpi = 300)
} else if(length(plot_list) == 3) {
  combined_plot <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 2, nrow = 2)
  ggsave("18S_combined_analysis_panel_final.png", combined_plot, width = 16, height = 12, dpi = 300)
}

# Save the hierarchical clustering plot
png("18S_07_hierarchical_clustering_final.png", width = 12, height = 8, units = "in", res = 300)
plot(hc, main = "Hierarchical Clustering: 18S rRNA Microbial Communities", 
     sub = paste("Bray-Curtis dissimilarity |", ncol(asv_t), "ASVs |", nrow(asv_t), "samples"),
     xlab = "", cex = 0.7)
habitat_colors_vec <- habitat_colors[coord_filtered$Habitat[match(hc$labels, coord_filtered$SampleID)]]
legend("topright", legend = names(habitat_colors), fill = habitat_colors, 
       title = "Habitat", cex = 0.8)
dev.off()

cat("Plots saved:\n")
cat("- 18S_01_sample_locations_final.png\n")
cat("- 18S_02_distance_decay_combined_final.png\n")
if(length(plot_list) >= 3) cat("- 18S_03_distance_decay_ponds_final.png\n")
if(length(plot_list) >= 4) cat("- 18S_04_distance_decay_mats_final.png\n")
if(length(plot_list) >= 5) cat("- 18S_05_NMDS_ordination_final.png\n")
cat("- 18S_combined_analysis_panel_final.png\n")
cat("- 18S_07_hierarchical_clustering_final.png\n")

# STATISTICAL SUMMARY
cat("\n", rep("=", 60), "\n")
cat("18S rRNA STATISTICAL SUMMARY\n")
cat(rep("=", 60), "\n")

cat("\nDataset Overview:")
cat("\n- Total samples analyzed:", nrow(coord_filtered))
cat("\n- Total 18S ASVs retained:", ncol(asv_t))
cat("\n- Habitat distribution:")
for(hab in unique(coord_filtered$Habitat)) {
  cat("\n  -", hab, ":", sum(coord_filtered$Habitat == hab), "samples")
}

cat("\n\nGeographic spread:")
cat("\n- Latitude range:", round(range(coord_filtered$Latitude), 4))
cat("\n- Longitude range:", round(range(coord_filtered$Longitude), 4))
cat("\n- Max geographic distance:", round(max(geo_vec), 2), "km")

cat("\n\nCommunity dissimilarity:")
cat("\n- Bray-Curtis range:", round(range(comm_vec), 3))
cat("\n- Mean dissimilarity:", round(mean(comm_vec), 3))

cat("\n\nMantel Test Results:")
for(i in 1:nrow(mantel_results)) {
  cat("\n", mantel_results$Analysis[i], ":")
  cat("\n  - Correlation (r):", mantel_results$Mantel_r[i])
  cat("\n  - P-value:", mantel_results$P_value[i])
  cat("\n  - Significant:", ifelse(mantel_results$Significant[i], "YES", "NO"))
  cat("\n  - Sample size:", mantel_results$N_samples[i])
}

# Interpretation guide
cat("\n", rep("=", 60), "\n")
cat("INTERPRETATION GUIDE\n")
cat(rep("=", 60), "\n")
cat("\n- Positive Mantel r: Communities become more dissimilar with distance")
cat("\n- Negative Mantel r: Communities become more similar with distance (rare)")
cat("\n- r ≈ 0: No relationship between distance and community similarity")
cat("\n- Significant p-value (< 0.05): Pattern is not due to chance")
cat("\n- NMDS stress < 0.2: Good representation of community relationships")

cat("\n\nEcological implications:")
if(combined_mantel$result$Significant) {
  cat("\n- Significant overall distance-decay suggests limited dispersal")
} else {
  cat("\n- No significant overall distance-decay suggests high dispersal or environmental filtering")
}

if(exists("pond_mantel") && exists("mat_mantel")) {
  if(pond_mantel$result$Significant && !mat_mantel$result$Significant) {
    cat("\n- Ponds show distance-decay but mats don't: different connectivity mechanisms")
  } else if(!pond_mantel$result$Significant && mat_mantel$result$Significant) {
    cat("\n- Mats show distance-decay but ponds don't: aquatic connectivity in ponds")
  } else if(pond_mantel$result$Significant && mat_mantel$result$Significant) {
    cat("\n- Both habitats show distance-decay: limited dispersal in both")
  }
}

cat("\n\n18S rRNA Analysis completed successfully!\n")

# Store results for further analysis
analysis_results_18S <- list(
  results_table = mantel_results,
  coordinate_data = coord_filtered,
  community_data = asv_t,
  distance_matrices = list(geographic = geo_dist_obj, community = comm_dist),
  ordination = if(exists("nmds")) nmds else NULL
)

# You can access results using:
# analysis_results_18S$results_table
# analysis_results_18S$coordinate_data
# etc.