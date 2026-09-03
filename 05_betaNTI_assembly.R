# Beta Nearest Taxon Index (betaNTI) Analysis for Bacterial Communities
# Analysis of phylogenetic community assembly across Pond, Lake, and Desiccated environments
# TUFTS EXCLUDED — filament samples (filament 1-5) removed from both 16S and 18S

# Load required libraries
library(picante)
library(iCAMP)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)

# Set working directory (adjust as needed)
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

cat("Starting betaNTI analysis for bacterial communities across environments\n")

# ===============================================
# DATA LOADING AND PREPARATION
# ===============================================

###16S
# Load rarefied ASV table (samples as rows, ASVs as columns)
cat("Loading rarefied ASV data...\n")
bac_asv <- read.csv("ASV_table_rare_notufts.csv", check.names = FALSE, row.names = 1)
# Transpose — samples as rows, ASVs as columns
bac_asv <- as.data.frame(t(bac_asv))

# Verify orientation
cat("Rows (should be samples):", rownames(bac_asv)[1:3], "\n")
cat("Cols (should be ESV IDs):", colnames(bac_asv)[1:3], "\n")
cat("Dimensions:", nrow(bac_asv), "samples x", ncol(bac_asv), "ASVs\n")

# Load metadata with environment information
cat("Loading metadata...\n")
metadata <- read.delim("Drymats_metadata_final.txt", check.names = FALSE)

# Check what the first column is called
head(metadata)

# Set sample names as row names
rownames(metadata) <- metadata$Sample  # adjust if column is named differently

print(table(metadata$Environment))

# Ensure sample names match between ASV table and metadata
common_samples <- intersect(rownames(bac_asv), rownames(metadata))
cat("Common samples found:", length(common_samples), "\n")

if(length(common_samples) == 0) {
  stop("No common samples found between ASV table and metadata. Check sample naming.")
}

# Filter to common samples
bac_asv <- bac_asv[common_samples, ]
metadata <- metadata[common_samples, ]

# Filter for target environments
target_envs <- c("Pond", "Lake", "Desiccated")
env_samples <- metadata$Environment %in% target_envs

if(sum(env_samples) == 0) {
  stop("No samples found for target environments. Check environment column.")
}

bac_asv <- bac_asv[env_samples, ]
metadata <- metadata[env_samples, ]

cat("Final dataset:", nrow(bac_asv), "samples,", ncol(bac_asv), "ASVs\n")
cat("Environment distribution:\n")
print(table(metadata$Environment))

# ===============================================
# PHYLOGENETIC TREE PREPARATION
# ===============================================

cat("Loading and preparing phylogenetic tree...\n")
tree <- read.tree("repset_aln.tre")

# Ensure ASV names match between table and tree
# Add quotes around ASV names if needed for tree matching
asv_names_tree <- paste0("'", colnames(bac_asv), "'")
tree_tips_available <- intersect(asv_names_tree, tree$tip.label)

if(length(tree_tips_available) < 10) {
  # Try without quotes
  tree_tips_available <- intersect(colnames(bac_asv), tree$tip.label)
  matching_asvs <- colnames(bac_asv)[colnames(bac_asv) %in% tree$tip.label]
} else {
  # Use quoted names
  matching_asvs <- gsub("'", "", tree_tips_available)
  colnames(bac_asv) <- paste0("'", colnames(bac_asv), "'")
}

cat("ASVs matching tree:", length(tree_tips_available), "\n")

if(length(tree_tips_available) < 10) {
  stop("Too few ASVs match the phylogenetic tree. Check ASV naming convention.")
}

# Prune tree and ASV table to matching taxa
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tree_tips_available])
bac_asv_pruned <- bac_asv[, colnames(bac_asv) %in% tree_pruned$tip.label]

# Reorder ASV table to match tree tip order
bac_asv_pruned <- bac_asv_pruned[, tree_pruned$tip.label]

cat("Final phylogenetic dataset:", ncol(bac_asv_pruned), "ASVs\n")

# Get the number of tips (ASVs)
n_tips <- length(tree_pruned$tip.label)
cat("Number of ASVs in tree:", n_tips, "\n")

# Calculate all pairwise distances (tips + internal nodes)
all_distances <- dist.nodes(tree_pruned)
all_dist_matrix <- as.matrix(all_distances)

# Extract only tip-to-tip distances (first n_tips rows and columns)
phy_dist <- all_dist_matrix[1:n_tips, 1:n_tips]

# Set proper row and column names
rownames(phy_dist) <- tree_pruned$tip.label
colnames(phy_dist) <- tree_pruned$tip.label

# Verify dimensions and content
cat("Phylogenetic distance matrix dimensions:", dim(phy_dist), "\n")
cat("Distance range:", range(phy_dist), "\n")
cat("Sample of distances:\n")
print(phy_dist[1:3, 1:3])

saveRDS(bac_asv_pruned, "bac_asv_pruned_notufts.rds")
saveRDS(phy_dist,       "phy_dist_notufts.rds")
cat("Objects saved. Ready to transfer to HPC.\n")

###VERY IMPORTANT: run this on the server. It is too slow on R studio
# ===============================================
# BETA NTI CALCULATION
# ===============================================

cat("Calculating beta Nearest Taxon Index (betaNTI)...\n")
cat("This may take several minutes...\n")

# Calculate betaNTI using iCAMP package
# betaNTI > +2: heterogeneous selection (deterministic)
# betaNTI < -2: homogeneous selection (deterministic) 
# betaNTI between -2 and +2: stochastic processes
betaNTI_result <- bNTIn.p(
  comm = bac_asv_pruned,
  dis = phy_dist,
  nworker = 4,
  memo.size.GB = 50,
  weighted = TRUE,
  exclude.consp = FALSE,
  rand = 1000,
  output.bMNTD = FALSE,
  sig.index = "bNTI",
  unit.sum = NULL,
  correct.special = FALSE,
  detail.null = FALSE,
  special.method = "MNTD",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975
)

cat("betaNTI calculation completed\n")

# Extract the betaNTI matrix and save as CSV
betaNTI_matrix <- as.matrix(betaNTI_result$index)
write.csv(betaNTI_matrix, "betaNTI_matrix_notufts.csv") #saved in the server and then transferred to my computer

# ===============================================
# DATA PROCESSING FOR VISUALIZATION
# ===============================================

cat("Processing betaNTI results for visualization...\n")
# Read betaNTI matrix from CSV file
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
betaNTI_matrix <- read.csv("betaNTI_matrix_notufts.csv", row.names = 1, check.names = FALSE)

# Extract betaNTI matrix
betaNTI_matrix <- as.matrix(betaNTI_matrix)

# Convert to long format for analysis
betaNTI_long <- expand.grid(
  Sample1 = rownames(betaNTI_matrix),
  Sample2 = colnames(betaNTI_matrix),
  stringsAsFactors = FALSE
)

# Add betaNTI values
betaNTI_long$betaNTI <- as.vector(betaNTI_matrix)

# Remove self-comparisons and duplicate pairs
betaNTI_long <- betaNTI_long[betaNTI_long$Sample1 != betaNTI_long$Sample2, ]

# Remove upper triangle duplicates
sample_pairs <- paste(pmin(betaNTI_long$Sample1, betaNTI_long$Sample2),
                      pmax(betaNTI_long$Sample1, betaNTI_long$Sample2), sep = "_")
betaNTI_long <- betaNTI_long[!duplicated(sample_pairs), ]

# Add environment information
betaNTI_long$Env1 <- metadata[betaNTI_long$Sample1, "Environment"]
betaNTI_long$Env2 <- metadata[betaNTI_long$Sample2, "Environment"]

# Create comparison categories
betaNTI_long$Comparison_Type <- ifelse(
  betaNTI_long$Env1 == betaNTI_long$Env2,
  paste("Within", betaNTI_long$Env1),
  "Between Environments"
)

# Focus on within-environment comparisons for the desired plot
within_env <- betaNTI_long[betaNTI_long$Env1 == betaNTI_long$Env2, ]
within_env$Environment <- within_env$Env1

cat("Data processing completed\n")
cat("Within-environment comparisons:\n")
print(table(within_env$Environment))

# ===============================================
# STATISTICAL ANALYSIS
# ===============================================

cat("Performing statistical analysis...\n")

# Test for differences between environments
if(nrow(within_env) > 0) {
  # ANOVA test
  aov_result <- aov(betaNTI ~ Environment, data = within_env)
  aov_summary <- summary(aov_result)
  
  cat("ANOVA Results:\n")
  print(aov_summary)
  
  # Tukey HSD for pairwise comparisons if significant
  if(aov_summary[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_result <- TukeyHSD(aov_result)
    cat("Tukey HSD Results:\n")
    print(tukey_result)
  }
  
  # Summary statistics by environment
  summary_stats <- within_env %>%
    dplyr::group_by(Environment) %>%
    dplyr::summarise(
      n = n(),
      mean_betaNTI = mean(betaNTI, na.rm = TRUE),
      sd_betaNTI = sd(betaNTI, na.rm = TRUE),
      median_betaNTI = median(betaNTI, na.rm = TRUE),
      prop_deterministic = mean(abs(betaNTI) > 2, na.rm = TRUE),
      prop_stochastic = mean(abs(betaNTI) <= 2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Summary statistics by environment:\n")
  print(summary_stats)
} else {
  cat("Warning: No within-environment comparisons found\n")
}

# ===============================================
# VISUALIZATION
# ===============================================

cat("Creating visualization...\n")

# Set factor levels for consistent ordering
within_env$Environment <- factor(within_env$Environment, 
                                 levels = c("Desiccated", "Pond", "Lake"))

# Create the plot matching your desired style
p <- ggplot(within_env, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  
  # Add significance threshold lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", size = 0.5) +
  
  # Labels and theme
  labs(
    x = "Environment",
    y = "Nearest Taxon Index (NTI)",
    title = "NTI Across Different Environments"
  ) +
  
  # Set colors to match your plot
  scale_color_manual(values = c("Desiccated" = "#F8766D", 
                                "Pond" = "#619CFF", 
                                "Lake" = "#00BA38")) +
  
  # Theme styling
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )

# Display the plot
print(p)

# Save the plot
ggsave("betaNTI_across_environments_notufts.pdf", plot = p, width = 8, height = 6, dpi = 300)
ggsave("betaNTI_across_environments_notufts.png", plot = p, width = 8, height = 6, dpi = 300)

cat("Plot saved as betaNTI_across_environments_notufts.pdf and .png\n")

# ===============================================
# SAVE RESULTS
# ===============================================

cat("Saving results...\n")

# Save processed data
write.csv(within_env, "betaNTI_within_environments_notufts.csv", row.names = FALSE)
write.csv(betaNTI_long, "betaNTI_all_comparisons_notufts.csv", row.names = FALSE)

# Save summary statistics
if(exists("summary_stats")) {
  write.csv(summary_stats, "betaNTI_summary_statistics_notufts.csv", row.names = FALSE)
}

# Save betaNTI matrix
write.csv(betaNTI_matrix, "betaNTI_matrix_notufts.csv")

cat("Analysis completed successfully!\n")
cat("Files created:\n")
cat("- betaNTI_across_environments_notufts.pdf/.png: Main plot\n")
cat("- betaNTI_within_environments_notufts.csv: Within-environment comparisons\n")
cat("- betaNTI_all_comparisons_notufts.csv: All pairwise comparisons\n")
cat("- betaNTI_summary_statistics_notufts.csv: Summary statistics\n")
cat("- betaNTI_matrix_notufts.csv: Full betaNTI matrix\n")

# ===============================================
# INTERPRETATION GUIDE
# ===============================================

cat("\n", rep("=", 60), "\n")
cat("INTERPRETATION GUIDE\n")
cat(rep("=", 60), "\n")
cat("betaNTI > +2: Heterogeneous selection (deterministic processes)\n")
cat("betaNTI < -2: Homogeneous selection (deterministic processes)\n") 
cat("betaNTI between -2 and +2: Stochastic processes dominate\n")
cat("\nDashed lines at ±2 indicate the threshold for deterministic vs stochastic processes\n")

###18S
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
cat("Loading rarefied ASV data...\n")
euk_asv <- read.csv("ASV_table_rare_notufts.csv", check.names = FALSE, row.names = 1)

# Transpose — samples as rows, ASVs as columns
euk_asv <- as.data.frame(t(euk_asv))
cat("Rows (should be samples):", rownames(euk_asv)[1:3], "\n")
cat("Cols (should be ESV IDs):", colnames(euk_asv)[1:3], "\n")
cat("Dimensions:", nrow(euk_asv), "samples x", ncol(euk_asv), "ASVs\n")

# Load metadata
cat("Loading metadata...\n")
metadata <- read.delim("Drymats_metadata_final.txt", check.names = FALSE)

# Fix row names — set sample names as row names
rownames(metadata) <- metadata$Sample
print(table(metadata$Environment))

print(table(metadata$Environment))

# Match samples
common_samples <- intersect(rownames(euk_asv), rownames(metadata))
cat("Common samples found:", length(common_samples), "\n")

if(length(common_samples) == 0) {
  stop("No common samples found. Check sample naming.")
}

# Filter to common samples
euk_asv <- euk_asv[common_samples, ]
metadata <- metadata[common_samples, ]

# Filter for target environments
target_envs <- c("Pond", "Lake", "Desiccated")
env_samples <- metadata$Environment %in% target_envs

if(sum(env_samples) == 0) {
  stop("No samples found for target environments. Check environment column.")
}

euk_asv <- euk_asv[env_samples, ]
metadata <- metadata[env_samples, ]

cat("Final dataset:", nrow(euk_asv), "samples,", ncol(euk_asv), "ASVs\n")
cat("Environment distribution:\n")
print(table(metadata$Environment))

# ===============================================
# PHYLOGENETIC TREE PREPARATION
# ===============================================

cat("Loading and preparing phylogenetic tree...\n")
tree <- read.tree("repset_aln.tre")

# Ensure ASV names match between table and tree
# Add quotes around ASV names if needed for tree matching
asv_names_tree <- paste0("'", colnames(euk_asv), "'")
tree_tips_available <- intersect(asv_names_tree, tree$tip.label)

if(length(tree_tips_available) < 10) {
  # Try without quotes
  tree_tips_available <- intersect(colnames(euk_asv), tree$tip.label)
  matching_asvs <- colnames(euk_asv)[colnames(euk_asv) %in% tree$tip.label]
} else {
  # Use quoted names
  matching_asvs <- gsub("'", "", tree_tips_available)
  colnames(euk_asv) <- paste0("'", colnames(euk_asv), "'")
}

cat("ASVs matching tree:", length(tree_tips_available), "\n")

if(length(tree_tips_available) < 10) {
  stop("Too few ASVs match the phylogenetic tree. Check ASV naming convention.")
}

# Prune tree and ASV table to matching taxa
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% tree_tips_available])
euk_asv_pruned <- euk_asv[, colnames(euk_asv) %in% tree_pruned$tip.label]

# Reorder ASV table to match tree tip order
euk_asv_pruned <- euk_asv_pruned[, tree_pruned$tip.label]

cat("Final phylogenetic dataset:", ncol(euk_asv_pruned), "ASVs\n")

# Calculate phylogenetic distance matrix
cat("Calculating phylogenetic distances...\n")
phy_dist <- cophenetic(tree_pruned)

saveRDS(euk_asv_pruned, "euk_asv_pruned_notufts.rds")
saveRDS(phy_dist,       "phy_dist_18S_notufts.rds")
cat("18S objects saved. Ready to transfer to HPC.\n")

# ===============================================
# BETA NTI CALCULATION
# ===============================================

cat("Calculating beta Nearest Taxon Index (betaNTI)...\n")
cat("This may take several minutes...\n")

# Calculate betaNTI using iCAMP package
# betaNTI > +2: heterogeneous selection (deterministic)
# betaNTI < -2: homogeneous selection (deterministic) 
# betaNTI between -2 and +2: stochastic processes
betaNTI_result <- bNTIn.p(
  comm = euk_asv_pruned,
  dis = phy_dist,
  nworker = 4,
  memo.size.GB = 50,
  weighted = TRUE,
  exclude.consp = FALSE,
  rand = 1000,
  output.bMNTD = FALSE,
  sig.index = "bNTI",
  unit.sum = NULL,
  correct.special = FALSE,
  detail.null = FALSE,
  special.method = "MNTD",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975
)

cat("betaNTI calculation completed\n")

# Extract the betaNTI matrix and save as CSV
betaNTI_matrix <- as.matrix(betaNTI_result$index)
write.csv(betaNTI_matrix, "18SbetaNTI_matrix_notufts.csv") #saved in the server and then transferred to my computer

# ===============================================
# DATA PROCESSING FOR VISUALIZATION
# ===============================================
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
cat("Processing betaNTI results for visualization...\n")
# Read betaNTI matrix from CSV file
betaNTI_matrix <- read.csv("betaNTI_matrix_notufts.csv", row.names = 1, check.names = FALSE)

# Extract betaNTI matrix
betaNTI_matrix <- as.matrix(betaNTI_matrix)

# Convert to long format for analysis
betaNTI_long <- expand.grid(
  Sample1 = rownames(betaNTI_matrix),
  Sample2 = colnames(betaNTI_matrix),
  stringsAsFactors = FALSE
)

# Add betaNTI values
betaNTI_long$betaNTI <- as.vector(betaNTI_matrix)

# Remove self-comparisons and duplicate pairs
betaNTI_long <- betaNTI_long[betaNTI_long$Sample1 != betaNTI_long$Sample2, ]

# Remove upper triangle duplicates
sample_pairs <- paste(pmin(betaNTI_long$Sample1, betaNTI_long$Sample2),
                      pmax(betaNTI_long$Sample1, betaNTI_long$Sample2), sep = "_")
betaNTI_long <- betaNTI_long[!duplicated(sample_pairs), ]

# Add environment information
betaNTI_long$Env1 <- metadata[betaNTI_long$Sample1, "Environment"]
betaNTI_long$Env2 <- metadata[betaNTI_long$Sample2, "Environment"]

# Create comparison categories
betaNTI_long$Comparison_Type <- ifelse(
  betaNTI_long$Env1 == betaNTI_long$Env2,
  paste("Within", betaNTI_long$Env1),
  "Between Environments"
)

# Focus on within-environment comparisons for the desired plot
within_env <- betaNTI_long[betaNTI_long$Env1 == betaNTI_long$Env2, ]
within_env$Environment <- within_env$Env1

cat("Data processing completed\n")
cat("Within-environment comparisons:\n")
print(table(within_env$Environment))

# ===============================================
# STATISTICAL ANALYSIS
# ===============================================

cat("Performing statistical analysis...\n")

# Test for differences between environments
if(nrow(within_env) > 0) {
  # ANOVA test
  aov_result <- aov(betaNTI ~ Environment, data = within_env)
  aov_summary <- summary(aov_result)
  
  cat("ANOVA Results:\n")
  print(aov_summary)
  
  # Tukey HSD for pairwise comparisons if significant
  if(aov_summary[[1]][["Pr(>F)"]][1] < 0.05) {
    tukey_result <- TukeyHSD(aov_result)
    cat("Tukey HSD Results:\n")
    print(tukey_result)
  }
  
  # Summary statistics by environment
  summary_stats <- within_env %>%
    dplyr::group_by(Environment) %>%
    dplyr::summarise(
      n = n(),
      mean_betaNTI = mean(betaNTI, na.rm = TRUE),
      sd_betaNTI = sd(betaNTI, na.rm = TRUE),
      median_betaNTI = median(betaNTI, na.rm = TRUE),
      prop_deterministic = mean(abs(betaNTI) > 2, na.rm = TRUE),
      prop_stochastic = mean(abs(betaNTI) <= 2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Summary statistics by environment:\n")
  print(summary_stats)
} else {
  cat("Warning: No within-environment comparisons found\n")
}

# ===============================================
# VISUALIZATION
# ===============================================

cat("Creating visualization...\n")

# Set factor levels for consistent ordering
within_env$Environment <- factor(within_env$Environment, 
                                 levels = c("Desiccated", "Pond", "Lake"))

# Create the plot matching your desired style
p <- ggplot(within_env, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  
  # Add significance threshold lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", size = 0.5) +
  
  # Labels and theme
  labs(
    x = "Environment",
    y = "Nearest Taxon Index (NTI)",
    title = "NTI Across Different Environments"
  ) +
  
  # Set colors to match your plot
  scale_color_manual(values = c("Desiccated" = "#F8766D", 
                                "Pond" = "#619CFF", 
                                "Lake" = "#00BA38")) +
  
  # Theme styling
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  )

# Display the plot
print(p)

# Save the plot
ggsave("betaNTI_across_environments_18S_notufts.pdf", plot = p, width = 8, height = 6, dpi = 300)
ggsave("betaNTI_across_environments_18S_notufts.png", plot = p, width = 8, height = 6, dpi = 300)

cat("Plot saved as betaNTI_across_environments_18S_notufts.pdf and .png\n")

# ===============================================
# SAVE RESULTS
# ===============================================

cat("Saving results...\n")

# Save processed data
write.csv(within_env, "betaNTI_within_environments_notufts.csv", row.names = FALSE)
write.csv(betaNTI_long, "betaNTI_all_comparisons_notufts.csv", row.names = FALSE)

# Save summary statistics
if(exists("summary_stats")) {
  write.csv(summary_stats, "betaNTI_summary_statistics_notufts.csv", row.names = FALSE)
}

# Save betaNTI matrix
write.csv(betaNTI_matrix, "18SbetaNTI_matrix_notufts.csv")

cat("Analysis completed successfully!\n")
cat("Files created:\n")
cat("- betaNTI_across_environments_18S_notufts.pdf/.png: Main plot\n")
cat("- betaNTI_within_environments_notufts.csv: Within-environment comparisons\n")
cat("- betaNTI_all_comparisons_notufts.csv: All pairwise comparisons\n")
cat("- betaNTI_summary_statistics_notufts.csv: Summary statistics\n")
cat("- 18SbetaNTI_matrix_notufts.csv: Full betaNTI matrix\n")

# ===============================================
# INTERPRETATION GUIDE
# ===============================================

cat("\n", rep("=", 60), "\n")
cat("INTERPRETATION GUIDE\n")
cat(rep("=", 60), "\n")
cat("betaNTI > +2: Heterogeneous selection (deterministic processes)\n")
cat("betaNTI < -2: Homogeneous selection (deterministic processes)\n") 
cat("betaNTI between -2 and +2: Stochastic processes dominate\n")
cat("\nDashed lines at ±2 indicate the threshold for deterministic vs stochastic processes\n")


### Combined 16S and 18S betaNTI Visualization
# Creates side-by-side plots with aligned y-axes for comparison

# Load required libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# ===============================================
# LOAD AND PREPARE DATA FROM BOTH ANALYSES
# ===============================================

# Load 16S betaNTI results
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
within_env_16s <- read.csv("betaNTI_within_environments_notufts.csv", stringsAsFactors = FALSE)

# Load 18S betaNTI results
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
within_env_18s <- read.csv("betaNTI_within_environments_notufts.csv", stringsAsFactors = FALSE)

# Set factor levels for consistent ordering
within_env_16s$Environment <- factor(within_env_16s$Environment, 
                                     levels = c("Desiccated", "Pond", "Lake"))
within_env_18s$Environment <- factor(within_env_18s$Environment, 
                                     levels = c("Desiccated", "Pond", "Lake"))

# ===============================================
# DETERMINE COMMON Y-AXIS RANGE
# ===============================================

# Find the range of betaNTI values across both datasets
y_min <- min(c(within_env_16s$betaNTI, within_env_18s$betaNTI), na.rm = TRUE)
y_max <- max(c(within_env_16s$betaNTI, within_env_18s$betaNTI), na.rm = TRUE)

# Add some padding and ensure -2 and +2 lines are visible
y_range <- c(min(y_min, -3), max(y_max, 3))

cat("Combined y-axis range:", y_range, "\n")

# ===============================================
# CREATE 16S PLOT
# ===============================================

plot_16s <- ggplot(within_env_16s, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  
  # Add significance threshold lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", size = 0.5) +
  
  # Set common y-axis range
  ylim(y_range) +
  
  # Labels and theme
  labs(
    x = "Environment",
    y = expression(beta*NTI),
    title = "16S Bacterial Communities"
  ) +
  
  # Set colors for 16S
  scale_color_manual(values = c("Desiccated" = "#F8766D", 
                                "Pond" = "#619CFF", 
                                "Lake" = "#00BA38")) +
  
  # Theme styling
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    plot.margin = margin(10, 10, 10, 20)
  )

# ===============================================
# CREATE 18S PLOT
# ===============================================

plot_18s <- ggplot(within_env_18s, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  
  # Add significance threshold lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", size = 0.5) +
  
  # Set common y-axis range
  ylim(y_range) +
  
  # Labels and theme
  labs(
    x = "Environment",
    y = expression(beta*NTI),
    title = "18S Eukaryotic Communities"
  ) +
  
  # Set colors for 18S (different from 16S)
  scale_color_manual(values = c("Desiccated" = "#E31A1C", 
                                "Pond" = "#1F78B4", 
                                "Lake" = "#33A02C")) +
  
  # Theme styling
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    plot.margin = margin(10, 20, 10, 10)
  )

# ===============================================
# COMBINE PLOTS SIDE BY SIDE
# ===============================================
library(grid)
library(gridExtra)
# Create combined plot
combined_plot <- grid.arrange(
  plot_16s, 
  plot_18s, 
  ncol = 2, 
  # Use expression() to force R to render the beta symbol mathematically
  top = grid::textGrob(expression(bold(beta*NTI~"Across Different Environments")), 
                       gp = grid::gpar(fontsize = 16))
)

# ===============================================
# SAVE COMBINED PLOT
# ===============================================

# Save to output directory
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/betaNTI_no_tufts/")

# Save as PDF
ggsave("combined_16S_18S_betaNTI_notufts.pdf", plot = combined_plot, 
       width = 12, height = 6, dpi = 300)

# Save as PNG
ggsave("combined_16S_18S_betaNTI_notufts.png", plot = combined_plot, 
       width = 12, height = 6, dpi = 300)

# Alternative method using cowplot for better control
library(cowplot)

combined_plot_cowplot <- plot_grid(plot_16s, plot_18s, 
                                   labels = c("A", "B"), 
                                   label_size = 16,
                                   ncol = 2, 
                                   align = "hv")

# Add main title
title <- ggdraw() + 
  draw_label(expression(bold(beta*NTI~"Across Different Environments")), 
             size = 16)

final_plot <- plot_grid(title, combined_plot_cowplot, 
                        ncol = 1, rel_heights = c(0.1, 1))

# Save cowplot version
ggsave("combined_16S_18S_betaNTI_cowplot_notufts.pdf", plot = final_plot, 
       width = 12, height = 6, dpi = 300)

ggsave("combined_16S_18S_betaNTI_cowplot_notufts.png", plot = final_plot, 
       width = 12, height = 6, dpi = 300)

# ===============================================
# PRINT SUMMARY STATISTICS
# ===============================================

cat("16S betaNTI Summary:\n")
summary_16s <- within_env_16s %>%
  dplyr::group_by(Environment) %>%
  dplyr::summarise(                  # <-- Explicitly tell R to use dplyr's summarise
    n = dplyr::n(),
    mean = round(mean(betaNTI, na.rm = TRUE), 3),
    median = round(median(betaNTI, na.rm = TRUE), 3),
    sd = round(sd(betaNTI, na.rm = TRUE), 3),
    min = round(min(betaNTI, na.rm = TRUE), 3),
    max = round(max(betaNTI, na.rm = TRUE), 3)
  )
print(summary_16s)

cat("\n18S betaNTI Summary:\n")
summary_18s <- within_env_18s %>%
  dplyr::group_by(Environment) %>%
  dplyr::summarise(                  # <-- Explicitly tell R to use dplyr's summarise
    n = dplyr::n(),
    mean = round(mean(betaNTI, na.rm = TRUE), 3),
    median = round(median(betaNTI, na.rm = TRUE), 3),
    sd = round(sd(betaNTI, na.rm = TRUE), 3),
    min = round(min(betaNTI, na.rm = TRUE), 3),
    max = round(max(betaNTI, na.rm = TRUE), 3)
  )
print(summary_18s)

cat("\nPlots saved as:\n")
cat("- combined_16S_18S_betaNTI_notufts.pdf/.png (gridExtra version)\n")
cat("- combined_16S_18S_betaNTI_cowplot_notufts.pdf/.png (cowplot version with labels)\n")

# Display the plot
print(final_plot)

# =============================================================================
# FIGURE 6 (FINAL): Combined βNTI plot — Bacteria (A) and Eukarya (B),
# square panels, shared y-axis, shared habitat colour palette, A)/B) tags
# positioned close to each panel title.
# =============================================================================

library(ggplot2)
library(cowplot)
library(dplyr)

# -----------------------------------------------------------------------
# STEP 1: Reload already-computed betaNTI results (fast — no need to
# re-run the slow bNTIn.p() calculation)
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")
within_env_16s <- read.csv("betaNTI_within_environments_notufts.csv", stringsAsFactors = FALSE)

setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")
within_env_18s <- read.csv("betaNTI_within_environments_notufts.csv", stringsAsFactors = FALSE)

# -----------------------------------------------------------------------
# Shared habitat colour palette (matches Figure 2 and other manuscript figures)
# -----------------------------------------------------------------------
env_colours <- c("Desiccated" = "#E31A1C",
                 "Pond"       = "#1F78B4",
                 "Lake"       = "#33A02C")

habitat_order <- c("Desiccated", "Pond", "Lake")

within_env_16s$Environment <- factor(within_env_16s$Environment, levels = habitat_order)
within_env_18s$Environment <- factor(within_env_18s$Environment, levels = habitat_order)

# -----------------------------------------------------------------------
# Shared y-axis range across both panels
# -----------------------------------------------------------------------
y_min <- min(c(within_env_16s$betaNTI, within_env_18s$betaNTI), na.rm = TRUE)
y_max <- max(c(within_env_16s$betaNTI, within_env_18s$betaNTI), na.rm = TRUE)
y_range <- c(min(y_min, -3), max(y_max, 3))

# -----------------------------------------------------------------------
# Panel A — Bacterial (16S rRNA gene)
# -----------------------------------------------------------------------
plot_16s <- ggplot(within_env_16s, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  geom_hline(yintercept = 2,  linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", linewidth = 0.5) +
  ylim(y_range) +
  labs(x = NULL, y = expression(bold(beta*NTI)),
       title = "Bacterial (16S rRNA gene) communities") +
  scale_color_manual(values = env_colours) +
  theme_bw() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.title.y  = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(size = 10, face = "bold"),
        axis.text.y   = element_text(size = 10, face = "bold"),
        legend.position = "none",
        aspect.ratio  = 1,
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        plot.margin = margin(2, 10, 10, 20))

# -----------------------------------------------------------------------
# Panel B — Eukaryotic (18S rRNA gene)
# -----------------------------------------------------------------------
plot_18s <- ggplot(within_env_18s, aes(x = Environment, y = betaNTI, color = Environment)) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = NA, lwd = 0.7) +
  geom_hline(yintercept = 2,  linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "black", linewidth = 0.5) +
  ylim(y_range) +
  labs(x = NULL, y = expression(bold(beta*NTI)),
       title = "Eukaryotic (18S rRNA gene) communities") +
  scale_color_manual(values = env_colours) +
  theme_bw() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        axis.title.y  = element_text(size = 12, face = "bold"),
        axis.text.x   = element_text(size = 10, face = "bold"),
        axis.text.y   = element_text(size = 10, face = "bold"),
        legend.position = "none",
        aspect.ratio  = 1,
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95"),
        plot.margin = margin(2, 20, 10, 10))

# -----------------------------------------------------------------------
# Combine with A)/B) panel labels, positioned close to each panel title
# -----------------------------------------------------------------------
combined_plot_cowplot <- plot_grid(plot_16s, plot_18s,
                                   labels = c("A)", "B)"),
                                   label_size = 16,
                                   label_fontface = "bold",
                                   label_x = 0.02,
                                   label_y = 0.90,
                                   hjust = 0,
                                   vjust = 1,
                                   ncol = 2)   # <-- removed align = "hv"

# -----------------------------------------------------------------------
# Save — PDF and high-res JPEG
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/betaNTI_no_tufts/")

ggsave("Figure6_Combined_betaNTI.pdf", plot = combined_plot_cowplot,
       width = 12, height = 6, dpi = 300)   # <-- height reverted to 6

ggsave("Figure6_Combined_betaNTI.jpeg", plot = combined_plot_cowplot,
       width = 12, height = 6, dpi = 300, quality = 100)

cat("Saved: Figure6_Combined_betaNTI.pdf and Figure6_Combined_betaNTI.jpeg\n")

# -----------------------------------------------------------------------
# Print summary statistics for reference
# -----------------------------------------------------------------------
cat("16S betaNTI Summary:\n")
summary_16s <- within_env_16s %>%
  dplyr::group_by(Environment) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = round(mean(betaNTI, na.rm = TRUE), 3),
    median = round(median(betaNTI, na.rm = TRUE), 3),
    sd = round(sd(betaNTI, na.rm = TRUE), 3),
    min = round(min(betaNTI, na.rm = TRUE), 3),
    max = round(max(betaNTI, na.rm = TRUE), 3)
  )
print(summary_16s)

cat("\n18S betaNTI Summary:\n")
summary_18s <- within_env_18s %>%
  dplyr::group_by(Environment) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = round(mean(betaNTI, na.rm = TRUE), 3),
    median = round(median(betaNTI, na.rm = TRUE), 3),
    sd = round(sd(betaNTI, na.rm = TRUE), 3),
    min = round(min(betaNTI, na.rm = TRUE), 3),
    max = round(max(betaNTI, na.rm = TRUE), 3)
  )
print(summary_18s)
