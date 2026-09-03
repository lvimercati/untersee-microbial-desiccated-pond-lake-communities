# Core Microbiome Analysis for Combined Pond-Lake-Desiccated Environments
# Updated version with improvements
# Based on microbiome R package methodology
# TUFTS EXCLUDED — filament samples (filament 1-5) removed from both 16S and 18S

# Set working directory
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

# Define output directory
output_dir <- getwd()

# Create core microbiome analysis directory — notufts version
core_dir <- file.path(output_dir, "core_microbiome_notufts")
if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)

# Load required libraries
library(microbiome)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(mctoolsr)

cat("\n", rep("=", 60), "\n")
cat("CORE MICROBIOME ANALYSIS - TUFTS EXCLUDED\n")
cat(rep("=", 60), "\n")

# Function to prepare data for core analysis
prepare_core_data <- function(data_loaded, metadata_loaded, taxonomy_loaded,
                              environments = c("Lake", "Pond", "Desiccated")) {
  
  cat("Preparing data for core microbiome analysis...\n")
  
  # Debug metadata structure
  cat("Metadata structure:\n")
  cat("Class:", class(metadata_loaded), "\n")
  cat("Dimensions:", dim(metadata_loaded), "\n")
  cat("Column names:", colnames(metadata_loaded), "\n")
  
  # Convert to data frame if needed
  if(!is.data.frame(metadata_loaded)) {
    metadata_loaded <- as.data.frame(metadata_loaded, stringsAsFactors = FALSE)
  }
  
  # The row names of metadata should be the sample names
  cat("Metadata row names (first 5):\n")
  print(head(rownames(metadata_loaded), 5))
  
  cat("Data column names (first 5):\n")
  print(head(colnames(data_loaded), 5))
  
  cat("Environment values:", unique(metadata_loaded$Environment), "\n")
  
  # Filter metadata for target environments
  meta_filtered <- metadata_loaded[metadata_loaded$Environment %in% environments, ]
  cat("Samples in target environments:", nrow(meta_filtered), "\n")
  cat("Environment distribution:\n")
  print(table(meta_filtered$Environment))
  
  # Get common samples between data and metadata
  common_samples <- intersect(colnames(data_loaded), rownames(meta_filtered))
  cat("Common samples found:", length(common_samples), "\n")
  
  if(length(common_samples) == 0) {
    cat("\nDebugging sample name mismatch:\n")
    cat("Data column names (first 10):\n")
    data_names <- colnames(data_loaded)[1:min(10, ncol(data_loaded))]
    for(i in seq_along(data_names)) {
      cat(i, ":", data_names[i], "\n")
    }
    
    cat("\nMetadata row names (first 10):\n")
    meta_names <- rownames(meta_filtered)[1:min(10, nrow(meta_filtered))]
    for(i in seq_along(meta_names)) {
      cat(i, ":", meta_names[i], "\n")
    }
    
    # Check if any data column names contain metadata row names as substrings
    cat("\nChecking for substring matches...\n")
    for(i in 1:min(5, length(meta_names))) {
      matching_cols <- grep(meta_names[i], colnames(data_loaded), value = TRUE)
      if(length(matching_cols) > 0) {
        cat("Metadata:", meta_names[i], "-> Data matches:", paste(matching_cols, collapse = ", "), "\n")
      }
    }
    
    stop("No common samples found between data and metadata")
  }
  
  # Filter data matrices
  count_matrix <- data_loaded[, common_samples]
  metadata_final <- meta_filtered[common_samples, ]
  
  # Remove ASVs with zero counts across all samples
  count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
  
  # Filter taxonomy to match ASVs
  taxonomy_final <- taxonomy_loaded[rownames(count_matrix), ]
  
  cat("Final dataset:\n")
  cat("- ASVs:", nrow(count_matrix), "\n")
  cat("- Samples:", ncol(count_matrix), "\n")
  cat("- Environments:", paste(unique(metadata_final$Environment), collapse = ", "), "\n")
  
  return(list(
    counts = count_matrix,
    metadata = metadata_final,
    taxonomy = taxonomy_final
  ))
}

# Function to create phyloseq object
create_phyloseq_object <- function(core_data) {
  
  cat("Creating phyloseq object...\n")
  
  # Create phyloseq components
  otu_table <- otu_table(core_data$counts, taxa_are_rows = TRUE)
  sample_data <- sample_data(core_data$metadata)
  tax_table <- tax_table(as.matrix(core_data$taxonomy))
  
  # Create phyloseq object
  ps <- phyloseq(otu_table, sample_data, tax_table)
  
  cat("Phyloseq object created with", ntaxa(ps), "taxa and", nsamples(ps), "samples\n")
  
  return(ps)
}

# Function to perform core microbiome analysis
perform_core_analysis <- function(ps, detection_thresholds = NULL, prevalence_thresholds = NULL) {
  
  cat("Performing core microbiome analysis...\n")
  
  # Transform to relative abundances (compositional)
  ps_rel <- microbiome::transform(ps, "compositional")
  
  # Set default thresholds if not provided
  if(is.null(detection_thresholds)) {
    detection_thresholds <- c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                              1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95)/100
  }
  
  if(is.null(prevalence_thresholds)) {
    prevalence_thresholds <- seq(0.05, 1, 0.05)
  }
  
  cat("Detection thresholds:", length(detection_thresholds), "(",
      round(min(detection_thresholds)*100, 3), "% to",
      round(max(detection_thresholds)*100, 1), "%)\n")
  
  # Calculate core microbiome across thresholds
  core_results <- list()
  
  for(i in seq_along(detection_thresholds)) {
    det_thresh <- detection_thresholds[i]
    
    prevalence_data <- c()
    
    for(asv in taxa_names(ps_rel)) {
      asv_abundances <- otu_table(ps_rel)[asv, ]
      n_samples_above_threshold <- sum(asv_abundances >= det_thresh)
      prevalence <- n_samples_above_threshold / nsamples(ps_rel)
      prevalence_data[asv] <- prevalence
    }
    
    core_results[[paste0("det_", round(det_thresh*100, 4), "pct")]] <- prevalence_data
  }
  
  return(list(
    ps_rel = ps_rel,
    core_results = core_results,
    detection_thresholds = detection_thresholds,
    prevalence_thresholds = prevalence_thresholds
  ))
}

# Function to create core microbiome heatmap
create_core_heatmap <- function(core_analysis, top_n_taxa = 15, save_plot = TRUE, 
                                dataset_name = "16S", title_suffix = "Combined Environments") {
  
  cat("Creating core microbiome heatmap...\n")
  
  core_matrix <- do.call(cbind, core_analysis$core_results)
  colnames(core_matrix) <- paste0(round(core_analysis$detection_thresholds*100, 4), "%")
  
  first_col_prevalence <- core_matrix[, 1]
  top_taxa <- names(sort(first_col_prevalence, decreasing = TRUE))[1:min(top_n_taxa, nrow(core_matrix))]
  
  core_matrix_subset <- core_matrix[top_taxa, ]
  
  ps_rel <- core_analysis$ps_rel
  tax_table_df <- as.data.frame(tax_table(ps_rel))
  
  genus_labels <- character(length(top_taxa))
  for(i in seq_along(top_taxa)) {
    asv_id <- top_taxa[i]
    if(asv_id %in% rownames(tax_table_df)) {
      best_label <- asv_id
      
      for(tax_level in 6:1) {
        tax_col <- paste0("taxonomy", tax_level)
        if(tax_col %in% colnames(tax_table_df)) {
          tax_name <- tax_table_df[asv_id, tax_col]
          if(!is.na(tax_name) && tax_name != "" && tax_name != "NA") {
            tax_name_clean <- gsub("_X$", "", tax_name)
            if(tax_name_clean != "") {
              best_label <- paste0(tax_name_clean, " (", asv_id, ")")
              break
            }
          }
        }
      }
      genus_labels[i] <- best_label
    } else {
      genus_labels[i] <- asv_id
    }
  }
  
  rownames(core_matrix_subset) <- genus_labels
  
  core_df <- reshape2::melt(core_matrix_subset)
  colnames(core_df) <- c("Taxon", "Detection_Threshold", "Prevalence")
  
  taxon_avg_prevalence <- aggregate(Prevalence ~ Taxon, core_df, mean)
  taxon_order <- taxon_avg_prevalence$Taxon[order(taxon_avg_prevalence$Prevalence, decreasing = TRUE)]
  
  core_df$Taxon <- factor(core_df$Taxon, levels = rev(taxon_order))
  
  p <- ggplot(core_df, aes(x = Detection_Threshold, y = Taxon, fill = Prevalence)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0.5, name = "Prevalence",
                         labels = function(x) paste0(round(x*100), "%"),
                         guide = guide_colorbar(title.position = "top")) +
    labs(title = paste0("Core Microbiome Analysis - ", dataset_name, " ", title_suffix,
                        " (tufts excluded)"),
         subtitle = "Prevalence of taxa across detection thresholds (ordered by prevalence)",
         x = "Detection Threshold (Relative Abundance %)",
         y = "Taxon (Genus)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "right",
          legend.title = element_text(angle = 0),
          panel.grid = element_blank()) +
    scale_x_discrete(breaks = colnames(core_matrix_subset)[seq(1, ncol(core_matrix_subset), 3)])
  
  if(save_plot) {
    ggsave(file.path(core_dir, paste0("core_microbiome_heatmap_", dataset_name, ".pdf")),
           p, width = 12, height = 8, dpi = 300)
    ggsave(file.path(core_dir, paste0("core_microbiome_heatmap_", dataset_name, ".png")),
           p, width = 12, height = 8, dpi = 300)
    cat("Heatmap saved to:", file.path(core_dir, paste0("core_microbiome_heatmap_", dataset_name, ".pdf/.png")), "\n")
  }
  
  return(p)
}

# Function to identify and summarize core taxa with mean abundance
identify_core_taxa <- function(core_analysis, taxonomy_data,
                               detection_threshold = 0.001, prevalence_threshold = 0.9) {
  
  cat("\nIdentifying core taxa...\n")
  cat("Criteria: >", detection_threshold*100, "% abundance in >", prevalence_threshold*100, "% of samples\n")
  
  det_thresh_col <- which.min(abs(core_analysis$detection_thresholds - detection_threshold))
  actual_threshold <- core_analysis$detection_thresholds[det_thresh_col]
  
  cat("Using detection threshold:", round(actual_threshold*100, 3), "%\n")
  
  prevalence_data <- core_analysis$core_results[[det_thresh_col]]
  
  core_taxa <- names(prevalence_data[prevalence_data >= prevalence_threshold])
  
  cat("Core taxa identified:", length(core_taxa), "\n")
  
  if(length(core_taxa) > 0) {
    core_taxonomy <- taxonomy_data[core_taxa, ]
    
    ps_rel <- core_analysis$ps_rel
    otu_mat <- as.matrix(otu_table(ps_rel))
    
    mean_abundance   <- apply(otu_mat[core_taxa, , drop = FALSE], 1, mean)
    median_abundance <- apply(otu_mat[core_taxa, , drop = FALSE], 1, median)
    max_abundance    <- apply(otu_mat[core_taxa, , drop = FALSE], 1, max)
    
    core_summary <- data.frame(
      ASV              = core_taxa,
      Prevalence       = round(prevalence_data[core_taxa], 3),
      Mean_Abundance   = round(mean_abundance[core_taxa], 4),
      Median_Abundance = round(median_abundance[core_taxa], 4),
      Max_Abundance    = round(max_abundance[core_taxa], 4),
      core_taxonomy,
      stringsAsFactors = FALSE
    )
    
    core_summary <- core_summary[order(core_summary$Prevalence, decreasing = TRUE), ]
    
    cat("\nCore taxa summary:\n")
    print(core_summary)
    
    return(core_summary)
  } else {
    cat("No taxa meet the core criteria\n")
    return(NULL)
  }
}

# Function to create prevalence vs detection threshold curve
create_prevalence_curve <- function(core_analysis, save_plot = TRUE, dataset_name = "16S") {
  
  cat("Creating prevalence curve...\n")
  
  n_taxa_data <- data.frame(
    detection_threshold = core_analysis$detection_thresholds,
    n_taxa_50pct        = numeric(length(core_analysis$detection_thresholds)),
    n_taxa_90pct        = numeric(length(core_analysis$detection_thresholds)),
    n_taxa_100pct       = numeric(length(core_analysis$detection_thresholds))
  )
  
  for(i in seq_along(core_analysis$detection_thresholds)) {
    prevalence_data <- core_analysis$core_results[[i]]
    n_taxa_data$n_taxa_50pct[i]  <- sum(prevalence_data >= 0.5)
    n_taxa_data$n_taxa_90pct[i]  <- sum(prevalence_data >= 0.9)
    n_taxa_data$n_taxa_100pct[i] <- sum(prevalence_data >= 1.0)
  }
  
  n_taxa_long <- reshape2::melt(n_taxa_data, id.vars = "detection_threshold")
  n_taxa_long$prevalence_threshold <- factor(n_taxa_long$variable,
                                             levels = c("n_taxa_50pct", "n_taxa_90pct", "n_taxa_100pct"),
                                             labels = c("≥50% prevalence", "≥90% prevalence", "100% prevalence"))
  
  p_curve <- ggplot(n_taxa_long, aes(x = detection_threshold*100, y = value,
                                     color = prevalence_threshold)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
    scale_x_log10() +
    labs(title = "Core Microbiome Size vs Detection Threshold (tufts excluded)",
         subtitle = paste0(dataset_name, " Combined Environments (Pond + Lake + Desiccated)"),
         x = "Detection Threshold (% Relative Abundance)",
         y = "Number of Taxa",
         color = "Prevalence Threshold") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 14, face = "bold"))
  
  if(save_plot) {
    ggsave(file.path(core_dir, paste0("core_prevalence_curve_", dataset_name, ".pdf")),
           p_curve, width = 10, height = 6, dpi = 300)
    cat("Prevalence curve saved to:", file.path(core_dir, paste0("core_prevalence_curve_", dataset_name, ".pdf")), "\n")
  }
  
  return(p_curve)
}

# Function to analyze core microbiome by environment
analyze_core_by_environment <- function(ps, dataset_name = "16S") {
  
  cat("\n", rep("=", 50), "\n")
  cat("ENVIRONMENT-SPECIFIC CORE MICROBIOME ANALYSIS\n")
  cat(rep("=", 50), "\n")
  
  environments <- c("Lake", "Pond", "Desiccated")
  env_results  <- list()
  
  for(environ in environments) {  # renamed from env to environ
    cat("\nAnalyzing environment:", environ, "\n")
    
    keep_samples <- sample_names(ps)[sample_data(ps)$Environment == environ]
    ps_subset    <- prune_samples(keep_samples, ps)
    
    cat("Samples in", environ, ":", nsamples(ps_subset), "\n")
    cat("Taxa in", environ, ":", ntaxa(ps_subset), "\n")
    
    core_analysis_env <- perform_core_analysis(ps_subset)
    
    heatmap_env <- create_core_heatmap(core_analysis_env,
                                       top_n_taxa = 15,
                                       save_plot = TRUE,
                                       dataset_name = paste0(dataset_name, "_", environ),
                                       title_suffix = paste0(environ, " Environment"))
    
    core_taxa_env <- identify_core_taxa(core_analysis_env,
                                        as.data.frame(tax_table(ps_subset)),
                                        detection_threshold = 0.001,
                                        prevalence_threshold = 0.9)
    
    if(!is.null(core_taxa_env)) {
      write.csv(core_taxa_env,
                file.path(core_dir, paste0("core_taxa_summary_", dataset_name, "_", environ, ".csv")),
                row.names = FALSE)
    }
    
    prevalence_curve_env <- create_prevalence_curve(core_analysis_env,
                                                    save_plot = TRUE,
                                                    dataset_name = paste0(dataset_name, "_", environ))
    
    cat("\n", environ, "Environment Summary:\n")
    thresholds_to_check <- c(0.001, 0.01, 0.1)
    for(thresh in thresholds_to_check) {
      thresh_idx      <- which.min(abs(core_analysis_env$detection_thresholds - thresh))
      prevalence_data <- core_analysis_env$core_results[[thresh_idx]]
      
      cat(paste0("At ", thresh*100, "% detection threshold:\n"))
      cat("- Taxa in ≥50% samples:", sum(prevalence_data >= 0.5), "\n")
      cat("- Taxa in ≥90% samples:", sum(prevalence_data >= 0.9), "\n")
      cat("- Taxa in 100% samples:", sum(prevalence_data >= 1.0), "\n")
    }
    
    env_results[[environ]] <- core_analysis_env
  }
  
  return(env_results)
}

# Function to export full core matrix
export_full_core_matrix <- function(core_analysis, dataset_name = "16S") {
  
  cat("\nExporting full core matrix...\n")
  
  core_matrix <- do.call(cbind, core_analysis$core_results)
  colnames(core_matrix) <- paste0(round(core_analysis$detection_thresholds*100, 4), "%")
  
  ps_rel <- core_analysis$ps_rel
  tax_table_df <- as.data.frame(tax_table(ps_rel))
  
  core_matrix_with_tax <- cbind(tax_table_df[rownames(core_matrix), ], core_matrix)
  
  write.csv(core_matrix_with_tax,
            file.path(core_dir, paste0("full_core_matrix_", dataset_name, ".csv")))
  
  cat("Full core matrix saved to:", file.path(core_dir, paste0("full_core_matrix_", dataset_name, ".csv")), "\n")
  
  return(core_matrix_with_tax)
}

################################################################################
# MAIN ANALYSIS WORKFLOW - 16S rRNA
################################################################################

cat("\n", rep("=", 60), "\n")
cat("CORE MICROBIOME ANALYSIS - 16S COMBINED ENVIRONMENTS (tufts excluded)\n")
cat(rep("=", 60), "\n")

# Step 1: Load notufts RDS and remove tufts
input_filt_rare_16s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/bac_input_filt_rare_notufts.rds")

# Safety filter — removes tufts if still present
input_filt_rare_16s <- filter_data(input_filt_rare_16s, 'Type', filter_vals = 'filament')
cat("16S samples after tuft removal:", ncol(input_filt_rare_16s$data_loaded), "\n")
print(table(input_filt_rare_16s$map_loaded$Environment))

# Check metadata component
cat("Checking data structure...\n")
cat("Available components:", names(input_filt_rare_16s), "\n")

if("map_loaded" %in% names(input_filt_rare_16s)) {
  cat("Using map_loaded for metadata...\n")
  metadata_to_use <- input_filt_rare_16s$map_loaded
} else {
  metadata_to_use <- input_filt_rare_16s$metadata_loaded
}

cat("Metadata dimensions:", dim(metadata_to_use), "\n")
cat("Metadata column names:", colnames(metadata_to_use), "\n")

core_data_16s <- prepare_core_data(input_filt_rare_16s$data_loaded,
                                   metadata_to_use,
                                   input_filt_rare_16s$taxonomy_loaded)

# Step 2: Create phyloseq object
ps_16s <- create_phyloseq_object(core_data_16s)

# Step 3: Perform core analysis (combined)
core_analysis_16s <- perform_core_analysis(ps_16s)

# Step 4: Create heatmap
heatmap_plot_16s <- create_core_heatmap(core_analysis_16s, top_n_taxa = 15, dataset_name = "16S")

# Step 5: Identify core taxa
core_taxa_summary_16s <- identify_core_taxa(core_analysis_16s, core_data_16s$taxonomy,
                                            detection_threshold = 0.001,
                                            prevalence_threshold = 0.9)

if(!is.null(core_taxa_summary_16s)) {
  write.csv(core_taxa_summary_16s,
            file.path(core_dir, "core_taxa_summary_16S.csv"),
            row.names = FALSE)
}

# Step 6: Create prevalence curve
prevalence_curve_16s <- create_prevalence_curve(core_analysis_16s, dataset_name = "16S")

# Step 7: Export full core matrix
full_core_matrix_16s <- export_full_core_matrix(core_analysis_16s, dataset_name = "16S")

# Step 8: Analyze by environment
env_results_16s <- analyze_core_by_environment(ps_16s, dataset_name = "16S")

# Step 9: Summary statistics
cat("\n", rep("=", 50), "\n")
cat("COMBINED ENVIRONMENTS - CORE MICROBIOME SUMMARY (16S)\n")
cat(rep("=", 50), "\n")
cat("Dataset overview:\n")
cat("- Total ASVs analyzed:", ntaxa(ps_16s), "\n")
cat("- Total samples:", nsamples(ps_16s), "\n")
cat("- Environments: Pond, Lake, Desiccated\n")

thresholds_to_check <- c(0.0001, 0.001, 0.01, 0.1)
for(thresh in thresholds_to_check) {
  thresh_idx      <- which.min(abs(core_analysis_16s$detection_thresholds - thresh))
  prevalence_data <- core_analysis_16s$core_results[[thresh_idx]]
  
  cat(paste0("\nAt ", thresh*100, "% detection threshold:\n"))
  cat("- Taxa in ≥50% samples:", sum(prevalence_data >= 0.5), "\n")
  cat("- Taxa in ≥90% samples:", sum(prevalence_data >= 0.9), "\n")
  cat("- Taxa in 100% samples:", sum(prevalence_data >= 1.0), "\n")
}

################################################################################
# MAIN ANALYSIS WORKFLOW - 18S rRNA
################################################################################

cat("\n", rep("=", 60), "\n")
cat("CORE MICROBIOME ANALYSIS - 18S COMBINED ENVIRONMENTS (tufts excluded)\n")
cat(rep("=", 60), "\n")

# Step 1: Load notufts RDS and remove tufts
input_filt_rare_18s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/euk_input_filt_rare_notufts.rds")

# Safety filter — removes tufts if still present
input_filt_rare_18s <- filter_data(input_filt_rare_18s, 'Type', filter_vals = 'filament')
cat("18S samples after tuft removal:", ncol(input_filt_rare_18s$data_loaded), "\n")
print(table(input_filt_rare_18s$map_loaded$Environment))

# Check metadata component
cat("Checking data structure...\n")
cat("Available components:", names(input_filt_rare_18s), "\n")

if("map_loaded" %in% names(input_filt_rare_18s)) {
  cat("Using map_loaded for metadata...\n")
  metadata_to_use <- input_filt_rare_18s$map_loaded
} else {
  metadata_to_use <- input_filt_rare_18s$metadata_loaded
}

cat("Metadata dimensions:", dim(metadata_to_use), "\n")
cat("Metadata column names:", colnames(metadata_to_use), "\n")

core_data_18s <- prepare_core_data(input_filt_rare_18s$data_loaded,
                                   metadata_to_use,
                                   input_filt_rare_18s$taxonomy_loaded)

# Step 2: Create phyloseq object
ps_18s <- create_phyloseq_object(core_data_18s)

# Step 3: Perform core analysis (combined)
core_analysis_18s <- perform_core_analysis(ps_18s)

# Step 4: Create heatmap
heatmap_plot_18s <- create_core_heatmap(core_analysis_18s, top_n_taxa = 15, dataset_name = "18S")

# Step 5: Identify core taxa
core_taxa_summary_18s <- identify_core_taxa(core_analysis_18s, core_data_18s$taxonomy,
                                            detection_threshold = 0.001,
                                            prevalence_threshold = 0.9)

if(!is.null(core_taxa_summary_18s)) {
  write.csv(core_taxa_summary_18s,
            file.path(core_dir, "core_taxa_summary_18S.csv"),
            row.names = FALSE)
}

# Step 6: Create prevalence curve
prevalence_curve_18s <- create_prevalence_curve(core_analysis_18s, dataset_name = "18S")

# Step 7: Export full core matrix
full_core_matrix_18s <- export_full_core_matrix(core_analysis_18s, dataset_name = "18S")

# Step 8: Analyze by environment
env_results_18s <- analyze_core_by_environment(ps_18s, dataset_name = "18S")

# Step 9: Summary statistics
cat("\n", rep("=", 50), "\n")
cat("COMBINED ENVIRONMENTS - CORE MICROBIOME SUMMARY (18S)\n")
cat(rep("=", 50), "\n")
cat("Dataset overview:\n")
cat("- Total ASVs analyzed:", ntaxa(ps_18s), "\n")
cat("- Total samples:", nsamples(ps_18s), "\n")
cat("- Environments: Pond, Lake, Desiccated\n")

thresholds_to_check <- c(0.0001, 0.001, 0.01, 0.1)
for(thresh in thresholds_to_check) {
  thresh_idx      <- which.min(abs(core_analysis_18s$detection_thresholds - thresh))
  prevalence_data <- core_analysis_18s$core_results[[thresh_idx]]
  
  cat(paste0("\nAt ", thresh*100, "% detection threshold:\n"))
  cat("- Taxa in ≥50% samples:", sum(prevalence_data >= 0.5), "\n")
  cat("- Taxa in ≥90% samples:", sum(prevalence_data >= 0.9), "\n")
  cat("- Taxa in 100% samples:", sum(prevalence_data >= 1.0), "\n")
}

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n", rep("=", 60), "\n")
cat("CORE MICROBIOME ANALYSIS COMPLETED (tufts excluded)!\n")
cat(rep("=", 60), "\n")
cat("Output files saved in:", core_dir, "\n")
cat("\nFiles created:\n")
cat("\nCOMBINED ENVIRONMENTS:\n")
cat("- core_microbiome_heatmap_16S.pdf/.png\n")
cat("- core_microbiome_heatmap_18S.pdf/.png\n")
cat("- core_prevalence_curve_16S.pdf\n")
cat("- core_prevalence_curve_18S.pdf\n")
cat("- core_taxa_summary_16S.csv\n")
cat("- core_taxa_summary_18S.csv\n")
cat("- full_core_matrix_16S.csv\n")
cat("- full_core_matrix_18S.csv\n")
cat("\nENVIRONMENT-SPECIFIC:\n")
cat("- core_microbiome_heatmap_16S_Lake/Pond/Desiccated.pdf/.png\n")
cat("- core_microbiome_heatmap_18S_Lake/Pond/Desiccated.pdf/.png\n")
cat("- core_prevalence_curve_16S_Lake/Pond/Desiccated.pdf\n")
cat("- core_prevalence_curve_18S_Lake/Pond/Desiccated.pdf\n")
cat("- core_taxa_summary_16S_Lake/Pond/Desiccated.csv\n")
cat("- core_taxa_summary_18S_Lake/Pond/Desiccated.csv\n")

cat("\nAnalysis complete!\n")

# =============================================================================
# COMBINED FIGURE 4: Core microbiome heatmaps — Bacteria (left) + Eukarya (right)
# Reuses create_core_heatmap()'s exact aesthetics (diverging red/white/blue
# fill scale, theme_minimal(), prevalence-based taxon ordering). Two-step
# process: (1) generate each domain's heatmap panel independently, exported
# without a title/subtitle so they align as clean side-by-side panels,
# (2) composite with precise positioning, following the same reliable
# approach used for Figure 2. Outputs both PDF and high-res JPEG.
# =============================================================================

library(magick)

# -----------------------------------------------------------------------
# Modified heatmap function — same core logic as create_core_heatmap(),
# but returns a clean panel (no title, smaller legend footprint, panel
# label instead) suitable for side-by-side combination
# -----------------------------------------------------------------------
create_core_heatmap_panel <- function(core_analysis, top_n_taxa = 15,
                                      panel_label = "A) Bacteria (16S rRNA gene)") {
  
  core_matrix <- do.call(cbind, core_analysis$core_results)
  colnames(core_matrix) <- paste0(round(core_analysis$detection_thresholds * 100, 4), "%")
  
  first_col_prevalence <- core_matrix[, 1]
  top_taxa <- names(sort(first_col_prevalence, decreasing = TRUE))[1:min(top_n_taxa, nrow(core_matrix))]
  
  core_matrix_subset <- core_matrix[top_taxa, ]
  
  ps_rel <- core_analysis$ps_rel
  tax_table_df <- as.data.frame(tax_table(ps_rel))
  
  genus_labels <- character(length(top_taxa))
  for (i in seq_along(top_taxa)) {
    asv_id <- top_taxa[i]
    if (asv_id %in% rownames(tax_table_df)) {
      best_label <- asv_id
      for (tax_level in 6:1) {
        tax_col <- paste0("taxonomy", tax_level)
        if (tax_col %in% colnames(tax_table_df)) {
          tax_name <- tax_table_df[asv_id, tax_col]
          if (!is.na(tax_name) && tax_name != "" && tax_name != "NA") {
            tax_name_clean <- gsub("_X$", "", tax_name)
            if (tax_name_clean != "") {
              best_label <- paste0(tax_name_clean, " (", asv_id, ")")
              break
            }
          }
        }
      }
      genus_labels[i] <- best_label
    } else {
      genus_labels[i] <- asv_id
    }
  }
  
  rownames(core_matrix_subset) <- genus_labels
  
  core_df <- reshape2::melt(core_matrix_subset)
  colnames(core_df) <- c("Taxon", "Detection_Threshold", "Prevalence")
  
  taxon_avg_prevalence <- aggregate(Prevalence ~ Taxon, core_df, mean)
  taxon_order <- taxon_avg_prevalence$Taxon[order(taxon_avg_prevalence$Prevalence, decreasing = TRUE)]
  
  core_df$Taxon <- factor(core_df$Taxon, levels = rev(taxon_order))
  
  p <- ggplot(core_df, aes(x = Detection_Threshold, y = Taxon, fill = Prevalence)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0.5, name = "Prevalence",
                         labels = function(x) paste0(round(x * 100), "%"),
                         guide = guide_colorbar(title.position = "top")) +
    labs(title = panel_label,
         x = "Detection Threshold (Relative Abundance %)",
         y = "Taxon") +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
          axis.text.y  = element_text(size = 9, face = "bold"),
          axis.title   = element_text(face = "bold", size = 11),
          plot.title   = element_text(size = 13, face = "bold", hjust = 0),
          legend.position = "none",   # legend extracted separately, as in Figure 2
          panel.grid   = element_blank()) +
    scale_x_discrete(breaks = colnames(core_matrix_subset)[seq(1, ncol(core_matrix_subset), 3)])
  
  return(p)
}

# -----------------------------------------------------------------------
# Legend-only version (for separate extraction, matching Figure 2 approach)
# -----------------------------------------------------------------------
create_core_heatmap_legend_source <- function(core_analysis, top_n_taxa = 15) {
  core_matrix <- do.call(cbind, core_analysis$core_results)
  colnames(core_matrix) <- paste0(round(core_analysis$detection_thresholds * 100, 4), "%")
  first_col_prevalence <- core_matrix[, 1]
  top_taxa <- names(sort(first_col_prevalence, decreasing = TRUE))[1:min(top_n_taxa, nrow(core_matrix))]
  core_matrix_subset <- core_matrix[top_taxa, ]
  core_df <- reshape2::melt(core_matrix_subset)
  colnames(core_df) <- c("Taxon", "Detection_Threshold", "Prevalence")
  
  ggplot(core_df, aes(x = Detection_Threshold, y = Taxon, fill = Prevalence)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0.5, name = "Prevalence",
                         labels = function(x) paste0(round(x * 100), "%"),
                         guide = guide_colorbar(title.position = "top")) +
    theme_minimal() +
    theme(legend.title = element_text(face = "bold", size = 11),
          legend.text  = element_text(size = 9, face = "bold"),
          legend.key.height = unit(0.6, "cm"))
}

# -----------------------------------------------------------------------
# STEP 1: Build both domain panels (assumes core_analysis_16s and
# core_analysis_18s already exist from running the main script above)
# -----------------------------------------------------------------------
panel_bac <- create_core_heatmap_panel(core_analysis_16s, top_n_taxa = 15,
                                       panel_label = "A) Bacteria (16S rRNA gene)")
panel_euk <- create_core_heatmap_panel(core_analysis_18s, top_n_taxa = 15,
                                       panel_label = "B) Eukarya (18S rRNA gene)")

legend_source_bac <- create_core_heatmap_legend_source(core_analysis_16s, top_n_taxa = 15)
legend_bac_grob   <- cowplot::get_legend(legend_source_bac)
legend_bac_plot   <- ggdraw() + draw_grob(legend_bac_grob)

# -----------------------------------------------------------------------
# STEP 2: Export panels (side by side) and legend — sized to match your
# reference proportions (wide, short — ~1.5:1 per panel)
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

panels_combined_core <- panel_bac | panel_euk

ggsave("Figure4_panels_only.png", panels_combined_core,
       width = 20, height = 8, dpi = 300, bg = "white")

ggsave("Figure4_legend.png", legend_bac_plot,
       width = 3, height = 5, dpi = 300, bg = "white")

# -----------------------------------------------------------------------
# STEP 3: Composite — legend placed to the right, vertically centered
# on the full combined panel height
# -----------------------------------------------------------------------
dpi <- 300
cm_to_px <- function(cm) round(cm * dpi / 2.54)

panels_img  <- image_read("Figure4_panels_only.png")
legend_core <- image_trim(image_read("Figure4_legend.png"))

panels_info      <- image_info(panels_img)
legend_core_info <- image_info(legend_core)

gap_cm   <- 1
x_offset <- panels_info$width + cm_to_px(gap_cm)
y_offset <- (panels_info$height - legend_core_info$height) / 2

canvas_width  <- panels_info$width + legend_core_info$width + cm_to_px(gap_cm) + cm_to_px(1)
canvas_height <- panels_info$height

final_figure4 <- image_blank(width = canvas_width, height = canvas_height, color = "white")
final_figure4 <- image_composite(final_figure4, panels_img, offset = "+0+0")
final_figure4 <- image_composite(final_figure4, legend_core,
                                 offset = paste0("+", round(x_offset), "+", round(y_offset)))

# -----------------------------------------------------------------------
# STEP 4: Save — high-res JPEG and PDF, both reflecting the updated size
# -----------------------------------------------------------------------
image_write(final_figure4, path = "Figure4_Combined_Final.jpeg", format = "jpeg", quality = 100)
image_write(final_figure4, path = "Figure4_Combined_Final.pdf", format = "pdf", density = dpi)

cat("Saved: Figure4_Combined_Final.jpeg and Figure4_Combined_Final.pdf\n")
cat("Final canvas size:", round(canvas_width/dpi, 2), "x", round(canvas_height/dpi, 2), "inches\n")