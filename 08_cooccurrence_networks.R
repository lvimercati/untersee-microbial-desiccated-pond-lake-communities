# =========================================================================
# 16S NETWORK ANALYSIS — CLR-transformed, significance-filtered, tuft-excluded
# =========================================================================

library(dplyr)
library(Hmisc)        # for rcorr
library(igraph)
library(RColorBrewer)
library(compositions)  # for clr
library(mctoolsr)

set.seed(500)

# -------------------------------------------------------------------------
# 1. LOAD DATA — rarefied, tuft-excluded dataset (consistent with rest of manuscript)
# -------------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/")

input_filt_rare <- readRDS("bac_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')  # excludes tuft samples
cat("Samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

output_dir <- "~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/Network_final_16S"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# -------------------------------------------------------------------------
# 2. GLOBAL COLOR MAP (built once, applied consistently across all habitat networks)
# -------------------------------------------------------------------------
input_filt_rare$taxonomy_loaded$taxonomy2 <- as.character(input_filt_rare$taxonomy_loaded$taxonomy2)
input_filt_rare$taxonomy_loaded$taxonomy2[is.na(input_filt_rare$taxonomy_loaded$taxonomy2) | input_filt_rare$taxonomy_loaded$taxonomy2 == "NA"] <- "Unclassified"
input_filt_rare$taxonomy_loaded$taxonomy2 <- as.factor(input_filt_rare$taxonomy_loaded$taxonomy2)

all_phyla <- levels(input_filt_rare$taxonomy_loaded$taxonomy2)
n_phyla   <- length(all_phyla)
palette_pool <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
if (n_phyla > length(palette_pool)) {
  warning("Not enough unique colors for all phyla; colors will repeat.")
}
GLOBAL_COLORS <- setNames(rep(palette_pool, ceiling(n_phyla / length(palette_pool)))[1:n_phyla], all_phyla)

# -------------------------------------------------------------------------
# 3. CLR TRANSFORMATION
# -------------------------------------------------------------------------
clr_transform_data <- function(count_matrix, pseudocount = 0.5) {
  count_matrix_pc <- count_matrix + pseudocount
  clr_data <- clr(t(count_matrix_pc))
  return(t(clr_data))
}

# -------------------------------------------------------------------------
# 4. NETWORK CONSTRUCTION — CLR correlations, |r| cutoff AND FDR-adjusted p-value cutoff
# -------------------------------------------------------------------------
build_network <- function(data_matrix, taxonomy, cor_cutoff = 0.8, alpha = 0.05, min_samples = 5) {
  
  if (ncol(data_matrix) < min_samples) {
    warning(paste("Only", ncol(data_matrix), "samples available — network not built (min =", min_samples, ")"))
    return(NULL)
  }
  
  cor_matrix <- rcorr(t(data_matrix), type = "spearman")
  r_mat <- cor_matrix$r
  p_mat <- cor_matrix$P
  diag(r_mat) <- 0
  diag(p_mat) <- 1
  r_mat[is.na(r_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  
  # FDR (Benjamini-Hochberg) correction across the unique upper-triangle of tested pairs
  upper_idx <- upper.tri(p_mat)
  p_adj_vec <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat),
                      dimnames = dimnames(p_mat))
  p_adj_mat[upper_idx] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]  # FIXED: proper mirror, not addition
  diag(p_adj_mat) <- 1
  
  # Keep edges that pass BOTH the correlation-strength cutoff AND FDR-adjusted significance
  keep_mask <- (abs(r_mat) >= cor_cutoff) & (p_adj_mat < alpha)
  r_mat_filtered <- r_mat
  r_mat_filtered[!keep_mask] <- 0
  
  net <- graph_from_adjacency_matrix(r_mat_filtered, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[weight == 0])
  net <- delete_vertices(net, which(degree(net) == 0))
  
  V(net)$phylum <- taxonomy[V(net)$name, "taxonomy2"]
  V(net)$color  <- GLOBAL_COLORS[V(net)$phylum]
  
  stats <- list(
    n_vertices      = vcount(net),
    n_edges         = ecount(net),
    n_pos_edges     = sum(E(net)$weight > 0),
    n_neg_edges     = sum(E(net)$weight < 0),
    density         = edge_density(net),
    avg_degree      = mean(degree(net)),
    transitivity    = transitivity(net),
    cor_cutoff      = cor_cutoff,
    fdr_alpha       = alpha,
    n_samples_used  = ncol(data_matrix)
  )
  
  return(list(network = net, r_mat = r_mat, p_adj_mat = p_adj_mat, stats = stats))
}

# -------------------------------------------------------------------------
# 5. PLOTTING
# -------------------------------------------------------------------------
plot_network <- function(net_result, environment, output_dir) {
  
  if (is.null(net_result) || vcount(net_result$network) == 0) {
    cat("No network to plot for", environment, "(insufficient data or no significant edges)\n")
    return(invisible(NULL))
  }
  
  net <- net_result$network
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  deg <- degree(net, mode = "all")
  abs_w <- abs(E(net)$weight)
  net_abs <- net; E(net_abs)$weight <- abs_w
  layout_cl <- layout_with_fr(net_abs)
  
  pdf(file.path(plot_dir, paste0(environment, "_16S_network_CLR_FDR.pdf")), width = 9, height = 6.5)
  par(mar = c(1, 1, 3, 1))
  plot(net,
       vertex.color = V(net)$color,
       vertex.size = pmax(deg * 0.6, 3),
       vertex.shape = "circle",
       vertex.frame.color = "black",
       vertex.label = NA,
       edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved = 0.2,
       edge.width = abs_w * 1.5,
       layout = layout_cl,
       main = paste(environment, "- 16S co-occurrence network\n(CLR, |r| >=",
                    net_result$stats$cor_cutoff, ", FDR-adj p <", net_result$stats$fdr_alpha, ")"))
  
  phyla_in_net <- sort(unique(V(net)$phylum))
  
  if (length(phyla_in_net) > 0) {
    legend("topright", legend = phyla_in_net, pch = 21, pt.bg = GLOBAL_COLORS[phyla_in_net],
           col = "black", pt.cex = 1.3, cex = 0.7, bty = "n", title = "Phylum")
  }
  
  legend("bottomright", legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 2, cex = 0.7, bty = "n", title = "Correlation")
  
  dev.off()
  cat("Plot saved:", file.path(plot_dir, paste0(environment, "_16S_network_CLR_FDR.pdf")), "\n")
}

# -------------------------------------------------------------------------
# 6. MAIN PIPELINE — build + filter + summarize for each habitat
# -------------------------------------------------------------------------
run_env_network <- function(env_name, input_filt_rare, top_n = 200, cor_cutoff = 0.8, alpha = 0.05) {
  
  cat("\n=== PROCESSING", toupper(env_name), "===\n")
  
  env_data <- filter_data(input_filt_rare, 'Environment', keep_vals = env_name)
  counts <- env_data$data_loaded
  cat("Samples in", env_name, ":", ncol(counts), "\n")
  
  # Prevalence filtering — retain ASVs present at reasonable abundance in >=20% of samples
  rel <- sweep(counts, 2, colSums(counts), "/")
  min_samples_prev <- ceiling(0.2 * ncol(counts))
  keep <- rowSums(rel > 0.0001) >= min_samples_prev
  counts_filt <- counts[keep, , drop = FALSE]
  cat("ASVs after prevalence filtering:", nrow(counts_filt), "\n")
  
  # Top-N by total abundance (consistent cutoff for all habitats)
  if (nrow(counts_filt) > top_n) {
    top_asvs <- names(sort(rowSums(counts_filt), decreasing = TRUE))[1:top_n]
    counts_final <- counts_filt[top_asvs, , drop = FALSE]
  } else {
    counts_final <- counts_filt
  }
  cat("ASVs after top-N filtering:", nrow(counts_final), "\n")
  
  # CLR transform
  clr_data <- clr_transform_data(counts_final)
  
  # Build significance-filtered network
  net_result <- build_network(clr_data, input_filt_rare$taxonomy_loaded,
                              cor_cutoff = cor_cutoff, alpha = alpha)
  
  if (!is.null(net_result)) {
    cat("\nNetwork stats for", env_name, ":\n")
    print(net_result$stats)
    plot_network(net_result, env_name, output_dir)
  }
  
  return(net_result)
}

environments <- c("Desiccated", "Pond", "Lake")  # kept in manuscript-standard order
all_results <- list()

for (env in environments) {
  all_results[[env]] <- run_env_network(env, input_filt_rare, top_n = 200, cor_cutoff = 0.8, alpha = 0.05)
}

# -------------------------------------------------------------------------
# 7. SUMMARY TABLE ACROSS HABITATS
# -------------------------------------------------------------------------
summary_data <- do.call(rbind, lapply(names(all_results), function(env) {
  s <- all_results[[env]]$stats
  if (is.null(s)) return(NULL)
  data.frame(Environment = env, N_Vertices = s$n_vertices, N_Edges = s$n_edges,
             Positive = s$n_pos_edges, Negative = s$n_neg_edges,
             Density = round(s$density, 4), Avg_Degree = round(s$avg_degree, 2),
             Transitivity = round(s$transitivity, 3), N_Samples = s$n_samples_used)
}))

print(summary_data)
write.csv(summary_data, file.path(output_dir, "network_summary_16S_CLR_FDR.csv"), row.names = FALSE)

cat("\nDone. Networks now use CLR-transformed data with BOTH |r| >= 0.8 AND BH-adjusted p < 0.05.\n")

# =========================================================================
# 18S NETWORK ANALYSIS — CLR-transformed, significance-filtered, tuft-excluded
# =========================================================================

library(dplyr)
library(Hmisc)        # for rcorr
library(igraph)
library(RColorBrewer)
library(compositions)  # for clr
library(mctoolsr)

set.seed(500)

# -------------------------------------------------------------------------
# 1. LOAD DATA — rarefied, tuft-excluded dataset (consistent with rest of manuscript)
# -------------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/")

input_filt_rare <- readRDS("euk_input_filt_rare_notufts.rds")
input_filt_rare <- filter_data(input_filt_rare, 'Type', filter_vals = 'filament')  # excludes tuft samples
cat("Samples after tuft removal:", ncol(input_filt_rare$data_loaded), "\n")
print(table(input_filt_rare$map_loaded$Environment))

output_dir <- "~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/Network_final_18S"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# -------------------------------------------------------------------------
# 2. GLOBAL COLOR MAP (built once, applied consistently across all habitat networks)
# -------------------------------------------------------------------------
input_filt_rare$taxonomy_loaded$taxonomy2 <- as.character(input_filt_rare$taxonomy_loaded$taxonomy2)
input_filt_rare$taxonomy_loaded$taxonomy2[is.na(input_filt_rare$taxonomy_loaded$taxonomy2) | input_filt_rare$taxonomy_loaded$taxonomy2 == "NA"] <- "Unclassified"
input_filt_rare$taxonomy_loaded$taxonomy2 <- as.factor(input_filt_rare$taxonomy_loaded$taxonomy2)

all_phyla <- levels(input_filt_rare$taxonomy_loaded$taxonomy2)
n_phyla   <- length(all_phyla)
# 18S typically has more phyla than 16S — colorRampPalette to ensure enough distinct colors
palette_pool <- colorRampPalette(brewer.pal(12, "Paired"))(n_phyla)
GLOBAL_COLORS <- setNames(palette_pool, all_phyla)

# -------------------------------------------------------------------------
# 3. CLR TRANSFORMATION
# -------------------------------------------------------------------------
clr_transform_data <- function(count_matrix, pseudocount = 0.5) {
  count_matrix_pc <- count_matrix + pseudocount
  clr_data <- clr(t(count_matrix_pc))
  return(t(clr_data))
}

# -------------------------------------------------------------------------
# 4. NETWORK CONSTRUCTION — CLR correlations, |r| cutoff AND FDR-adjusted p-value cutoff
# -------------------------------------------------------------------------
build_network <- function(data_matrix, taxonomy, cor_cutoff = 0.8, alpha = 0.05, min_samples = 5) {
  
  if (ncol(data_matrix) < min_samples) {
    warning(paste("Only", ncol(data_matrix), "samples available — network not built (min =", min_samples, ")"))
    return(NULL)
  }
  
  cor_matrix <- rcorr(t(data_matrix), type = "spearman")
  r_mat <- cor_matrix$r
  p_mat <- cor_matrix$P
  diag(r_mat) <- 0
  diag(p_mat) <- 1
  r_mat[is.na(r_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  
  # FDR (Benjamini-Hochberg) correction across the unique upper-triangle of tested pairs
  upper_idx <- upper.tri(p_mat)
  p_adj_vec <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat),
                      dimnames = dimnames(p_mat))
  p_adj_mat[upper_idx] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]  # correct mirror
  diag(p_adj_mat) <- 1
  
  # Keep edges that pass BOTH the correlation-strength cutoff AND FDR-adjusted significance
  keep_mask <- (abs(r_mat) >= cor_cutoff) & (p_adj_mat < alpha)
  r_mat_filtered <- r_mat
  r_mat_filtered[!keep_mask] <- 0
  
  net <- graph_from_adjacency_matrix(r_mat_filtered, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[weight == 0])
  net <- delete_vertices(net, which(degree(net) == 0))
  
  V(net)$phylum <- taxonomy[V(net)$name, "taxonomy2"]
  V(net)$color  <- GLOBAL_COLORS[V(net)$phylum]
  
  stats <- list(
    n_vertices      = vcount(net),
    n_edges         = ecount(net),
    n_pos_edges     = sum(E(net)$weight > 0),
    n_neg_edges     = sum(E(net)$weight < 0),
    density         = edge_density(net),
    avg_degree      = mean(degree(net)),
    transitivity    = transitivity(net),
    cor_cutoff      = cor_cutoff,
    fdr_alpha       = alpha,
    n_samples_used  = ncol(data_matrix)
  )
  
  return(list(network = net, r_mat = r_mat, p_adj_mat = p_adj_mat, stats = stats))
}

# -------------------------------------------------------------------------
# 5. PLOTTING
# -------------------------------------------------------------------------
plot_network <- function(net_result, environment, output_dir) {
  
  if (is.null(net_result) || vcount(net_result$network) == 0) {
    cat("No network to plot for", environment, "(insufficient data or no significant edges)\n")
    return(invisible(NULL))
  }
  
  net <- net_result$network
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  deg <- degree(net, mode = "all")
  abs_w <- abs(E(net)$weight)
  net_abs <- net; E(net_abs)$weight <- abs_w
  layout_cl <- layout_with_fr(net_abs)
  
  pdf(file.path(plot_dir, paste0(environment, "_18S_network_CLR_FDR.pdf")), width = 9, height = 6.5)
  par(mar = c(1, 1, 3, 1))
  plot(net,
       vertex.color = V(net)$color,
       vertex.size = pmax(deg * 0.6, 3),
       vertex.shape = "circle",
       vertex.frame.color = "black",
       vertex.label = NA,
       edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved = 0.2,
       edge.width = abs_w * 1.5,
       layout = layout_cl,
       main = paste(environment, "- 18S co-occurrence network\n(CLR, |r| >=",
                    net_result$stats$cor_cutoff, ", FDR-adj p <", net_result$stats$fdr_alpha, ")"))
  
  phyla_in_net <- sort(unique(V(net)$phylum))
  
  if (length(phyla_in_net) > 0) {
    legend("topright", legend = phyla_in_net, pch = 21, pt.bg = GLOBAL_COLORS[phyla_in_net],
           col = "black", pt.cex = 1.3, cex = 0.7, bty = "n", title = "Phylum")
  }
  
  legend("bottomright", legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 2, cex = 0.7, bty = "n", title = "Correlation")
  
  dev.off()
  cat("Plot saved:", file.path(plot_dir, paste0(environment, "_18S_network_CLR_FDR.pdf")), "\n")
}

# -------------------------------------------------------------------------
# 6. MAIN PIPELINE — build + filter + summarize for each habitat
# -------------------------------------------------------------------------
run_env_network <- function(env_name, input_filt_rare, top_n = 200, cor_cutoff = 0.8, alpha = 0.05) {
  
  cat("\n=== PROCESSING", toupper(env_name), "===\n")
  
  env_data <- filter_data(input_filt_rare, 'Environment', keep_vals = env_name)
  counts <- env_data$data_loaded
  cat("Samples in", env_name, ":", ncol(counts), "\n")
  
  # Prevalence filtering — retain ASVs present at reasonable abundance in >=20% of samples
  rel <- sweep(counts, 2, colSums(counts), "/")
  min_samples_prev <- ceiling(0.2 * ncol(counts))
  keep <- rowSums(rel > 0.0001) >= min_samples_prev
  counts_filt <- counts[keep, , drop = FALSE]
  cat("ASVs after prevalence filtering:", nrow(counts_filt), "\n")
  
  # Top-N by total abundance — SAME top_n as 16S for cross-domain comparability
  if (nrow(counts_filt) > top_n) {
    top_asvs <- names(sort(rowSums(counts_filt), decreasing = TRUE))[1:top_n]
    counts_final <- counts_filt[top_asvs, , drop = FALSE]
  } else {
    counts_final <- counts_filt
    cat("NOTE: fewer than", top_n, "ASVs available after prevalence filtering —",
        "using all", nrow(counts_final), "(not top-N capped)\n")
  }
  cat("ASVs after top-N filtering:", nrow(counts_final), "\n")
  
  # CLR transform
  clr_data <- clr_transform_data(counts_final)
  
  # Build significance-filtered network
  net_result <- build_network(clr_data, input_filt_rare$taxonomy_loaded,
                              cor_cutoff = cor_cutoff, alpha = alpha)
  
  if (!is.null(net_result)) {
    cat("\nNetwork stats for", env_name, ":\n")
    print(net_result$stats)
    plot_network(net_result, env_name, output_dir)
  }
  
  return(net_result)
}

environments <- c("Desiccated", "Pond", "Lake")  # kept in manuscript-standard order
all_results <- list()

for (env in environments) {
  all_results[[env]] <- run_env_network(env, input_filt_rare, top_n = 200, cor_cutoff = 0.8, alpha = 0.05)
}

# -------------------------------------------------------------------------
# 7. SUMMARY TABLE ACROSS HABITATS
# -------------------------------------------------------------------------
summary_data <- do.call(rbind, lapply(names(all_results), function(env) {
  s <- all_results[[env]]$stats
  if (is.null(s)) return(NULL)
  data.frame(Environment = env, N_Vertices = s$n_vertices, N_Edges = s$n_edges,
             Positive = s$n_pos_edges, Negative = s$n_neg_edges,
             Density = round(s$density, 4), Avg_Degree = round(s$avg_degree, 2),
             Transitivity = round(s$transitivity, 3), N_Samples = s$n_samples_used)
}))

print(summary_data)
write.csv(summary_data, file.path(output_dir, "network_summary_18S_CLR_FDR.csv"), row.names = FALSE)

cat("\nDone. Networks now use CLR-transformed data with BOTH |r| >= 0.8 AND BH-adjusted p < 0.05.\n")

# =========================================================================
# COMBINED 16S-18S NETWORK ANALYSIS — CLR-transformed, significance-filtered,
# tuft-excluded, cross-domain co-occurrence networks
# =========================================================================

library(dplyr)
library(Hmisc)
library(igraph)
library(RColorBrewer)
library(compositions)
library(mctoolsr)

set.seed(500)

# -------------------------------------------------------------------------
# 1. LOAD BOTH DOMAINS — rarefied, tuft-excluded datasets
# -------------------------------------------------------------------------
input_filt_16s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/bac_input_filt_rare_notufts.rds")
input_filt_16s <- filter_data(input_filt_16s, 'Type', filter_vals = 'filament')
cat("16S samples after tuft removal:", ncol(input_filt_16s$data_loaded), "\n")

input_filt_18s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/euk_input_filt_rare_notufts.rds")
input_filt_18s <- filter_data(input_filt_18s, 'Type', filter_vals = 'filament')
cat("18S samples after tuft removal:", ncol(input_filt_18s$data_loaded), "\n")

output_dir <- "~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/Network_final_Combined"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# -------------------------------------------------------------------------
# 2. GLOBAL COLOR MAPS — built once per domain, reused across all habitat networks
# -------------------------------------------------------------------------
prep_phylum_factor <- function(input_obj) {
  input_obj$taxonomy_loaded$taxonomy2 <- as.character(input_obj$taxonomy_loaded$taxonomy2)
  input_obj$taxonomy_loaded$taxonomy2[is.na(input_obj$taxonomy_loaded$taxonomy2) | input_obj$taxonomy_loaded$taxonomy2 == "NA"] <- "Unclassified"
  input_obj$taxonomy_loaded$taxonomy2 <- as.factor(input_obj$taxonomy_loaded$taxonomy2)
  input_obj
}
input_filt_16s <- prep_phylum_factor(input_filt_16s)
input_filt_18s <- prep_phylum_factor(input_filt_18s)

phyla_16s <- levels(input_filt_16s$taxonomy_loaded$taxonomy2)
GLOBAL_16S_COLORS <- setNames(
  rep(c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")),
      length.out = length(phyla_16s)),
  phyla_16s
)

phyla_18s <- levels(input_filt_18s$taxonomy_loaded$taxonomy2)
GLOBAL_18S_COLORS <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(phyla_18s)),
  phyla_18s
)

# -------------------------------------------------------------------------
# 3. CLR TRANSFORMATION
# -------------------------------------------------------------------------
clr_transform_data <- function(count_matrix, pseudocount = 0.5) {
  count_matrix_pc <- count_matrix + pseudocount
  clr_data <- clr(t(count_matrix_pc))
  return(t(clr_data))
}

# -------------------------------------------------------------------------
# 4. NETWORK CONSTRUCTION — shared FDR-corrected builder (single or combined domain)
# -------------------------------------------------------------------------
build_network_core <- function(data_matrix, cor_cutoff = 0.8, alpha = 0.05, min_samples = 5) {
  
  if (ncol(data_matrix) < min_samples) {
    warning(paste("Only", ncol(data_matrix), "samples available — network not built"))
    return(NULL)
  }
  
  cor_matrix <- rcorr(t(data_matrix), type = "spearman")
  r_mat <- cor_matrix$r
  p_mat <- cor_matrix$P
  diag(r_mat) <- 0
  diag(p_mat) <- 1
  r_mat[is.na(r_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  
  upper_idx <- upper.tri(p_mat)
  p_adj_vec <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
  p_adj_mat[upper_idx] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  diag(p_adj_mat) <- 1
  
  keep_mask <- (abs(r_mat) >= cor_cutoff) & (p_adj_mat < alpha)
  r_mat_filtered <- r_mat
  r_mat_filtered[!keep_mask] <- 0
  
  net <- graph_from_adjacency_matrix(r_mat_filtered, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[weight == 0])
  net <- delete_vertices(net, which(degree(net) == 0))
  
  list(network = net, r_mat = r_mat, p_adj_mat = p_adj_mat,
       cor_cutoff = cor_cutoff, alpha = alpha, n_samples_used = ncol(data_matrix))
}

build_single_network <- function(data_matrix, taxonomy, colors, cor_cutoff = 0.8, alpha = 0.05) {
  res <- build_network_core(data_matrix, cor_cutoff, alpha)
  if (is.null(res) || vcount(res$network) == 0) return(res)
  
  V(res$network)$phylum <- taxonomy[V(res$network)$name, "taxonomy2"]
  V(res$network)$color  <- colors[V(res$network)$phylum]
  
  res$stats <- list(
    n_vertices = vcount(res$network), n_edges = ecount(res$network),
    n_pos_edges = sum(E(res$network)$weight > 0), n_neg_edges = sum(E(res$network)$weight < 0),
    density = edge_density(res$network), avg_degree = mean(degree(res$network)),
    transitivity = transitivity(res$network),
    cor_cutoff = cor_cutoff, fdr_alpha = alpha, n_samples_used = res$n_samples_used
  )
  res
}

build_combined_network <- function(combined_clr, combined_taxonomy, cor_cutoff = 0.8, alpha = 0.05) {
  res <- build_network_core(combined_clr, cor_cutoff, alpha)
  if (is.null(res) || vcount(res$network) == 0) return(res)
  
  net <- res$network
  V(net)$phylum <- combined_taxonomy[V(net)$name, "taxonomy2"]
  V(net)$source <- combined_taxonomy[V(net)$name, "source"]
  V(net)$phylum_source <- paste(V(net)$phylum, V(net)$source, sep = "_")
  
  bac_phyla <- unique(V(net)$phylum[V(net)$source == "16S"])
  euk_phyla <- unique(V(net)$phylum[V(net)$source == "18S"])
  
  bac_color_lookup <- setNames(GLOBAL_16S_COLORS[bac_phyla], paste(bac_phyla, "16S", sep = "_"))
  euk_color_lookup <- setNames(GLOBAL_18S_COLORS[euk_phyla], paste(euk_phyla, "18S", sep = "_"))
  
  all_colors <- c(bac_color_lookup, euk_color_lookup)
  V(net)$color <- all_colors[V(net)$phylum_source]
  
  bac_v <- which(V(net)$source == "16S")
  euk_v <- which(V(net)$source == "18S")
  cross_e <- E(net)[bac_v %--% euk_v]
  within_bac_e <- E(net)[bac_v %--% bac_v]
  within_euk_e <- E(net)[euk_v %--% euk_v]
  
  res$network <- net
  res$stats <- list(
    n_vertices = vcount(net), n_edges = ecount(net),
    n_16s = length(bac_v), n_18s = length(euk_v),
    cross_domain = length(cross_e), within_bac = length(within_bac_e), within_euk = length(within_euk_e),
    cross_pos = sum(E(net)[cross_e]$weight > 0), cross_neg = sum(E(net)[cross_e]$weight < 0),
    density = edge_density(net), cor_cutoff = cor_cutoff, fdr_alpha = alpha,
    n_samples_used = res$n_samples_used
  )
  res
}

# -------------------------------------------------------------------------
# 5. DOUBLE-RING LAYOUT FOR COMBINED NETWORK PLOT
# -------------------------------------------------------------------------
create_double_ring_layout <- function(net, inner_radius = 0.6, outer_radius = 1.0) {
  bac_idx <- which(V(net)$source == "16S")
  euk_idx <- which(V(net)$source == "18S")
  
  layout_matrix <- matrix(0, nrow = vcount(net), ncol = 2)
  
  if (length(bac_idx) > 0) {
    ord <- bac_idx[order(V(net)$phylum[bac_idx])]
    ang <- seq(0, 2 * pi, length.out = length(bac_idx) + 1)[-(length(bac_idx) + 1)]
    layout_matrix[ord, 1] <- inner_radius * cos(ang)
    layout_matrix[ord, 2] <- inner_radius * sin(ang)
  }
  if (length(euk_idx) > 0) {
    ord <- euk_idx[order(V(net)$phylum[euk_idx])]
    ang <- seq(0, 2 * pi, length.out = length(euk_idx) + 1)[-(length(euk_idx) + 1)]
    layout_matrix[ord, 1] <- outer_radius * cos(ang)
    layout_matrix[ord, 2] <- outer_radius * sin(ang)
  }
  layout_matrix
}

# -------------------------------------------------------------------------
# 6. PLOTTING — SINGLE-HABITAT COMBINED NETWORK (double ring, split-panel legend)
#    Legend now uses a clamped y-position so Correlations can NEVER be pushed
#    off-panel, regardless of how many phyla are listed above it.
# -------------------------------------------------------------------------
plot_combined_network <- function(net_result, environment, output_dir) {
  
  if (is.null(net_result) || vcount(net_result$network) == 0) {
    cat("No combined network to plot for", environment, "\n")
    return(invisible(NULL))
  }
  
  net <- net_result$network
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  layout_matrix <- create_double_ring_layout(net)
  
  pdf(file.path(plot_dir, paste0(environment, "_combined_network_CLR_FDR.pdf")), width = 14, height = 10)
  layout(matrix(c(1, 2), ncol = 2), widths = c(2.5, 1))
  
  par(mar = c(1, 1, 3, 0))
  plot(net,
       vertex.color = V(net)$color,
       vertex.size = pmax(degree(net) * 0.5, 3),
       vertex.shape = ifelse(V(net)$source == "16S", "circle", "square"),
       vertex.frame.color = "black",
       vertex.label = NA,
       edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved = 0.3,
       edge.width = abs(E(net)$weight) * 2,
       layout = layout_matrix,
       xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), rescale = FALSE, asp = 1,
       main = paste(environment, "Environment: 16S-18S Network (CLR, |r|>=",
                    net_result$stats$cor_cutoff, ", FDR-adj p<", net_result$stats$fdr_alpha, ")"))
  
  text(0.7, 0, "16S", font = 2, cex = 1.2)
  text(1.15, 0, "18S", font = 2, cex = 1.2)
  
  par(mar = c(1, 0, 3, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Legend", cex.main = 1.2)
  
  bac_phyla_in_net <- sort(unique(V(net)$phylum[V(net)$source == "16S"]))
  euk_phyla_in_net <- sort(unique(V(net)$phylum[V(net)$source == "18S"]))
  
  n_legend_rows <- length(bac_phyla_in_net) + length(euk_phyla_in_net) + 4
  y_top <- 0.95
  line_height <- min(1 / n_legend_rows, 0.04)
  legend_cex <- ifelse(n_legend_rows > 25, 0.7, 0.9)
  
  if (length(bac_phyla_in_net) > 0) {
    legend(0.02, y_top, legend = bac_phyla_in_net, pch = 21, pt.bg = GLOBAL_16S_COLORS[bac_phyla_in_net],
           col = "black", pt.cex = 1.3, cex = legend_cex, title = "16S Taxa", title.font = 2, bty = "n")
    y_top <- max(y_top - (length(bac_phyla_in_net) + 2) * line_height, 0.15)
  }
  if (length(euk_phyla_in_net) > 0) {
    legend(0.02, y_top, legend = euk_phyla_in_net, pch = 22, pt.bg = GLOBAL_18S_COLORS[euk_phyla_in_net],
           col = "black", pt.cex = 1.3, cex = legend_cex, title = "18S Taxa", title.font = 2, bty = "n")
    y_top <- max(y_top - (length(euk_phyla_in_net) + 2) * line_height, 0.08)
  }
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 2, cex = legend_cex, title = "Correlations",
         title.font = 2, bty = "n")
  
  dev.off()
  cat("Plot saved:", file.path(plot_dir, paste0(environment, "_combined_network_CLR_FDR.pdf")), "\n")
}

# -------------------------------------------------------------------------
# 7. MAIN PIPELINE — per habitat: prevalence filter, top-N (domain-specific for
#    single networks, EQUAL for the combined plot), CLR, build, plot
# -------------------------------------------------------------------------
process_environment <- function(env_name, input_16s, input_18s,
                                top_n_16s = 200, top_n_18s = 100,
                                top_n_combined = 75,   # SAME value applied to BOTH domains for the combined plot
                                cor_cutoff = 0.8, alpha = 0.05) {
  
  cat("\n=== PROCESSING", toupper(env_name), "===\n")
  
  env_16s <- filter_data(input_16s, 'Environment', keep_vals = env_name)
  env_18s <- filter_data(input_18s, 'Environment', keep_vals = env_name)
  
  common_samples <- intersect(colnames(env_16s$data_loaded), colnames(env_18s$data_loaded))
  cat("Common samples:", length(common_samples), "\n")
  if (length(common_samples) < 5) {
    cat("WARNING: insufficient common samples — skipping\n")
    return(NULL)
  }
  
  bac_counts <- env_16s$data_loaded[, common_samples, drop = FALSE]
  euk_counts <- env_18s$data_loaded[, common_samples, drop = FALSE]
  
  filter_and_topn <- function(counts, top_n) {
    rel <- sweep(counts, 2, colSums(counts), "/")
    min_s <- ceiling(0.2 * ncol(counts))
    keep <- rowSums(rel > 0.0001) >= min_s
    counts_f <- counts[keep, , drop = FALSE]
    if (nrow(counts_f) > top_n) {
      top_ids <- names(sort(rowSums(counts_f), decreasing = TRUE))[1:top_n]
      counts_f <- counts_f[top_ids, , drop = FALSE]
    }
    counts_f
  }
  
  # --- FULL, domain-appropriate pools, for standalone single-domain networks ---
  bac_final <- filter_and_topn(bac_counts, top_n_16s)
  euk_final <- filter_and_topn(euk_counts, top_n_18s)
  cat("16S ASVs (single):", nrow(bac_final), "| 18S ASVs (single):", nrow(euk_final), "\n")
  
  bac_clr <- clr_transform_data(bac_final)
  euk_clr <- clr_transform_data(euk_final)
  
  results <- list()
  
  if (nrow(bac_final) > 5) {
    results$net_16s <- build_single_network(bac_clr, input_16s$taxonomy_loaded, GLOBAL_16S_COLORS,
                                            cor_cutoff, alpha)
  }
  if (nrow(euk_final) > 5) {
    results$net_18s <- build_single_network(euk_clr, input_18s$taxonomy_loaded, GLOBAL_18S_COLORS,
                                            cor_cutoff, alpha)
  }
  
  # --- EQUAL pools (top_n_combined for BOTH domains), specifically for the combined plot ---
  bac_final_comb <- filter_and_topn(bac_counts, top_n_combined)
  euk_final_comb <- filter_and_topn(euk_counts, top_n_combined)
  cat("16S ASVs (combined):", nrow(bac_final_comb), "| 18S ASVs (combined):", nrow(euk_final_comb), "\n")
  
  if (nrow(bac_final_comb) > 5 && nrow(euk_final_comb) > 5) {
    bac_clr_comb <- clr_transform_data(bac_final_comb)
    euk_clr_comb <- clr_transform_data(euk_final_comb)
    combined_clr <- rbind(bac_clr_comb, euk_clr_comb)
    
    bac_tax <- input_16s$taxonomy_loaded[rownames(bac_clr_comb), , drop = FALSE]
    euk_tax <- input_18s$taxonomy_loaded[rownames(euk_clr_comb), , drop = FALSE]
    
    bac_tax$taxonomy2 <- as.character(bac_tax$taxonomy2)
    bac_tax$taxonomy2[is.na(bac_tax$taxonomy2) | bac_tax$taxonomy2 == "NA"] <- "Unclassified"
    euk_tax$taxonomy2 <- as.character(euk_tax$taxonomy2)
    euk_tax$taxonomy2[is.na(euk_tax$taxonomy2) | euk_tax$taxonomy2 == "NA"] <- "Unclassified"
    
    bac_tax$source <- "16S"
    euk_tax$source <- "18S"
    combined_tax <- rbind(bac_tax[, c("taxonomy2", "source")], euk_tax[, c("taxonomy2", "source")])
    
    results$net_combined <- build_combined_network(combined_clr, combined_tax, cor_cutoff, alpha)
    plot_combined_network(results$net_combined, env_name, output_dir)
  }
  
  if (!is.null(results$net_16s)) { cat("16S network:\n"); print(results$net_16s$stats) }
  if (!is.null(results$net_18s)) { cat("18S network:\n"); print(results$net_18s$stats) }
  if (!is.null(results$net_combined)) { cat("Combined network:\n"); print(results$net_combined$stats) }
  
  results
}

# -------------------------------------------------------------------------
# 8. RUN FOR ALL HABITATS (Desiccated -> Pond -> Lake)
#    Combined-plot cap now set to 75/75 (equal, per your request), while
#    single-domain networks stay at their own domain-appropriate caps.
# -------------------------------------------------------------------------
environments <- c("Desiccated", "Pond", "Lake")
all_results <- list()
for (env in environments) {
  all_results[[env]] <- process_environment(env, input_filt_16s, input_filt_18s,
                                            top_n_16s = 200, top_n_18s = 100,
                                            top_n_combined = 75,
                                            cor_cutoff = 0.8, alpha = 0.05)
}

# -------------------------------------------------------------------------
# 9. SUMMARY TABLE — cross-domain comparison across habitats
# -------------------------------------------------------------------------
summary_rows <- list()
for (env in names(all_results)) {
  res <- all_results[[env]]
  if (is.null(res)) next
  if (!is.null(res$net_combined) && !is.null(res$net_combined$stats)) {
    s <- res$net_combined$stats
    cross_pct <- if (s$cross_domain > 0) round(100 * s$cross_pos / s$cross_domain, 1) else NA
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Environment = env, N_16S = s$n_16s, N_18S = s$n_18s,
      Within_16S_Edges = s$within_bac, Within_18S_Edges = s$within_euk,
      Cross_Domain_Edges = s$cross_domain, Cross_Domain_Pct_Positive = cross_pct,
      Density = round(s$density, 4)
    )
  }
}
summary_data <- do.call(rbind, summary_rows)
print(summary_data)
write.csv(summary_data, file.path(output_dir, "network_summary_combined_CLR_FDR.csv"), row.names = FALSE)

cat("\nDone. Combined 16S-18S networks now use CLR-transformed data,\n",
    "BOTH |r| >= 0.8 AND BH-adjusted p < 0.05, correct FDR mirroring, tuft exclusion,\n",
    "and an equal top_n=75 cap per domain for combined-plot readability.\n")

# =========================================================================
# 10. COMBINED ROW FIGURE — Desiccated | Pond | Lake, side by side
#     (single, cleaned-up version — legend now clamped to avoid overlap)
# =========================================================================
plot_combined_row_figure <- function(all_results, output_dir,
                                     env_order = c("Desiccated", "Pond", "Lake"),
                                     inner_radius = 0.6, outer_radius = 1.0) {
  
  valid_envs <- env_order[sapply(env_order, function(e) {
    !is.null(all_results[[e]]) &&
      !is.null(all_results[[e]]$net_combined) &&
      !is.null(all_results[[e]]$net_combined$network) &&
      vcount(all_results[[e]]$net_combined$network) > 0
  })]
  
  if (length(valid_envs) == 0) {
    cat("No valid combined networks available to plot.\n")
    return(invisible(NULL))
  }
  
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  all_bac_phyla <- character(0)
  all_euk_phyla <- character(0)
  for (env in valid_envs) {
    net <- all_results[[env]]$net_combined$network
    all_bac_phyla <- union(all_bac_phyla, unique(V(net)$phylum[V(net)$source == "16S"]))
    all_euk_phyla <- union(all_euk_phyla, unique(V(net)$phylum[V(net)$source == "18S"]))
  }
  all_bac_phyla <- sort(all_bac_phyla)
  all_euk_phyla <- sort(all_euk_phyla)
  
  n_legend_rows <- length(all_bac_phyla) + length(all_euk_phyla) + 4
  fig_height <- max(6.5, n_legend_rows * 0.28)
  
  pdf(file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row.pdf"),
      width = 6 * length(valid_envs) + 3.8, height = fig_height)
  
  layout(matrix(1:(length(valid_envs) + 1), nrow = 1),
         widths = c(rep(5, length(valid_envs)), 2.8))
  
  for (env in valid_envs) {
    net_result <- all_results[[env]]$net_combined
    net <- net_result$network
    layout_matrix <- create_double_ring_layout(net, inner_radius, outer_radius)
    
    par(mar = c(1, 1, 5, 1))
    plot(net,
         vertex.color = V(net)$color,
         vertex.size = pmax(degree(net) * 0.5, 3),
         vertex.shape = ifelse(V(net)$source == "16S", "circle", "square"),
         vertex.frame.color = "black",
         vertex.label = NA,
         edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
         edge.curved = 0.3,
         edge.width = abs(E(net)$weight) * 2,
         layout = layout_matrix,
         main = "",
         xlim = c(-outer_radius * 1.3, outer_radius * 1.3),
         ylim = c(-outer_radius * 1.3, outer_radius * 1.3),
         rescale = FALSE, asp = 1)
    
    title(main = env, cex.main = 2.2, font.main = 2, line = 2.5)
    
    text(inner_radius * 1.15, 0, "16S", font = 2, cex = 1.4, col = "black")
    text(outer_radius * 1.15, 0, "18S", font = 2, cex = 1.4, col = "black")
  }
  
  par(mar = c(1, 0, 5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Legend", cex.main = 2.0, font.main = 2, line = 2.5)
  
  y_top <- 0.95
  line_height <- min(1 / n_legend_rows, 0.035)
  
  if (length(all_bac_phyla) > 0) {
    legend(0.02, y_top, legend = all_bac_phyla, pch = 21, pt.bg = GLOBAL_16S_COLORS[all_bac_phyla],
           col = "black", pt.cex = 1.6, cex = 0.95, title = "16S Taxa",
           title.font = 2, bty = "n")
    y_top <- max(y_top - (length(all_bac_phyla) + 2) * line_height, 0.15)
  }
  
  if (length(all_euk_phyla) > 0) {
    legend(0.02, y_top, legend = all_euk_phyla, pch = 22, pt.bg = GLOBAL_18S_COLORS[all_euk_phyla],
           col = "black", pt.cex = 1.6, cex = 0.95, title = "18S Taxa",
           title.font = 2, bty = "n")
    y_top <- max(y_top - (length(all_euk_phyla) + 2) * line_height, 0.08)
  }
  
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 3, cex = 1.0, title = "Correlations",
         title.font = 2, bty = "n")
  
  dev.off()
  cat("Combined row figure saved:", file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row.pdf"), "\n")
}

plot_combined_row_figure(all_results, output_dir, env_order = c("Desiccated", "Pond", "Lake"))

# =========================================================================
# 11. SINGLE-DOMAIN ROW FIGURES — Desiccated | Pond | Lake, one row per domain
#     (legend now clamped to avoid overlap, same fix as above)
# =========================================================================
plot_single_domain_row_figure <- function(all_results, domain, output_dir,
                                          env_order = c("Desiccated", "Pond", "Lake")) {
  
  net_field <- if (domain == "16S") "net_16s" else "net_18s"
  colors    <- if (domain == "16S") GLOBAL_16S_COLORS else GLOBAL_18S_COLORS
  
  valid_envs <- env_order[sapply(env_order, function(e) {
    !is.null(all_results[[e]]) &&
      !is.null(all_results[[e]][[net_field]]) &&
      !is.null(all_results[[e]][[net_field]]$network) &&
      vcount(all_results[[e]][[net_field]]$network) > 0
  })]
  
  if (length(valid_envs) == 0) {
    cat("No valid", domain, "networks available to plot.\n")
    return(invisible(NULL))
  }
  
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  all_phyla <- character(0)
  for (env in valid_envs) {
    net <- all_results[[env]][[net_field]]$network
    all_phyla <- union(all_phyla, unique(V(net)$phylum))
  }
  all_phyla <- sort(all_phyla)
  
  n_legend_rows <- length(all_phyla) + 3
  fig_height <- max(6.5, n_legend_rows * 0.28)
  
  pdf(file.path(plot_dir, paste0(domain, "_Network_AllEnvironments_Row.pdf")),
      width = 6 * length(valid_envs) + 3.5, height = fig_height)
  
  layout(matrix(1:(length(valid_envs) + 1), nrow = 1),
         widths = c(rep(5, length(valid_envs)), 2.6))
  
  for (env in valid_envs) {
    net_result <- all_results[[env]][[net_field]]
    net <- net_result$network
    
    deg <- degree(net, mode = "all")
    abs_w <- abs(E(net)$weight)
    net_abs <- net; E(net_abs)$weight <- abs_w
    layout_cl <- layout_with_fr(net_abs)
    
    par(mar = c(1, 1, 5, 1))
    plot(net,
         vertex.color = V(net)$color,
         vertex.size = pmax(deg * 0.6, 3),
         vertex.shape = "circle",
         vertex.frame.color = "black",
         vertex.label = NA,
         edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
         edge.curved = 0.2,
         edge.width = abs_w * 1.5,
         layout = layout_cl,
         main = "")
    
    title(main = env, cex.main = 2.2, font.main = 2, line = 2.5)
  }
  
  par(mar = c(1, 0, 5, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Legend", cex.main = 2.0, font.main = 2, line = 2.5)
  
  y_top <- 0.95
  line_height <- min(1 / n_legend_rows, 0.035)
  
  if (length(all_phyla) > 0) {
    legend(0.02, y_top, legend = all_phyla, pch = 21, pt.bg = colors[all_phyla],
           col = "black", pt.cex = 1.6, cex = 0.95, title = paste(domain, "Taxa"),
           title.font = 2, bty = "n")
    y_top <- max(y_top - (length(all_phyla) + 2) * line_height, 0.1)
  }
  
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 3, cex = 1.0, title = "Correlations",
         title.font = 2, bty = "n")
  
  dev.off()
  cat("Row figure saved:", file.path(plot_dir, paste0(domain, "_Network_AllEnvironments_Row.pdf")), "\n")
}

plot_single_domain_row_figure(all_results, "16S", output_dir, env_order = c("Desiccated", "Pond", "Lake"))
plot_single_domain_row_figure(all_results, "18S", output_dir, env_order = c("Desiccated", "Pond", "Lake"))

# =============================================================================
# FIGURE 9 (STANDALONE, FINAL): Combined 16S-18S network row figure
# Desiccated | Pond | Lake, cross-domain double-ring networks.
# This script is self-contained: it loads all required data, defines every
# helper function, rebuilds all_results from scratch, and produces the
# final corrected row figure. No need to run any other script first.
#
# Fixes applied vs. the original figure:
#   - "16S"/"18S" ring labels moved much closer to their respective rings
#   - "Legend" title removed
#   - Legend taxon-name text enlarged
#   - 16S/18S legend blocks placed closer together
# =============================================================================

library(dplyr)
library(Hmisc)        # for rcorr
library(igraph)
library(RColorBrewer)
library(compositions)  # for clr
library(mctoolsr)

set.seed(500)

# -------------------------------------------------------------------------
# 1. LOAD BOTH DOMAINS — rarefied, tuft-excluded datasets
# -------------------------------------------------------------------------
input_filt_16s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_16S/bac_input_filt_rare_notufts.rds")
input_filt_16s <- filter_data(input_filt_16s, 'Type', filter_vals = 'filament')
cat("16S samples after tuft removal:", ncol(input_filt_16s$data_loaded), "\n")

input_filt_18s <- readRDS("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_18S/euk_input_filt_rare_notufts.rds")
input_filt_18s <- filter_data(input_filt_18s, 'Type', filter_vals = 'filament')
cat("18S samples after tuft removal:", ncol(input_filt_18s$data_loaded), "\n")

output_dir <- "~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/Network_final_Combined"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

# -------------------------------------------------------------------------
# 2. GLOBAL COLOR MAPS — built once per domain, reused across all habitat networks
# -------------------------------------------------------------------------
prep_phylum_factor <- function(input_obj) {
  input_obj$taxonomy_loaded$taxonomy2 <- as.character(input_obj$taxonomy_loaded$taxonomy2)
  input_obj$taxonomy_loaded$taxonomy2[is.na(input_obj$taxonomy_loaded$taxonomy2) | input_obj$taxonomy_loaded$taxonomy2 == "NA"] <- "Unclassified"
  input_obj$taxonomy_loaded$taxonomy2 <- as.factor(input_obj$taxonomy_loaded$taxonomy2)
  input_obj
}
input_filt_16s <- prep_phylum_factor(input_filt_16s)
input_filt_18s <- prep_phylum_factor(input_filt_18s)

phyla_16s <- levels(input_filt_16s$taxonomy_loaded$taxonomy2)
GLOBAL_16S_COLORS <- setNames(
  rep(c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")),
      length.out = length(phyla_16s)),
  phyla_16s
)

phyla_18s <- levels(input_filt_18s$taxonomy_loaded$taxonomy2)
GLOBAL_18S_COLORS <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(phyla_18s)),
  phyla_18s
)

# -------------------------------------------------------------------------
# 3. CLR TRANSFORMATION
# -------------------------------------------------------------------------
clr_transform_data <- function(count_matrix, pseudocount = 0.5) {
  count_matrix_pc <- count_matrix + pseudocount
  clr_data <- clr(t(count_matrix_pc))
  return(t(clr_data))
}

# -------------------------------------------------------------------------
# 4. NETWORK CONSTRUCTION — shared FDR-corrected builder (single or combined domain)
# -------------------------------------------------------------------------
build_network_core <- function(data_matrix, cor_cutoff = 0.8, alpha = 0.05, min_samples = 5) {
  
  if (ncol(data_matrix) < min_samples) {
    warning(paste("Only", ncol(data_matrix), "samples available — network not built"))
    return(NULL)
  }
  
  cor_matrix <- rcorr(t(data_matrix), type = "spearman")
  r_mat <- cor_matrix$r
  p_mat <- cor_matrix$P
  diag(r_mat) <- 0
  diag(p_mat) <- 1
  r_mat[is.na(r_mat)] <- 0
  p_mat[is.na(p_mat)] <- 1
  
  upper_idx <- upper.tri(p_mat)
  p_adj_vec <- p.adjust(p_mat[upper_idx], method = "BH")
  p_adj_mat <- matrix(1, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
  p_adj_mat[upper_idx] <- p_adj_vec
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  diag(p_adj_mat) <- 1
  
  keep_mask <- (abs(r_mat) >= cor_cutoff) & (p_adj_mat < alpha)
  r_mat_filtered <- r_mat
  r_mat_filtered[!keep_mask] <- 0
  
  net <- graph_from_adjacency_matrix(r_mat_filtered, mode = "lower", weighted = TRUE)
  net <- delete_edges(net, E(net)[weight == 0])
  net <- delete_vertices(net, which(degree(net) == 0))
  
  list(network = net, r_mat = r_mat, p_adj_mat = p_adj_mat,
       cor_cutoff = cor_cutoff, alpha = alpha, n_samples_used = ncol(data_matrix))
}

build_single_network <- function(data_matrix, taxonomy, colors, cor_cutoff = 0.8, alpha = 0.05) {
  res <- build_network_core(data_matrix, cor_cutoff, alpha)
  if (is.null(res) || vcount(res$network) == 0) return(res)
  
  V(res$network)$phylum <- taxonomy[V(res$network)$name, "taxonomy2"]
  V(res$network)$color  <- colors[V(res$network)$phylum]
  
  res$stats <- list(
    n_vertices = vcount(res$network), n_edges = ecount(res$network),
    n_pos_edges = sum(E(res$network)$weight > 0), n_neg_edges = sum(E(res$network)$weight < 0),
    density = edge_density(res$network), avg_degree = mean(degree(res$network)),
    transitivity = transitivity(res$network),
    cor_cutoff = cor_cutoff, fdr_alpha = alpha, n_samples_used = res$n_samples_used
  )
  res
}

build_combined_network <- function(combined_clr, combined_taxonomy, cor_cutoff = 0.8, alpha = 0.05) {
  res <- build_network_core(combined_clr, cor_cutoff, alpha)
  if (is.null(res) || vcount(res$network) == 0) return(res)
  
  net <- res$network
  V(net)$phylum <- combined_taxonomy[V(net)$name, "taxonomy2"]
  V(net)$source <- combined_taxonomy[V(net)$name, "source"]
  V(net)$phylum_source <- paste(V(net)$phylum, V(net)$source, sep = "_")
  
  bac_phyla <- unique(V(net)$phylum[V(net)$source == "16S"])
  euk_phyla <- unique(V(net)$phylum[V(net)$source == "18S"])
  
  bac_color_lookup <- setNames(GLOBAL_16S_COLORS[bac_phyla], paste(bac_phyla, "16S", sep = "_"))
  euk_color_lookup <- setNames(GLOBAL_18S_COLORS[euk_phyla], paste(euk_phyla, "18S", sep = "_"))
  
  all_colors <- c(bac_color_lookup, euk_color_lookup)
  V(net)$color <- all_colors[V(net)$phylum_source]
  
  bac_v <- which(V(net)$source == "16S")
  euk_v <- which(V(net)$source == "18S")
  cross_e <- E(net)[bac_v %--% euk_v]
  within_bac_e <- E(net)[bac_v %--% bac_v]
  within_euk_e <- E(net)[euk_v %--% euk_v]
  
  res$network <- net
  res$stats <- list(
    n_vertices = vcount(net), n_edges = ecount(net),
    n_16s = length(bac_v), n_18s = length(euk_v),
    cross_domain = length(cross_e), within_bac = length(within_bac_e), within_euk = length(within_euk_e),
    cross_pos = sum(E(net)[cross_e]$weight > 0), cross_neg = sum(E(net)[cross_e]$weight < 0),
    density = edge_density(net), cor_cutoff = cor_cutoff, fdr_alpha = alpha,
    n_samples_used = res$n_samples_used
  )
  res
}

# -------------------------------------------------------------------------
# 5. DOUBLE-RING LAYOUT FOR COMBINED NETWORK PLOT
# -------------------------------------------------------------------------
create_double_ring_layout <- function(net, inner_radius = 0.6, outer_radius = 1.0) {
  bac_idx <- which(V(net)$source == "16S")
  euk_idx <- which(V(net)$source == "18S")
  
  layout_matrix <- matrix(0, nrow = vcount(net), ncol = 2)
  
  if (length(bac_idx) > 0) {
    ord <- bac_idx[order(V(net)$phylum[bac_idx])]
    ang <- seq(0, 2 * pi, length.out = length(bac_idx) + 1)[-(length(bac_idx) + 1)]
    layout_matrix[ord, 1] <- inner_radius * cos(ang)
    layout_matrix[ord, 2] <- inner_radius * sin(ang)
  }
  if (length(euk_idx) > 0) {
    ord <- euk_idx[order(V(net)$phylum[euk_idx])]
    ang <- seq(0, 2 * pi, length.out = length(euk_idx) + 1)[-(length(euk_idx) + 1)]
    layout_matrix[ord, 1] <- outer_radius * cos(ang)
    layout_matrix[ord, 2] <- outer_radius * sin(ang)
  }
  layout_matrix
}

# -------------------------------------------------------------------------
# 6. PLOTTING — SINGLE-HABITAT COMBINED NETWORK (kept for completeness;
#    not used directly by the final row figure below, but process_environment()
#    calls it as a side effect to save per-habitat individual plots)
# -------------------------------------------------------------------------
plot_combined_network <- function(net_result, environment, output_dir) {
  
  if (is.null(net_result) || vcount(net_result$network) == 0) {
    cat("No combined network to plot for", environment, "\n")
    return(invisible(NULL))
  }
  
  net <- net_result$network
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  layout_matrix <- create_double_ring_layout(net)
  
  pdf(file.path(plot_dir, paste0(environment, "_combined_network_CLR_FDR.pdf")), width = 14, height = 10)
  layout(matrix(c(1, 2), ncol = 2), widths = c(2.5, 1))
  
  par(mar = c(1, 1, 3, 0))
  plot(net,
       vertex.color = V(net)$color,
       vertex.size = pmax(degree(net) * 0.5, 3),
       vertex.shape = ifelse(V(net)$source == "16S", "circle", "square"),
       vertex.frame.color = "black",
       vertex.label = NA,
       edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
       edge.curved = 0.3,
       edge.width = abs(E(net)$weight) * 2,
       layout = layout_matrix,
       xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), rescale = FALSE, asp = 1,
       main = paste(environment, "Environment: 16S-18S Network (CLR, |r|>=",
                    net_result$stats$cor_cutoff, ", FDR-adj p<", net_result$stats$fdr_alpha, ")"))
  
  text(0.7, 0, "16S", font = 2, cex = 1.2)
  text(1.15, 0, "18S", font = 2, cex = 1.2)
  
  par(mar = c(1, 0, 3, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Legend", cex.main = 1.2)
  
  bac_phyla_in_net <- sort(unique(V(net)$phylum[V(net)$source == "16S"]))
  euk_phyla_in_net <- sort(unique(V(net)$phylum[V(net)$source == "18S"]))
  
  n_legend_rows <- length(bac_phyla_in_net) + length(euk_phyla_in_net) + 4
  y_top <- 0.95
  line_height <- min(1 / n_legend_rows, 0.04)
  legend_cex <- ifelse(n_legend_rows > 25, 0.7, 0.9)
  
  if (length(bac_phyla_in_net) > 0) {
    legend(0.02, y_top, legend = bac_phyla_in_net, pch = 21, pt.bg = GLOBAL_16S_COLORS[bac_phyla_in_net],
           col = "black", pt.cex = 1.3, cex = legend_cex, title = "16S Taxa", title.font = 2, bty = "n")
    y_top <- max(y_top - (length(bac_phyla_in_net) + 2) * line_height, 0.15)
  }
  if (length(euk_phyla_in_net) > 0) {
    legend(0.02, y_top, legend = euk_phyla_in_net, pch = 22, pt.bg = GLOBAL_18S_COLORS[euk_phyla_in_net],
           col = "black", pt.cex = 1.3, cex = legend_cex, title = "18S Taxa", title.font = 2, bty = "n")
    y_top <- max(y_top - (length(euk_phyla_in_net) + 2) * line_height, 0.08)
  }
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 2, cex = legend_cex, title = "Correlations",
         title.font = 2, bty = "n")
  
  dev.off()
  cat("Plot saved:", file.path(plot_dir, paste0(environment, "_combined_network_CLR_FDR.pdf")), "\n")
}

# -------------------------------------------------------------------------
# 7. MAIN PIPELINE — per habitat: prevalence filter, top-N, CLR, build, plot
# -------------------------------------------------------------------------
process_environment <- function(env_name, input_16s, input_18s,
                                top_n_16s = 200, top_n_18s = 100,
                                top_n_combined = 75,
                                cor_cutoff = 0.8, alpha = 0.05) {
  
  cat("\n=== PROCESSING", toupper(env_name), "===\n")
  
  env_16s <- filter_data(input_16s, 'Environment', keep_vals = env_name)
  env_18s <- filter_data(input_18s, 'Environment', keep_vals = env_name)
  
  common_samples <- intersect(colnames(env_16s$data_loaded), colnames(env_18s$data_loaded))
  cat("Common samples:", length(common_samples), "\n")
  if (length(common_samples) < 5) {
    cat("WARNING: insufficient common samples — skipping\n")
    return(NULL)
  }
  
  bac_counts <- env_16s$data_loaded[, common_samples, drop = FALSE]
  euk_counts <- env_18s$data_loaded[, common_samples, drop = FALSE]
  
  filter_and_topn <- function(counts, top_n) {
    rel <- sweep(counts, 2, colSums(counts), "/")
    min_s <- ceiling(0.2 * ncol(counts))
    keep <- rowSums(rel > 0.0001) >= min_s
    counts_f <- counts[keep, , drop = FALSE]
    if (nrow(counts_f) > top_n) {
      top_ids <- names(sort(rowSums(counts_f), decreasing = TRUE))[1:top_n]
      counts_f <- counts_f[top_ids, , drop = FALSE]
    }
    counts_f
  }
  
  bac_final <- filter_and_topn(bac_counts, top_n_16s)
  euk_final <- filter_and_topn(euk_counts, top_n_18s)
  cat("16S ASVs (single):", nrow(bac_final), "| 18S ASVs (single):", nrow(euk_final), "\n")
  
  bac_clr <- clr_transform_data(bac_final)
  euk_clr <- clr_transform_data(euk_final)
  
  results <- list()
  
  if (nrow(bac_final) > 5) {
    results$net_16s <- build_single_network(bac_clr, input_16s$taxonomy_loaded, GLOBAL_16S_COLORS,
                                            cor_cutoff, alpha)
  }
  if (nrow(euk_final) > 5) {
    results$net_18s <- build_single_network(euk_clr, input_18s$taxonomy_loaded, GLOBAL_18S_COLORS,
                                            cor_cutoff, alpha)
  }
  
  bac_final_comb <- filter_and_topn(bac_counts, top_n_combined)
  euk_final_comb <- filter_and_topn(euk_counts, top_n_combined)
  cat("16S ASVs (combined):", nrow(bac_final_comb), "| 18S ASVs (combined):", nrow(euk_final_comb), "\n")
  
  if (nrow(bac_final_comb) > 5 && nrow(euk_final_comb) > 5) {
    bac_clr_comb <- clr_transform_data(bac_final_comb)
    euk_clr_comb <- clr_transform_data(euk_final_comb)
    combined_clr <- rbind(bac_clr_comb, euk_clr_comb)
    
    bac_tax <- input_16s$taxonomy_loaded[rownames(bac_clr_comb), , drop = FALSE]
    euk_tax <- input_18s$taxonomy_loaded[rownames(euk_clr_comb), , drop = FALSE]
    
    bac_tax$taxonomy2 <- as.character(bac_tax$taxonomy2)
    bac_tax$taxonomy2[is.na(bac_tax$taxonomy2) | bac_tax$taxonomy2 == "NA"] <- "Unclassified"
    euk_tax$taxonomy2 <- as.character(euk_tax$taxonomy2)
    euk_tax$taxonomy2[is.na(euk_tax$taxonomy2) | euk_tax$taxonomy2 == "NA"] <- "Unclassified"
    
    bac_tax$source <- "16S"
    euk_tax$source <- "18S"
    combined_tax <- rbind(bac_tax[, c("taxonomy2", "source")], euk_tax[, c("taxonomy2", "source")])
    
    results$net_combined <- build_combined_network(combined_clr, combined_tax, cor_cutoff, alpha)
    plot_combined_network(results$net_combined, env_name, output_dir)
  }
  
  if (!is.null(results$net_16s)) { cat("16S network:\n"); print(results$net_16s$stats) }
  if (!is.null(results$net_18s)) { cat("18S network:\n"); print(results$net_18s$stats) }
  if (!is.null(results$net_combined)) { cat("Combined network:\n"); print(results$net_combined$stats) }
  
  results
}

# -------------------------------------------------------------------------
# 8. RUN FOR ALL HABITATS (Desiccated -> Pond -> Lake)
# -------------------------------------------------------------------------
environments <- c("Desiccated", "Pond", "Lake")
all_results <- list()
for (env in environments) {
  all_results[[env]] <- process_environment(env, input_filt_16s, input_filt_18s,
                                            top_n_16s = 200, top_n_18s = 100,
                                            top_n_combined = 75,
                                            cor_cutoff = 0.8, alpha = 0.05)
}

# -------------------------------------------------------------------------
# 9. FINAL COMBINED ROW FIGURE — Desiccated | Pond | Lake
#    Fixes: ring labels closer to rings; "Legend" title removed; legend
#    taxon-name text enlarged; 16S/18S legend blocks closer together.
# -------------------------------------------------------------------------
plot_combined_row_figure_final <- function(all_results, output_dir,
                                           env_order = c("Desiccated", "Pond", "Lake"),
                                           inner_radius = 0.6, outer_radius = 1.0) {
  
  valid_envs <- env_order[sapply(env_order, function(e) {
    !is.null(all_results[[e]]) &&
      !is.null(all_results[[e]]$net_combined) &&
      !is.null(all_results[[e]]$net_combined$network) &&
      vcount(all_results[[e]]$net_combined$network) > 0
  })]
  
  if (length(valid_envs) == 0) {
    cat("No valid combined networks available to plot.\n")
    return(invisible(NULL))
  }
  
  plot_dir <- file.path(output_dir, "network_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  all_bac_phyla <- character(0)
  all_euk_phyla <- character(0)
  for (env in valid_envs) {
    net <- all_results[[env]]$net_combined$network
    all_bac_phyla <- union(all_bac_phyla, unique(V(net)$phylum[V(net)$source == "16S"]))
    all_euk_phyla <- union(all_euk_phyla, unique(V(net)$phylum[V(net)$source == "18S"]))
  }
  all_bac_phyla <- sort(all_bac_phyla)
  all_euk_phyla <- sort(all_euk_phyla)
  
  n_legend_rows <- length(all_bac_phyla) + length(all_euk_phyla) + 4
  fig_height <- 6.5
  
  pdf(file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row_FINAL.pdf"),
      width = 6 * length(valid_envs) + 3.8, height = fig_height)
  
  layout(matrix(1:(length(valid_envs) + 1), nrow = 1),
         widths = c(rep(5, length(valid_envs)), 2.8))
  
  # Precise 1cm title-to-plot-box gap, computed from this device's actual
  # line height (par("csi"), in inches) rather than a guessed multiplier.
  # Top margin sized to gap + enough lines for the title text itself
  # (cex.main=2.2 needs ~2.5 lines of height) so it is never clipped.
  gap_cm      <- 1
  gap_lines   <- (gap_cm / 2.54) / par("csi")
  top_margin  <- gap_lines + 2.5
  
  for (env in valid_envs) {
    net_result <- all_results[[env]]$net_combined
    net <- net_result$network
    layout_matrix <- create_double_ring_layout(net, inner_radius, outer_radius)
    
    par(mar = c(1, 1, top_margin, 1))
    plot(net,
         vertex.color = V(net)$color,
         vertex.size = pmax(degree(net) * 0.5, 3),
         vertex.shape = ifelse(V(net)$source == "16S", "circle", "square"),
         vertex.frame.color = "black",
         vertex.label = NA,
         edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
         edge.curved = 0.3,
         edge.width = abs(E(net)$weight) * 2,
         layout = layout_matrix,
         main = "",
         xlim = c(-outer_radius * 1.18, outer_radius * 1.18),
         ylim = c(-outer_radius * 1.18, outer_radius * 1.18),
         rescale = FALSE, asp = 1)
    
    title(main = env, cex.main = 2.2, font.main = 2, line = gap_lines)
    
    # Ring labels — pushed further right; xlim widened above so they are
    # never clipped at this distance
    text(inner_radius * 1.18, 0, "16S", font = 2, cex = 2, col = "black")
    text(outer_radius * 1.14, 0, "18S", font = 2, cex = 2, col = "black")
  }
  
  par(mar = c(1, 0, top_margin, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  # "Legend" title intentionally removed
  
  # Legend block spacing uses the ACTUAL bounding box returned by each
  # legend() call ($rect$top, $rect$h) instead of an estimated line_height —
  # this guarantees no overlap regardless of text size or entry count
  y_top <- 0.98
  
  if (length(all_bac_phyla) > 0) {
    l1 <- legend(0.02, y_top, legend = all_bac_phyla, pch = 21, pt.bg = GLOBAL_16S_COLORS[all_bac_phyla],
                 col = "black", pt.cex = 1.9, cex = 1.35, title = "16S Taxa", ncol = 2,
                 title.font = 2, bty = "n")
    y_top <- l1$rect$top - l1$rect$h - 0.03
  }
  
  if (length(all_euk_phyla) > 0) {
    l2 <- legend(0.02, y_top, legend = all_euk_phyla, pch = 22, pt.bg = GLOBAL_18S_COLORS[all_euk_phyla],
                 col = "black", pt.cex = 1.9, cex = 1.35, title = "18S Taxa", ncol = 2,
                 title.font = 2, bty = "n")
    y_top <- l2$rect$top - l2$rect$h - 0.03
  }
  
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 3, cex = 1.2, title = "Correlations",
         title.font = 2, bty = "n")
  
  dev.off()
  cat("Final combined row figure saved:",
      file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row_FINAL.pdf"), "\n")
  
  # -----------------------------------------------------------------------
  # Also save as high-res JPEG (re-runs the same plotting code into a
  # different graphics device, since this is base-R plotting, not ggplot2)
  # -----------------------------------------------------------------------
  jpeg(file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row_FINAL.jpeg"),
       width = 6 * length(valid_envs) + 3.8, height = fig_height,
       units = "in", res = 300, quality = 100)
  
  layout(matrix(1:(length(valid_envs) + 1), nrow = 1),
         widths = c(rep(5, length(valid_envs)), 2.8))
  
  # Recompute the precise 1cm gap for THIS device — par("csi") can differ
  # between the PDF (vector) and JPEG (raster, 300dpi) devices
  gap_cm      <- 0.5
  gap_lines   <- (gap_cm / 2.54) / par("csi")
  top_margin  <- gap_lines + 2.5
  
  for (env in valid_envs) {
    net_result <- all_results[[env]]$net_combined
    net <- net_result$network
    layout_matrix <- create_double_ring_layout(net, inner_radius, outer_radius)
    
    par(mar = c(1, 1, top_margin, 1))
    plot(net,
         vertex.color = V(net)$color,
         vertex.size = pmax(degree(net) * 0.5, 3),
         vertex.shape = ifelse(V(net)$source == "16S", "circle", "square"),
         vertex.frame.color = "black",
         vertex.label = NA,
         edge.color = ifelse(E(net)$weight > 0, "#619CFF", "#F8766D"),
         edge.curved = 0.3,
         edge.width = abs(E(net)$weight) * 2,
         layout = layout_matrix,
         main = "",
         xlim = c(-outer_radius * 1.18, outer_radius * 1.18),
         ylim = c(-outer_radius * 1.18, outer_radius * 1.18),
         rescale = FALSE, asp = 1)
    
    title(main = env, cex.main = 2.2, font.main = 2, line = gap_lines)
    text(inner_radius * 1.18, 0, "16S", font = 2, cex = 2, col = "black")
    text(outer_radius * 1.14, 0, "18S", font = 2, cex = 2, col = "black")
  }
  
  par(mar = c(1, 0, top_margin, 1))
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  y_top <- 0.98
  
  if (length(all_bac_phyla) > 0) {
    l1 <- legend(0.02, y_top, legend = all_bac_phyla, pch = 21, pt.bg = GLOBAL_16S_COLORS[all_bac_phyla],
                 col = "black", pt.cex = 1.9, cex = 1.35, title = "16S Taxa", ncol = 2,
                 title.font = 2, title.cex = 1.6, bty = "n")
    y_top <- l1$rect$top - l1$rect$h - 0.03
  }
  
  if (length(all_euk_phyla) > 0) {
    l2 <- legend(0.02, y_top, legend = all_euk_phyla, pch = 22, pt.bg = GLOBAL_18S_COLORS[all_euk_phyla],
                 col = "black", pt.cex = 1.9, cex = 1.35, title = "18S Taxa", ncol = 2,
                 title.font = 2, title.cex = 1.6, bty = "n")
    y_top <- l2$rect$top - l2$rect$h - 0.03
  }
  
  legend(0.02, y_top, legend = c("Positive", "Negative"),
         col = c("#619CFF", "#F8766D"), lty = 1, lwd = 3, cex = 1.2,title = "Correlations",
         title.font = 3, bty = "n")
  
  dev.off()
  cat("Final combined row figure (JPEG) saved:",
      file.path(plot_dir, "Combined_16S_18S_Network_AllEnvironments_Row_FINAL.jpeg"), "\n")
}

plot_combined_row_figure_final(all_results, output_dir, env_order = c("Desiccated", "Pond", "Lake"))
