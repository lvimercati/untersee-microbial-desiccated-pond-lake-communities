# =============================================================================
# FEAST SOURCE-TRACKING ANALYSIS — runs all six sink/domain combinations
# (16S and 18S, each with Desiccated/Pond/Lake as sink), aggregates results
# by source habitat, and prints mean compositional overlap per comparison.
# Input: FEAST_counts_16S.csv, FEAST_metadata_16S.csv, FEAST_counts_18S.csv,
#        FEAST_metadata_18S.csv (see 08_feast_figure.R for how these are
#        generated from the rarefied ASV tables)
# =============================================================================
library(FEAST)

# =============================================================================
# LOAD BOTH DOMAINS (adjust paths if needed — should already exist from
# your first successful run)
# =============================================================================
counts_16s <- read.csv("FEAST_counts_16S.csv", row.names = 1, check.names = FALSE)
meta_16s   <- read.csv("FEAST_metadata_16S.csv")
otus_16s   <- t(counts_16s)

counts_18s <- read.csv("FEAST_counts_18S.csv", row.names = 1, check.names = FALSE)
meta_18s   <- read.csv("FEAST_metadata_18S.csv")
otus_18s   <- t(counts_18s)

# =============================================================================
# HELPER FUNCTION 1 — builds metadata for a given sink habitat, runs FEAST
# =============================================================================
run_feast_sink <- function(otus, meta, sink_env, domain_label) {
  
  metadata <- data.frame(
    Env = meta$Environment,
    SourceSink = ifelse(meta$Environment == sink_env, "Sink", "Source"),
    id = ifelse(meta$Environment == sink_env, meta$SampleID, NA)
  )
  rownames(metadata) <- meta$SampleID
  
  otus <- otus[rownames(metadata), ]
  stopifnot(identical(rownames(otus), rownames(metadata)))
  
  outfile_name <- paste0(domain_label, "_", tolower(sink_env))
  cat("\n=== Running FEAST:", outfile_name, "===\n")
  cat("Sink samples:", sum(metadata$SourceSink == "Sink"),
      "| Source samples:", sum(metadata$SourceSink == "Source"), "\n")
  
  FEAST(
    C = otus,
    metadata = metadata,
    different_sources_flag = 0,
    dir_path = "~/results/",
    outfile = outfile_name
  )
  
  # Read back the saved output file and aggregate by source category
  result_file <- paste0("~/results/", outfile_name, "_source_contributions_matrix.txt")
  result_matrix <- read.table(result_file, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Identify the OTHER two habitats (not the sink) as source categories
  other_envs <- setdiff(c("Desiccated", "Pond", "Lake"), sink_env)
  
  summary_df <- data.frame(SampleID = rownames(result_matrix))
  for (env in other_envs) {
    cols <- grep(paste0("_", env, "$"), colnames(result_matrix))
    summary_df[[env]] <- rowSums(result_matrix[, cols, drop = FALSE])
  }
  summary_df$Unknown <- result_matrix[, "Unknown"]
  
  cat("\nMean across", nrow(summary_df), domain_label, sink_env, "samples:\n")
  for (env in other_envs) {
    cat(" ", env, ":", round(mean(summary_df[[env]]) * 100, 1), "%\n")
  }
  cat("  Unknown:", round(mean(summary_df$Unknown) * 100, 1), "%\n")
  
  summary_df
}

# =============================================================================
# RUN ALL SIX
# =============================================================================
feast_16s_desiccated <- run_feast_sink(otus_16s, meta_16s, "Desiccated", "16S")
feast_16s_pond       <- run_feast_sink(otus_16s, meta_16s, "Pond", "16S")
feast_16s_lake       <- run_feast_sink(otus_16s, meta_16s, "Lake", "16S")

feast_18s_desiccated <- run_feast_sink(otus_18s, meta_18s, "Desiccated", "18S")
feast_18s_pond       <- run_feast_sink(otus_18s, meta_18s, "Pond", "18S")
feast_18s_lake       <- run_feast_sink(otus_18s, meta_18s, "Lake", "18S")

cat("\n\n=== ALL SIX FEAST RUNS COMPLETE ===\n")
