# =============================================================================
# FEAST SOURCE CONTRIBUTION FIGURE (2x3): A) 16S (Bacteria), B) 18S (Eukarya)
# Standalone, self-contained. Reads the six FEAST "_source_contributions_
# matrix.txt" output files (transferred from the HPC), aggregates by source
# habitat, and builds a pie-chart figure matching the reference SourceTracker
# figure's style: same title/subtitle format, same colour scheme, same
# 2x3 (Desiccated | Pond | Lake) x (16S | 18S) layout.
# Outputs both PDF and high-res JPEG.
# =============================================================================

# -----------------------------------------------------------------------
# Set this to wherever you saved the six FEAST result files on your laptop
# -----------------------------------------------------------------------
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_Sourcetracker/FEAST/results/")

# -----------------------------------------------------------------------
# Shared colour palette — matches the reference figure and the rest of
# the manuscript's habitat colour convention
# -----------------------------------------------------------------------
habitat_colours <- c("Desiccated" = "#E31A1C", "Pond" = "#1F78B4",
                      "Lake" = "#33A02C", "Unknown" = "grey80")

# -----------------------------------------------------------------------
# FUNCTION: read a FEAST result file, aggregate contributions by source
# habitat (summing across all individual source samples within each
# category), and return mean percentages + sample count
# -----------------------------------------------------------------------
load_feast_summary <- function(filepath, sink_env) {
  result_matrix <- read.table(filepath, header = TRUE, sep = "\t", check.names = FALSE)
  other_envs <- setdiff(c("Desiccated", "Pond", "Lake"), sink_env)

  means <- c()
  for (env in other_envs) {
    cols <- grep(paste0("_", env, "$"), colnames(result_matrix))
    means[env] <- mean(rowSums(result_matrix[, cols, drop = FALSE]))
  }
  means["Unknown"] <- mean(result_matrix[, "Unknown"])

  list(means = means, n = nrow(result_matrix))
}

# -----------------------------------------------------------------------
# FUNCTION: draw one pie chart panel, matching the reference figure style
# (title + "Based on N samples" subtitle, percentage labels inside
# slices, external "Source Environment" legend)
# -----------------------------------------------------------------------
draw_source_pie <- function(means, n, sink_env) {
  # Order slices consistently: Desiccated, Pond, Lake, Unknown (whichever present)
  order_ref <- c("Desiccated", "Pond", "Lake", "Unknown")
  means <- means[order_ref[order_ref %in% names(means)]]

  cols <- habitat_colours[names(means)]
  labels <- paste0(round(means * 100, 1), "%")

  # Bold, larger percentage labels: pie()'s internal text() calls respect
  # the current par(font)/par(cex) settings, so set them beforehand
  old_par <- par(font = 2)
  pie(means, labels = labels, col = cols, cex = 1.4, main = "")
  par(old_par)

  title(main = paste0("Average Source Contributions to ", sink_env, " Samples"),
        cex.main = 1.3, font.main = 2, line = 2.6)
  title(main = paste0("Based on ", n, " samples"),
        cex.main = 1.0, font.main = 1, col.main = "grey40", line = 1.0)

  legend("topright", inset = c(-0.12, 0), legend = names(means),
         fill = cols, title = "Source Environment", title.font = 2,
         cex = 1.1, bty = "n", xpd = TRUE)
}

# =============================================================================
# LOAD ALL SIX FEAST RESULTS
# =============================================================================
res_16s_desic <- load_feast_summary("16S_desiccated_source_contributions_matrix.txt", "Desiccated")
res_16s_pond  <- load_feast_summary("16S_pond_source_contributions_matrix.txt", "Pond")
res_16s_lake  <- load_feast_summary("16S_lake_source_contributions_matrix.txt", "Lake")

res_18s_desic <- load_feast_summary("18S_desiccated_source_contributions_matrix.txt", "Desiccated")
res_18s_pond  <- load_feast_summary("18S_pond_source_contributions_matrix.txt", "Pond")
res_18s_lake  <- load_feast_summary("18S_lake_source_contributions_matrix.txt", "Lake")

# =============================================================================
# BUILD THE 2x3 FIGURE
# =============================================================================
setwd("~/Desktop/Research_Projects/Lake_Untersee/Dry_mats/No tufts/no_tufts_Sourcetracker/FEAST/")

build_figure <- function() {
  par(mfrow = c(2, 3), mar = c(1, 5, 4, 9), oma = c(1, 4, 0, 0))

  draw_source_pie(res_16s_desic$means, res_16s_desic$n, "Desiccated")
  draw_source_pie(res_16s_pond$means,  res_16s_pond$n,  "Pond")
  draw_source_pie(res_16s_lake$means,  res_16s_lake$n,  "Lake")

  draw_source_pie(res_18s_desic$means, res_18s_desic$n, "Desiccated")
  draw_source_pie(res_18s_pond$means,  res_18s_pond$n,  "Pond")
  draw_source_pie(res_18s_lake$means,  res_18s_lake$n,  "Lake")

  mtext("A) Bacterial (16S rRNA gene)", side = 2, line = 1, outer = TRUE,
        at = 0.75, cex = 1.4, font = 2)
  mtext("B) Eukaryotic (18S rRNA gene)", side = 2, line = 1, outer = TRUE,
        at = 0.25, cex = 1.4, font = 2)
}

pdf("Figure_FEAST_SourceContributions.pdf", width = 15, height = 9)
build_figure()
dev.off()

jpeg("Figure_FEAST_SourceContributions.jpeg", width = 15, height = 9,
     units = "in", res = 300, quality = 100)
build_figure()
dev.off()

cat("Saved: Figure_FEAST_SourceContributions.pdf and .jpeg\n")
