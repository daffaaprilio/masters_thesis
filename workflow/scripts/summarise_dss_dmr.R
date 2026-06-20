#!/usr/bin/env Rscript
# Summarise DMR results across all pairwise comparisons.
#
# Run after all run_dss_dmr.R jobs have finished:
#   ./docker/run.sh Rscript workflow/scripts/summarise_dss_dmr.R

suppressPackageStartupMessages({
  library(data.table)
})

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
OUT_DIR <- "results/DMR"
LOG_DIR <- "analysis/logs"
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(LOG_DIR, format(Sys.time(), "summarise_dss_dmr_%Y%m%d_%H%M%S.log"))
log_con  <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)
cat(sprintf("Log: %s\n", log_path))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

PAIRS <- c("SBC10_vs_SBC4", "SBC10_vs_SBC11", "SBC10_vs_SBC23",
           "SBC11_vs_SBC4", "SBC11_vs_SBC23", "SBC4_vs_SBC23")

sep_line <- strrep("─", 60)

# --------------------------------------------------------------------------- #
# Load per-pair DMR files
# --------------------------------------------------------------------------- #
summary_rows <- list()
all_dmrs     <- list()

for (pair_name in PAIRS) {
  dmr_path <- file.path(OUT_DIR, paste0(pair_name, ".5mC.DMR.tsv"))

  if (!file.exists(dmr_path)) {
    cat(sprintf("*** Missing: %s — skipping\n", dmr_path))
    next
  }

  dmr <- fread(dmr_path)
  n_dmr <- nrow(dmr)
  cat(sprintf("%s\n  DMRs loaded: %d\n", pair_name, n_dmr))

  all_dmrs[[pair_name]] <- dmr

  if (n_dmr == 0) {
    summary_rows[[pair_name]] <- data.frame(
      pair           = pair_name,
      n_DMR          = 0L,
      n_hyper_a      = NA_integer_,
      n_hyper_b      = NA_integer_,
      mean_diff      = NA_real_,
      mean_length_bp = NA_real_,
      mean_nCG       = NA_real_
    )
    next
  }

  sA <- sub("_vs_.*", "", pair_name)
  sB <- sub(".*_vs_", "", pair_name)

  n_hyper_a <- sum(dmr$direction == paste0("hyper_", sA), na.rm = TRUE)
  n_hyper_b <- sum(dmr$direction == paste0("hyper_", sB), na.rm = TRUE)
  mean_diff <- mean(dmr$diff.Methy, na.rm = TRUE)

  cat(sprintf("  hyper_%s: %d  |  hyper_%s: %d  |  mean diff: %.3f\n",
              sA, n_hyper_a, sB, n_hyper_b, mean_diff))

  summary_rows[[pair_name]] <- data.frame(
    pair           = pair_name,
    n_DMR          = n_dmr,
    n_hyper_a      = n_hyper_a,
    n_hyper_b      = n_hyper_b,
    mean_diff      = round(mean_diff, 4),
    mean_length_bp = round(mean(dmr$length, na.rm = TRUE), 1),
    mean_nCG       = round(mean(dmr$nCG,    na.rm = TRUE), 1)
  )
}

# --------------------------------------------------------------------------- #
# Write combined outputs
# --------------------------------------------------------------------------- #
cat(sprintf("\n%s\n  Summary\n%s\n", sep_line, sep_line))

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  fwrite(summary_df, file.path(OUT_DIR, "DMR_summary.tsv"), sep = "\t")
  cat("  DMR_summary.tsv written\n\n")
  print(summary_df, row.names = FALSE)
}

if (length(all_dmrs) > 0) {
  combined <- rbindlist(all_dmrs, fill = TRUE)
  fwrite(combined, file.path(OUT_DIR, "DMR_all_combined.tsv"), sep = "\t")
  cat(sprintf("\n  DMR_all_combined.tsv written (%d total DMRs across %d pairs)\n",
              nrow(combined), length(all_dmrs)))
} else {
  cat("\n  No DMR files found.\n")
}

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Done.\n")
