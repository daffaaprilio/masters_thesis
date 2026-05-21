#!/usr/bin/env Rscript
# DSS DMR analysis: SBC10 vs each other accession (TAA contrast)
#
# SBC10 is the high-TAA-producer (TAA +++, low secretion) and is always
# group1 so that diff.Methy > 0 means hypermethylated in SBC10.
#
# Results are exploratory — no biological replicates available (n=1 per group).
#
# Run via:
#   ./docker/run.sh Rscript analysis/scripts/run_dss_dmr_taa.R

suppressPackageStartupMessages({
  library(DSS)
  library(data.table)
})

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
DSS_DIR <- "analysis/data/taa_DSS"
OUT_DIR <- "analysis/data/taa_DMR"
LOG_DIR <- "analysis/logs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
log_path <- file.path(LOG_DIR, format(Sys.time(), "run_dss_dmr_taa_%Y%m%d_%H%M%S.log"))
log_con  <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)
cat(sprintf("Log: %s\n", log_path))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# --------------------------------------------------------------------------- #
# DSS parameters
# --------------------------------------------------------------------------- #
SMOOTH_SPAN   <- 500
DMR_DELTA     <- 0.1   # minimum methylation difference
DMR_P_THRESH  <- 0.2   # relaxed — justified by n=1 design
DMR_MINLEN    <- 50    # bp
DMR_MINCG     <- 20    # CpG sites per DMR
DMR_DIS_MERGE <- 300   # bp gap to merge neighbouring DMRs

# SBC10 is fixed as group1; comparators are group2
PAIRS <- list(
  c("SBC10", "SBC4"),
  c("SBC10", "SBC11"),
  c("SBC10", "SBC23")
)

# --------------------------------------------------------------------------- #
# Load all DSS input files
# --------------------------------------------------------------------------- #
cat("Loading 5mC DSS input files...\n")

samples_needed <- unique(unlist(PAIRS))
dss_data <- lapply(samples_needed, function(s) {
  path <- file.path(DSS_DIR, paste0(s, ".5mC.dss.txt"))
  if (!file.exists(path)) stop(sprintf("Missing DSS file: %s", path))
  df <- fread(path, header = TRUE, sep = "\t",
              colClasses = c(chr = "character", pos = "integer",
                             N   = "integer",   X   = "integer"))
  cat(sprintf("  %-8s %d sites\n", s, nrow(df)))
  as.data.frame(df)
})
names(dss_data) <- samples_needed

# --------------------------------------------------------------------------- #
# SBC10 vs others loop
# --------------------------------------------------------------------------- #
summary_rows <- list()
all_dmrs     <- list()

sep_line <- strrep("─", 60)

for (pair in PAIRS) {
  sA        <- pair[1]   # always SBC10
  sB        <- pair[2]
  pair_name <- paste0(sA, "_vs_", sB)

  cat(sprintf("\n%s\n  %s\n%s\n", sep_line, pair_name, sep_line))

  bs <- tryCatch(
    makeBSseqData(list(dss_data[[sA]], dss_data[[sB]]),
                  sampleNames = c(sA, sB)),
    error = function(e) {
      cat(sprintf("  ERROR building BSseq: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(bs)) next

  cat("  DMLtest (smoothing=TRUE, span=", SMOOTH_SPAN, ")...\n", sep = "")
  dml <- tryCatch(
    DMLtest(bs,
            group1         = 1L,
            group2         = 2L,
            smoothing      = TRUE,
            smoothing.span = SMOOTH_SPAN),
    error = function(e) {
      cat(sprintf("  ERROR in DMLtest: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dml)) next

  n_tested <- nrow(dml)
  pval_col <- intersect(c("pvalue", "pval"), names(dml))[1]
  n_sig    <- if (!is.na(pval_col)) sum(dml[[pval_col]] < DMR_P_THRESH, na.rm = TRUE) else NA
  cat(sprintf("  DMLs tested: %d  |  p < %.2f: %s\n",
              n_tested, DMR_P_THRESH,
              ifelse(is.na(n_sig), "n/a", format(n_sig, big.mark = ","))))

  dml_path <- file.path(OUT_DIR, paste0(pair_name, ".5mC.DML.tsv"))
  fwrite(as.data.table(dml), dml_path, sep = "\t")
  cat(sprintf("  DML table → %s\n", basename(dml_path)))

  cat(sprintf("  callDMR (delta=%.2f, p<%.2f, minlen=%d, minCG=%d, merge=%d)...\n",
              DMR_DELTA, DMR_P_THRESH, DMR_MINLEN, DMR_MINCG, DMR_DIS_MERGE))
  dmr <- tryCatch(
    callDMR(dml,
            delta       = DMR_DELTA,
            p.threshold = DMR_P_THRESH,
            minlen      = DMR_MINLEN,
            minCG       = DMR_MINCG,
            dis.merge   = DMR_DIS_MERGE),
    error = function(e) {
      cat(sprintf("  ERROR in callDMR: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(dmr) || nrow(dmr) == 0) {
    cat(sprintf("  *** No DMRs found for %s\n", pair_name))
    summary_rows[[pair_name]] <- data.frame(
      pair              = pair_name,
      n_DML_tested      = n_tested,
      n_DML_significant = ifelse(is.na(n_sig), NA_integer_, as.integer(n_sig)),
      n_DMR             = 0L,
      n_hyper_SBC10     = NA_integer_,
      n_hyper_other     = NA_integer_,
      mean_diff         = NA_real_,
      mean_length_bp    = NA_real_,
      mean_nCG          = NA_real_
    )
    next
  }

  n_dmr <- nrow(dmr)
  cat(sprintf("  DMRs called: %d\n", n_dmr))

  if (any(dmr$nCG < DMR_MINCG, na.rm = TRUE))
    cat(sprintf("  *** WARNING: %d DMR(s) have nCG < %d\n",
                sum(dmr$nCG < DMR_MINCG, na.rm = TRUE), DMR_MINCG))

  # diff.Methy > 0 → hypermethylated in SBC10 (group1)
  dmr$sample_a  <- sA
  dmr$sample_b  <- sB
  dmr$direction <- ifelse(dmr$diff.Methy > 0,
                          paste0("hyper_", sA),
                          paste0("hyper_", sB))

  mean_diff  <- mean(dmr$diff.Methy, na.rm = TRUE)
  diff_range <- range(dmr$diff.Methy, na.rm = TRUE)
  cat(sprintf("  diff.Methy: mean=%.3f  range=[%.3f, %.3f]\n",
              mean_diff, diff_range[1], diff_range[2]))

  want_cols <- c("chr", "start", "end", "length", "nCG",
                 "meanMethy1", "meanMethy2", "diff.Methy", "areaStat",
                 "sample_a", "sample_b", "direction")
  dmr_out <- dmr[, intersect(want_cols, names(dmr))]

  dmr_path <- file.path(OUT_DIR, paste0(pair_name, ".5mC.DMR.tsv"))
  fwrite(as.data.table(dmr_out), dmr_path, sep = "\t")
  cat(sprintf("  DMR table → %s\n", basename(dmr_path)))

  all_dmrs[[pair_name]] <- dmr_out

  n_hyper_sbc10 <- sum(dmr$direction == paste0("hyper_", sA), na.rm = TRUE)
  n_hyper_other <- sum(dmr$direction == paste0("hyper_", sB), na.rm = TRUE)

  summary_rows[[pair_name]] <- data.frame(
    pair              = pair_name,
    n_DML_tested      = n_tested,
    n_DML_significant = ifelse(is.na(n_sig), NA_integer_, as.integer(n_sig)),
    n_DMR             = n_dmr,
    n_hyper_SBC10     = n_hyper_sbc10,
    n_hyper_other     = n_hyper_other,
    mean_diff         = round(mean_diff, 4),
    mean_length_bp    = round(mean(dmr$length, na.rm = TRUE), 1),
    mean_nCG          = round(mean(dmr$nCG,    na.rm = TRUE), 1)
  )
}

# --------------------------------------------------------------------------- #
# Summary outputs
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
  combined <- rbindlist(lapply(all_dmrs, as.data.table), fill = TRUE)
  fwrite(combined, file.path(OUT_DIR, "DMR_all_combined.tsv"), sep = "\t")
  cat(sprintf("\n  DMR_all_combined.tsv written (%d total DMRs across %d pairs)\n",
              nrow(combined), length(all_dmrs)))
} else {
  cat("\n  No DMRs found in any pair.\n")
}

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Done.\n")
