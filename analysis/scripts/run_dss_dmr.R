#!/usr/bin/env Rscript
# Pairwise DSS DMR analysis on sorgoleone 5mC loci
# All six pairwise comparisons across four sorghum accessions (n=1 per group)
# Results are exploratory — no biological replicates available

suppressPackageStartupMessages({
  library(DSS)
  library(data.table)
})

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
DSS_DIR <- "analysis/data/sorgoleone_DSS"
OUT_DIR <- "analysis/data/sorgoleone_DMR"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Logging — tee all cat() output to a timestamped log file
# --------------------------------------------------------------------------- #
log_path <- file.path(OUT_DIR, format(Sys.time(), "run_dss_dmr_%Y%m%d_%H%M%S.log"))
log_con  <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)
cat(sprintf("Log: %s\n", log_path))
cat(sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

SAMPLES <- c("SBC4", "SBC10", "SBC11", "SBC23")
PAIRS   <- combn(SAMPLES, 2, simplify = FALSE)

# --------------------------------------------------------------------------- #
# DSS parameters
# --------------------------------------------------------------------------- #
SMOOTH_SPAN   <- 500
DMR_DELTA     <- 0.1
DMR_P_THRESH  <- 0.2
DMR_MINLEN    <- 50
DMR_MINCG     <- 20
DMR_DIS_MERGE <- 300

# --------------------------------------------------------------------------- #
# Load all DSS input files
# --------------------------------------------------------------------------- #
cat("Loading 5mC DSS input files...\n")

dss_data <- lapply(SAMPLES, function(s) {
  path <- file.path(DSS_DIR, paste0(s, ".5mC.dss.txt"))
  if (!file.exists(path)) stop(sprintf("Missing: %s", path))
  df <- fread(path, header = TRUE, sep = "\t",
              colClasses = c(chr = "character", pos = "integer",
                             N   = "integer",   X   = "integer"))
  cat(sprintf("  %-8s %d sites\n", s, nrow(df)))
  as.data.frame(df)
})
names(dss_data) <- SAMPLES

# --------------------------------------------------------------------------- #
# Pairwise loop
# --------------------------------------------------------------------------- #
summary_rows <- list()
all_dmrs     <- list()

sep_line <- strrep("─", 60)

for (pair in PAIRS) {
  sA        <- pair[1]
  sB        <- pair[2]
  pair_name <- paste0(sA, "_vs_", sB)

  cat(sprintf("\n%s\n  %s\n%s\n", sep_line, pair_name, sep_line))

  # Build BSseq object
  bs <- tryCatch(
    makeBSseqData(list(dss_data[[sA]], dss_data[[sB]]),
                  sampleNames = c(sA, sB)),
    error = function(e) {
      cat(sprintf("  ERROR building BSseq: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(bs)) next

  # DML test (smoothed)
  cat("  DMLtest (smoothing=TRUE, span=", SMOOTH_SPAN, ")...\n", sep = "")
  dml <- tryCatch(
    DMLtest(bs,
            group1          = 1L,
            group2          = 2L,
            smoothing       = TRUE,
            smoothing.span  = SMOOTH_SPAN),
    error = function(e) {
      cat(sprintf("  ERROR in DMLtest: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dml)) next

  n_tested <- nrow(dml)
  # DSS returns pvalue column; guard against version naming differences
  pval_col  <- intersect(c("pvalue", "pval"), names(dml))[1]
  n_sig     <- if (!is.na(pval_col)) sum(dml[[pval_col]] < DMR_P_THRESH, na.rm = TRUE) else NA
  cat(sprintf("  DMLs tested: %d  |  p < %.2f: %s\n",
              n_tested, DMR_P_THRESH, ifelse(is.na(n_sig), "n/a", n_sig)))

  # Write full DML table
  dml_path <- file.path(OUT_DIR, paste0(pair_name, ".5mC.DML.tsv"))
  fwrite(as.data.table(dml), dml_path, sep = "\t")

  # Call DMRs
  cat(sprintf("  callDMR (delta=%.1f, p<%.2f, minlen=%d, minCG=%d)...\n",
              DMR_DELTA, DMR_P_THRESH, DMR_MINLEN, DMR_MINCG))
  dmr <- tryCatch(
    callDMR(dml,
            delta      = DMR_DELTA,
            p.threshold = DMR_P_THRESH,
            minlen     = DMR_MINLEN,
            minCG      = DMR_MINCG,
            dis.merge  = DMR_DIS_MERGE),
    error = function(e) {
      cat(sprintf("  ERROR in callDMR: %s — skipping\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dmr)) {
    cat(sprintf("  *** WARNING: DMR calling failed for %s\n", pair_name))
    next
  }

  n_dmr <- nrow(dmr)
  cat(sprintf("  DMRs called: %d\n", n_dmr))

  if (n_dmr == 0) {
    cat(sprintf("  *** WARNING: No DMRs for %s — check coverage at target loci\n", pair_name))
    summary_rows[[pair_name]] <- data.frame(
      pair              = pair_name,
      n_DML_tested      = n_tested,
      n_DML_significant = ifelse(is.na(n_sig), NA_integer_, n_sig),
      n_DMR             = 0L,
      n_hyper_A         = NA_integer_,
      n_hyper_B         = NA_integer_,
      mean_diff         = NA_real_,
      mean_length       = NA_real_,
      mean_nCG          = NA_real_
    )
    next
  }

  # Sanity checks
  if (any(dmr$nCG < DMR_MINCG, na.rm = TRUE))
    cat(sprintf("  *** WARNING: %d DMR(s) have nCG < %d\n",
                sum(dmr$nCG < DMR_MINCG, na.rm = TRUE), DMR_MINCG))

  # Annotate direction and identity
  dmr$sample_a  <- sA
  dmr$sample_b  <- sB
  dmr$direction <- ifelse(dmr$diff.Methy > 0,
                          "hypermethylated_in_A",
                          "hypermethylated_in_B")

  mean_diff  <- mean(dmr$diff.Methy, na.rm = TRUE)
  diff_range <- range(dmr$diff.Methy, na.rm = TRUE)
  cat(sprintf("  diff.Methy: mean=%.3f  range=[%.3f, %.3f]\n",
              mean_diff, diff_range[1], diff_range[2]))

  # Select output columns (guard against DSS version differences)
  want_cols <- c("chr", "start", "end", "length", "nCG",
                 "meanMethy1", "meanMethy2", "diff.Methy", "areaStat",
                 "sample_a", "sample_b", "direction")
  dmr_out <- dmr[, intersect(want_cols, names(dmr))]

  dmr_path <- file.path(OUT_DIR, paste0(pair_name, ".5mC.DMR.tsv"))
  fwrite(as.data.table(dmr_out), dmr_path, sep = "\t")

  all_dmrs[[pair_name]] <- dmr_out

  n_hyper_A <- sum(dmr$direction == "hypermethylated_in_A", na.rm = TRUE)
  n_hyper_B <- sum(dmr$direction == "hypermethylated_in_B", na.rm = TRUE)

  summary_rows[[pair_name]] <- data.frame(
    pair              = pair_name,
    n_DML_tested      = n_tested,
    n_DML_significant = ifelse(is.na(n_sig), NA_integer_, as.integer(n_sig)),
    n_DMR             = n_dmr,
    n_hyper_A         = n_hyper_A,
    n_hyper_B         = n_hyper_B,
    mean_diff         = round(mean_diff, 4),
    mean_length       = round(mean(dmr$length,  na.rm = TRUE), 1),
    mean_nCG          = round(mean(dmr$nCG,     na.rm = TRUE), 1)
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
} else {
  cat("  No pairs completed successfully.\n")
}

if (length(all_dmrs) > 0) {
  combined <- rbindlist(lapply(all_dmrs, as.data.table), fill = TRUE)
  fwrite(combined, file.path(OUT_DIR, "DMR_all_pairs_combined.tsv"), sep = "\t")
  cat(sprintf("\n  DMR_all_pairs_combined.tsv written (%d total DMRs across %d pairs)\n",
              nrow(combined), length(all_dmrs)))
} else {
  cat("\n  No DMRs found in any pair — combined file not written.\n")
}

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("Done.\n")
