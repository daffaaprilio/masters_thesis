#!/usr/bin/env Rscript
# DSS DMR analysis: one pairwise comparison between two accessions.
#
# For a pair sA_vs_sB, sA is group1 so that diff.Methy > 0 means
# hypermethylated in sA (sample_a).
#
# Results are exploratory — no biological replicates available (n=1 per group).
#
# Run one comparison at a time (can be launched in parallel):
#   ./docker/run.sh Rscript workflow/scripts/run_dss_dmr.R SBC10_vs_SBC4
#   ./docker/run.sh Rscript workflow/scripts/run_dss_dmr.R SBC11_vs_SBC23
#   ./docker/run.sh Rscript workflow/scripts/run_dss_dmr.R SBC4_vs_SBC23

suppressPackageStartupMessages({
  library(DSS)
  library(data.table)
})

# --------------------------------------------------------------------------- #
# Pair selection via command-line argument  (must come before logging)
# --------------------------------------------------------------------------- #
ALL_PAIRS <- list(
  SBC10_vs_SBC4  = c("SBC10", "SBC4"),
  SBC10_vs_SBC11 = c("SBC10", "SBC11"),
  SBC10_vs_SBC23 = c("SBC10", "SBC23"),
  SBC11_vs_SBC4  = c("SBC11", "SBC4"),
  SBC11_vs_SBC23 = c("SBC11", "SBC23"),
  SBC4_vs_SBC23  = c("SBC4",  "SBC23")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript run_dss_dmr.R <pair>\n",
       "  where <pair> is one of: ", paste(names(ALL_PAIRS), collapse = ", "))
}
if (!args[1] %in% names(ALL_PAIRS)) {
  stop("Unknown pair '", args[1], "'. Choose from: ",
       paste(names(ALL_PAIRS), collapse = ", "))
}

PAIR_NAME <- args[1]

# Only the 10 main sorghum chromosomes — excludes unplaced scaffolds (NW_*)
# which would inflate runtime without contributing to gene-level analysis.
MAIN_CHROMS <- paste0("NC_01287", 0:9, ".2")

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
DSS_DIR <- "results/DSS"
OUT_DIR <- "results/DMR"
LOG_DIR <- "analysis/logs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Logging
# --------------------------------------------------------------------------- #
log_path <- file.path(LOG_DIR, sprintf("run_dss_dmr_%s_%s.log",
                                       PAIR_NAME, format(Sys.time(), "%Y%m%d_%H%M%S")))
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

PAIRS <- ALL_PAIRS[PAIR_NAME]
cat(sprintf("Running comparison: %s\n\n", PAIR_NAME))

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
  df <- df[df$chr %in% MAIN_CHROMS, ]
  cat(sprintf("  %-8s %d sites (main chromosomes only)\n", s, nrow(df)))
  as.data.frame(df)
})
names(dss_data) <- samples_needed

# --------------------------------------------------------------------------- #
# Pairwise DMR loop
# --------------------------------------------------------------------------- #
sep_line <- strrep("─", 60)

for (pair in PAIRS) {
  sA        <- pair[1]   # group1 (sample_a)
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
    next
  }

  n_dmr <- nrow(dmr)
  cat(sprintf("  DMRs called: %d\n", n_dmr))

  if (any(dmr$nCG < DMR_MINCG, na.rm = TRUE))
    cat(sprintf("  *** WARNING: %d DMR(s) have nCG < %d\n",
                sum(dmr$nCG < DMR_MINCG, na.rm = TRUE), DMR_MINCG))

  # diff.Methy > 0 → hypermethylated in group1 (sample_a)
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

}

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Done. Run summarise_dss_dmr.R once all comparisons complete.\n"))
