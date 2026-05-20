# Methylation Landscape EDA
Covers:
-   Modification type across each sample, i.e., 5mC, 6mA, 4mC (stacked bar graph)
-   Methylation context, i.e., CpG, CHG, CHH (stacked bar graph)
-   Detailed status of CpG 5mC over the contig (each sample, similar to the read depth plot)

## Running scripts for visualization

1.  Modification type, `modkit_summary_samples.sh`: Run `modkit summary` on each per-sample BAM in `resources/align_bam_sample/` to confirm modification types present (5mC, 6mA) and collect per-base methylation statistics. Outputs a timestamped TSV (`analysis/data/modkit_summary/modkit_summary_*.tsv`) with columns: sample, base, strand, n_called, n_mod, pct_modified, n_canonical, n_other, n_delete, n_fail, n_diff, n_nocall. <br>
    ```shell
    ./docker/run.sh bash analysis/scripts/modkit_summary_samples.sh
    ```