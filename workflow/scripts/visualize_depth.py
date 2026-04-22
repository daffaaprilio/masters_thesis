#!/home/daffa/Work/2026/thesis/.venv/bin/python

import argparse
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Visualize coverage depth per chromosome.")
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--input", help="Path to the .depth file")
    source_group.add_argument("--load-pickle", help="Path to load a DataFrame pickle from a previous run")
    parser.add_argument("--output", required=True, help="Path to the output plot image")
    parser.add_argument("--library", required=True, help="Library name (used as plot title)")
    parser.add_argument("--save-pickle", default=None, help="Optional path to save the DataFrame as a pickle file")
    parser.add_argument("--bin-size", type=int, default=10_000, help="Bin size in bp for averaging depth (default: 10000)")
    args = parser.parse_args()

    if args.load_pickle:
        df = load_df_pickle(args.load_pickle)
    else:
        df = read_depth_file(args.input)

    if args.save_pickle:
        save_df_pickle(df, args.save_pickle)

    visualize_depth(df, args.library, args.output, args.bin_size)

def read_depth_file(file_path):
    print(f"Reading {file_path} as pandas DataFrame")
    df = pd.read_csv(file_path, sep="\t", header=None,
                     names=["chrom", "pos", "depth"],
                     dtype={"chrom": str, "pos": int, "depth": int})
    return df

def save_df_pickle(df, output_path):
    df.to_pickle(output_path)
    print(f"DataFrame saved to {output_path}")

def load_df_pickle(pickle_path):
    print(f"Loading DataFrame from {pickle_path}")
    return pd.read_pickle(pickle_path)

def calc_coverage_statistics(df):
    return df["depth"].mean(), df["depth"].median(), df["depth"].std() 

def calc_coverage_statistics_chr(df):
    valid_contigs = [f"NC_01287{n}.2" for n in range(0, 10)]
    df_valid_contigs = df[df["chrom"].isin(valid_contigs)].copy()
    return df_valid_contigs.groupby("chrom")["depth"].mean(), df_valid_contigs.groupby("chrom").median(), df_valid_contigs.groupby("chrom").std()

def visualize_depth(df, library, filename, bin_size=10_000):
    print(f"Plotting depth per base for {library} (bin size: {bin_size} bp)")
    valid_contigs = [f"NC_01287{n}.2" for n in range(0, 10)]
    df_valid = df[df["chrom"].isin(valid_contigs)].copy()

    present = [c for c in valid_contigs if c in df_valid["chrom"].values]
    ncols = 2
    nrows = (len(present) + ncols - 1) // ncols

    # reserve bottom space for the global stats box
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 3 + 1.5), sharex=True, sharey=True)
    axes = axes.flatten()

    # global statistics (exclusive to 10 contigs/chromosomes) (whole-genome row)
    g_mean, g_median, g_std = calc_coverage_statistics(df_valid)

    # global statistics (include all contigs)
    ga_mean, ga_median, ga_std = calc_coverage_statistics(df)

    # per-chromosome statistics
    chr_mean, chr_median, chr_std = calc_coverage_statistics_chr(df_valid)

    for i, chrom in enumerate(present):
        ax = axes[i]
        chrom_df = df_valid[df_valid["chrom"] == chrom].copy()

        chrom_df["bin"] = (chrom_df["pos"] // bin_size) * bin_size
        binned = chrom_df.groupby("bin")["depth"].mean().reset_index()

        ax.scatter(binned["bin"] / 1_000_000, binned["depth"], s=1, alpha=0.5)
        ax.axhline(chrom_df["depth"].mean(), color="red", linewidth=0.8, linestyle="--", label="Mean")
        ax.axhline(chrom_df["depth"].median(), color="orange", linewidth=0.8, linestyle=":", label="Median")
        ax.set_title(chrom, fontsize=9)
        ax.set_xlabel("Position (Mbp)", fontsize=8)
        ax.set_ylabel("Depth", fontsize=8)
        ax.set_ylim(0, 60)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7, loc="upper right")

        # per-chromosome stats annotation (top-left corner of each subplot)
        c_mean   = chr_mean.get(chrom, float("nan"))
        c_median = chr_median["depth"].get(chrom, float("nan"))
        c_std    = chr_std["depth"].get(chrom, float("nan"))
        note = f"mean={c_mean:.1f}  median={c_median:.1f}  sd={c_std:.1f}"
        ax.text(0.02, 0.97, note, transform=ax.transAxes,
                fontsize=6.5, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    for j in range(len(present), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"{library} Coverage Depth per Chromosome")

    # global stats box anchored to the bottom of the figure
    global_note = (
        f"Whole-genome statistics (only 10 largest contigs/chromosomes)\n"
        f"Mean: {g_mean:.2f}    Median: {g_median:.2f}    Std: {g_std:.2f}\n"
        "\n"
        f"Whole-genome statistics (all contigs)\n"
        f"Mean: {ga_mean:.2f}    Median: {ga_median:.2f}    Std: {ga_std:.2f}\n"
    )
    fig.text(0.5, 0.01, global_note, ha="center", va="bottom", fontsize=8,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.8))

    fig.tight_layout(rect=[0, 0.06, 1, 1])

    base = os.path.splitext(filename)[0]
    for ext in (".png", ".svg", ".pdf"):
        out = base + ext
        print(f"Saving plot as {out}")
        plt.savefig(out, dpi=150)

if __name__ == "__main__":
    main()