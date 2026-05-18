#!/home/daffa/Work/2026/thesis/.venv/bin/python

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SAMPLES = {
    "SBC4": "resources/depth/r0074_depth.pkl",
    "SBC10": "resources/depth/r0066_depth.pkl",
    "SBC11": "resources/depth/SBC11_depth.pkl",
    "SBC23": "resources/depth/r0076_depth.pkl",
}

COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

VALID_CONTIGS = [f"NC_01287{n}.2" for n in range(0, 10)]

BIN_SIZE = 1000_000


def load_sample(name, path):
    print(f"Loading {name} from {path}")
    df = pd.read_pickle(path)
    return df


def calc_stats(df):
    return df["depth"].mean(), df["depth"].median(), df["depth"].std()


def bin_chrom(chrom_df, bin_size):
    chrom_df = chrom_df.copy()
    chrom_df["bin"] = (chrom_df["pos"] // bin_size) * bin_size
    return chrom_df.groupby("bin")["depth"].mean().reset_index()


def main():
    samples_data = {}
    for name, path in SAMPLES.items():
        df = load_sample(name, path)
        df_valid = df[df["chrom"].isin(VALID_CONTIGS)].copy()
        samples_data[name] = {"all": df, "valid": df_valid}

    present = [c for c in VALID_CONTIGS if any(
        c in d["valid"]["chrom"].values for d in samples_data.values()
    )]

    ncols = 2
    nrows = (len(present) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, nrows * 3.2 + 2.5), sharex=False, sharey=True)
    axes = axes.flatten()

    for i, chrom in enumerate(present):
        ax = axes[i]

        for (name, color) in zip(SAMPLES.keys(), COLORS):
            df_valid = samples_data[name]["valid"]
            chrom_df = df_valid[df_valid["chrom"] == chrom]
            if chrom_df.empty:
                continue

            binned = bin_chrom(chrom_df, BIN_SIZE)
            ax.plot(binned["bin"] / 1_000_000, binned["depth"],
                    linewidth=0.9, alpha=0.75, color=color, label=name)
            ax.axhline(chrom_df["depth"].mean(),
                       color=color, linewidth=0.8, linestyle="--", alpha=0.8)

        ax.set_title(chrom, fontsize=9)
        ax.set_xlabel("Position (Mbp)", fontsize=8)
        ax.set_ylabel("Depth", fontsize=8)
        ax.set_ylim(0, 80)
        ax.tick_params(labelsize=7)

        # per-chromosome stats annotation (top-left)
        stat_lines = []
        for name in SAMPLES:
            df_valid = samples_data[name]["valid"]
            chrom_df = df_valid[df_valid["chrom"] == chrom]
            if chrom_df.empty:
                continue
            m, med, s = calc_stats(chrom_df)
            stat_lines.append(f"{name}: mean={m:.1f}  med={med:.1f}  sd={s:.1f}")
        note = "\n".join(stat_lines)
        ax.text(0.02, 0.98, note, transform=ax.transAxes,
                fontsize=5.5, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.75))

    # hide unused axes
    for j in range(len(present), len(axes)):
        axes[j].set_visible(False)

    # shared legend
    handles = [mpatches.Patch(color=c, label=n) for n, c in zip(SAMPLES.keys(), COLORS)]
    fig.legend(handles=handles, loc="upper right", fontsize=8, title="Sample",
               title_fontsize=8, framealpha=0.9)

    fig.suptitle("Coverage Depth per Chromosome — All Samples", fontsize=12, y=1.005)

    # global stats table at bottom
    rows = []
    for name in SAMPLES:
        gm, gmed, gsd = calc_stats(samples_data[name]["valid"])
        am, amed, asd = calc_stats(samples_data[name]["all"])
        rows.append(
            f"{name}: 10-chr mean={gm:.1f} med={gmed:.1f} sd={gsd:.1f} | "
            f"all-contig mean={am:.1f} med={amed:.1f} sd={asd:.1f}"
        )
    global_note = "Whole-genome statistics\n" + "\n".join(rows)
    fig.text(0.5, 0.005, global_note, ha="center", va="bottom", fontsize=7,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.85))

    fig.tight_layout(rect=[0, 0.10, 1, 1])

    output_base = "analysis/00_data_quality/all_samples"
    for ext in (".png", ".svg", ".pdf"):
        out = output_base + ext
        print(f"Saving {out}")
        plt.savefig(out, dpi=150, bbox_inches="tight")

    print("Done.")


if __name__ == "__main__":
    main()