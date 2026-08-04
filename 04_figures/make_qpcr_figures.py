#!/usr/bin/env python3
"""Draw Figures 8-11 (qPCR relative expression) as vector PDFs.

Replaces the original raster panels, which sat at 86-120 dpi at print size and
could not meet the journal's 300/600 dpi requirement.  Every value here comes
from ``supplementary_qpcr_data.csv`` through the same mean-dCt pipeline that
produces the log2FC tables, so the figures and the tables cannot disagree.

Each figure is one chart: log2 fold change per gene on a zero baseline, with the
two per-replicate values drawn on top of the bar.  Replicates are paired by well
group, which cancels the systematic offset between the two groups of the plate,
and their mean is exactly the tabulated log2FC.  No significance markers are
drawn, since no inferential statistics were performed.

Usage:  python3 04_figures/make_qpcr_figures.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from qpcr_common import (  # noqa: E402
    HUBS, MARKERS, MATURE, MODEL, VORINO,
    load, log2fc, paired_log2fc,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
OUT = os.path.join(ROOT, "results", "figures")

MM = 1 / 25.4
FULL_W = 174 * MM          # Springer double-column width
HALF_W = 84 * MM           # Springer single-column width

UP_FACE = "#b6b6b6"
DOWN_FACE = "#7d7d7d"

plt.rcParams.update({
    "font.family": "sans-serif",
    # Liberation Sans is metrically identical to Arial, which the journal asks for.
    "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "pdf.fonttype": 42,      # embed as TrueType, not Type 3
    "ps.fonttype": 42,
})

def figure(genes, test, calibrator, filename, width, height_mm, legend_cols):
    data = load()
    fig, ax = plt.subplots(figsize=(width, height_mm * MM))

    values = [log2fc(data, g, test, calibrator) for g in genes]
    x = list(range(len(genes)))

    ax.bar(x, values, width=0.6, linewidth=0.9, edgecolor="black", zorder=2,
           color=[UP_FACE if v >= 0 else DOWN_FACE for v in values])

    for i, gene in enumerate(genes):
        points = paired_log2fc(data, gene, test, calibrator)
        offsets = [-0.13, 0.13][: len(points)] or []
        for dx, (value, censored) in zip(offsets, points):
            ax.plot(i + dx, value, marker="o", markersize=3.6,
                    markerfacecolor="white" if censored else "black",
                    markeredgecolor="black", markeredgewidth=0.8,
                    linestyle="none", zorder=4, clip_on=False)

    ax.axhline(0, color="black", linewidth=1.0, zorder=3)
    ax.set_xticks(x)
    # style="italic" rather than mathtext: mathtext leaves the trailing digits upright.
    ax.set_xticklabels(genes, style="italic")
    ax.set_xlim(-0.7, len(genes) - 0.3)
    # The comparison itself is named in the caption, not on the axis.
    ax.set_ylabel("log$_2$ fold change", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", length=0, pad=6, labelsize=8)
    ax.tick_params(axis="y", length=3, labelsize=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    ax.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor="black", markeredgecolor="black",
               label="replicate"),
        Line2D([], [], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor="white", markeredgecolor="black",
               label="replicate involving a Cq 35 ceiling well"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=legend_cols,
       frameon=False, fontsize=7, handletextpad=0.4, columnspacing=1.6)

    fig.tight_layout()
    path = os.path.join(OUT, filename)
    fig.savefig(path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(path)}  ({len(genes)} genes)")


def main():
    figure(MARKERS, MODEL, MATURE, "Fig8.pdf", HALF_W, 62, 1)
    figure(HUBS, MODEL, MATURE, "Fig9.pdf", FULL_W, 70, 2)
    figure(MARKERS, VORINO, MODEL, "Fig10.pdf", HALF_W, 62, 1)
    figure(HUBS, VORINO, MODEL, "Fig11.pdf", FULL_W, 70, 2)


if __name__ == "__main__":
    main()
