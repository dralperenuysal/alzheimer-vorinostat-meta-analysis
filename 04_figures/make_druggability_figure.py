#!/usr/bin/env python3
"""Draw Figure 5 (target selection) as a vector PDF.

The selection argument has two legs and the figure shows both.

Panel (a) puts network centrality against chemical tractability, one point per
hub gene. Being central is not enough — the mitochondrial genes sit high on
degree and on the floor for ligands — and having ligands is not enough either.
HDAC1 is the only gene high on both axes, which is the selection rule made
visible rather than asserted.

Panel (b) then asks why vorinostat among HDAC1's ligands. It is not the most
potent compound in the set, and the figure says so; it was chosen because it is
clinically approved with characterised CNS distribution, and the distribution
shows it still sits well inside the potent tail.

Ligand counts measure how much attention a protein has had as a drug target as
much as anything intrinsic. For a repurposing study that is the point — an
established chemical toolkit is what makes a target actionable — but the
caption says what is being counted.

Usage:  python3 04_figures/make_druggability_figure.py
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
PPI = os.path.join(ROOT, "results", "ppi")
DRUG = os.path.join(ROOT, "results", "druggability")
OUT = os.path.join(ROOT, "results", "figures")

MM = 1 / 25.4
FULL_W = 174 * MM

UP, DOWN = "#D55E00", "#0072B2"       # the direction colours of Figure 1
SELECTED = "CHEMBL98"                  # vorinostat
TARGET = "HDAC1"


def degrees():
    """Degree per gene, scaled by the top degree of its own network.

    The two networks are different sizes (146 proteins and 185 interactions
    against 402 and 1,109), so raw degrees are not comparable between them:
    the down-regulated hubs reach 30 where the up-regulated top is 17, which
    reflects network density rather than greater centrality. Expressing each
    degree as a fraction of the highest degree in its own network puts both on
    one axis honestly, and places the most central protein of each network at
    1.0.
    """
    out = {}
    for direction in ("up", "down"):
        path = os.path.join(PPI, "intersected_genes_{}_alz_degree.tsv".format(direction))
        with open(path, newline="") as fh:
            rows = [(r["gene"], int(r["degree"])) for r in csv.DictReader(fh, delimiter="\t")]
        top = max(d for _, d in rows)
        for gene, d in rows:
            out[gene] = (d / top, direction)
    return out


def hub_genes():
    out = []
    for direction in ("up", "down"):
        path = os.path.join(PPI, "intersected_genes_{}_alz_hub.txt".format(direction))
        with open(path) as fh:
            out += [l.strip() for l in fh if l.strip() and l.strip() != "gene"]
    return out


def ligands(genes):
    counts, potency = {g: set() for g in genes}, []
    with open(os.path.join(DRUG, "hub_gene_ligands.tsv"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = r["gene"]
            counts.setdefault(gene, set())
            if r["molecule_chembl_id"]:
                counts[gene].add(r["molecule_chembl_id"])
            if gene == TARGET and r["pchembl_value"]:
                potency.append((r["molecule_chembl_id"], float(r["pchembl_value"])))
    return {g: len(m) for g, m in counts.items()}, potency


def main():
    os.makedirs(OUT, exist_ok=True)
    deg = degrees()
    genes = [g for g in hub_genes() if g in deg]
    counts, potency = ligands(genes)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(FULL_W, FULL_W * 0.40))

    # ---- (a) centrality against tractability
    for direction, colour, label in ((("up"), UP, "Upregulated"),
                                     (("down"), DOWN, "Downregulated")):
        xs = [deg[g][0] for g in genes if deg[g][1] == direction]
        ys = [counts[g] for g in genes if deg[g][1] == direction]
        ax.scatter(xs, ys, s=14, c=colour, alpha=0.85, edgecolors="white",
                   linewidths=0.4, label=label, zorder=3)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(-0.35, 6e4)
    ax.set_xlim(0.05, 1.22)
    ax.set_xlabel("Degree centrality (fraction of its network's highest)", fontsize=7)
    ax.set_ylabel("Distinct ligands with an IC$_{50}$", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=6, frameon=False, loc="upper left", handletextpad=0.3)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Only the genes that carry the argument are named. The floor is a row of
    # zeros whose identities are not the point; the count below says how many.
    named = {"HDAC1": (-8, 8), "CREBBP": (-6, 7), "IKBKB": (5, 2),
             "FYN": (4, 4), "SIRT1": (5, -7), "PSMD14": (5, -1),
             "KAT2B": (4, -8), "ATP5F1A": (-10, 7)}
    for g, off in named.items():
        if g not in counts or g not in deg:
            continue
        ax.annotate(g, (deg[g][0], counts[g]), fontsize=5.6,
                    fontweight="bold" if g == TARGET else "normal",
                    textcoords="offset points", xytext=off, zorder=4)
    n_zero = sum(1 for g in genes if counts[g] == 0)
    ax.text(0.03, 0.15, "{} of {} genes have no annotated ligand".format(n_zero, len(genes)),
            transform=ax.transAxes, fontsize=5.6, color="#5F6368")

    # ---- (b) where the selected compound sits among that target's ligands
    vals = [v for _, v in potency]
    bx.hist(vals, bins=40, color=UP, alpha=0.75, edgecolor="white", linewidth=0.2)
    sel = next((v for m, v in potency if m == SELECTED), None)
    if sel is not None:
        pct = 100.0 * sum(1 for v in vals if v < sel) / len(vals)
        bx.axvline(sel, color="#202124", linewidth=0.9, linestyle="--", zorder=4)
        bx.annotate("vorinostat\npChEMBL {:.2f} ({:.0f}th pct)".format(sel, pct),
                    xy=(sel, bx.get_ylim()[1] * 0.80), fontsize=5.8, ha="left",
                    xytext=(6, 0), textcoords="offset points")
    bx.set_xlabel("pChEMBL (median per compound)", fontsize=7)
    bx.set_ylabel("Compounds", fontsize=7)
    bx.set_title("{} ligands (n = {:,})".format(TARGET, len(vals)),
                 fontsize=7, loc="left")
    bx.tick_params(labelsize=6)
    for spine in ("top", "right"):
        bx.spines[spine].set_visible(False)

    for axis, tag in ((ax, "(a)"), (bx, "(b)")):
        axis.text(-0.10, 1.06, tag, transform=axis.transAxes,
                  fontsize=8, fontweight="bold", va="top")

    fig.tight_layout(pad=0.5)
    fig.savefig(os.path.join(OUT, "Fig5.pdf"), format="pdf")
    plt.close(fig)
    print("wrote results/figures/Fig5.pdf  "
          "({} genes, {} {} ligands)".format(len(genes), len(vals), TARGET))


if __name__ == "__main__":
    main()
