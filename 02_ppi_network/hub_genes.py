#!/usr/bin/env python3
"""Rank the proteins in a STRING network by degree centrality.

This replaces the manual Cytoscape + CytoHubba step. CytoHubba's "Degree"
method is the number of edges incident on a node in the network it is given,
so it is reproduced exactly by counting edges here -- no layout, no plugin,
no GUI session to save.

The one thing a GUI cannot tell you is where the ranking is arbitrary. Degree
is an integer over a few hundred proteins, so ties are common, and a tie that
straddles the top-N boundary decides which genes get carried forward on
nothing more than sort order. Every output here therefore carries the tie
size, and --top warns when the cut falls inside a tied group.

Outputs, named after the input edge file:

  <name>_degree.tsv    every connected protein, ranked, with its tie group
  <name>_hub.txt       the top N gene symbols

Usage:
    python3 02_ppi_network/hub_genes.py \
        --edges results/ppi/intersected_genes_up_alz_edges.tsv \
        --top 14 --out results/ppi
"""

import argparse
import collections
import os
import sys


def read_edges(path):
    """Read a STRING edge table, or a two/three-column SIF."""
    pairs = set()
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        sys.exit("No edges in {}".format(path))

    header = lines[0].split("\t")
    if "preferredName_A" in header and "preferredName_B" in header:
        a_index = header.index("preferredName_A")
        b_index = header.index("preferredName_B")
        body = lines[1:]
    else:
        # SIF: source <tab> interaction <tab> target, no header
        a_index, b_index = 0, 2
        body = lines

    for line in body:
        fields = line.split("\t")
        if len(fields) <= max(a_index, b_index):
            continue
        a, b = fields[a_index], fields[b_index]
        if a and b and a != b:
            pairs.add(tuple(sorted((a, b))))
    if not pairs:
        sys.exit("Could not parse any interactions out of {}".format(path))
    return pairs


def rank_by_degree(pairs):
    """Return [(gene, degree, rank, tie_size)], highest degree first.

    Rank is competition-style: every member of a tied group gets the same
    rank, and the next distinct degree skips ahead. Within a tie, genes are
    ordered alphabetically purely so that repeated runs agree with each other
    -- alphabetical order carries no biological meaning.
    """
    degree = collections.Counter()
    for a, b in pairs:
        degree[a] += 1
        degree[b] += 1

    tie_size = collections.Counter(degree.values())
    ordered = sorted(degree.items(), key=lambda item: (-item[1], item[0]))

    ranked = []
    rank = 0
    previous_degree = None
    for position, (gene, value) in enumerate(ordered, start=1):
        if value != previous_degree:
            rank = position
            previous_degree = value
        ranked.append((gene, value, rank, tie_size[value]))
    return ranked


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--edges", required=True,
                        help="STRING edge table or SIF from string_network.py")
    parser.add_argument("--top", type=int, default=14,
                        help="how many hub genes to write out (default: %(default)s)")
    parser.add_argument("--out", default="results/ppi",
                        help="output directory (default: results/ppi)")
    args = parser.parse_args()

    pairs = read_edges(args.edges)
    ranked = rank_by_degree(pairs)

    name = os.path.basename(args.edges)
    for suffix in ("_edges.tsv", "_network.sif", ".tsv", ".sif"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    os.makedirs(args.out, exist_ok=True)

    degree_path = os.path.join(args.out, name + "_degree.tsv")
    with open(degree_path, "w") as handle:
        handle.write("gene\tdegree\trank\ttied_with\n")
        for gene, value, rank, tie in ranked:
            handle.write("{}\t{}\t{}\t{}\n".format(gene, value, rank, tie))

    top = ranked[:args.top]
    hub_path = os.path.join(args.out, name + "_hub.txt")
    with open(hub_path, "w") as handle:
        handle.write("gene\n")
        for gene, _, _, _ in top:
            handle.write(gene + "\n")

    print("{}: {} interactions among {} connected proteins".format(
        name, len(pairs), len(ranked)))
    print("  wrote {}".format(degree_path))
    print("  wrote {}".format(hub_path))
    print()
    print("  rank  degree  gene")
    for gene, value, rank, tie in top:
        marker = "   (tied with {} others)".format(tie - 1) if tie > 1 else ""
        print("  {:>4}  {:>6}  {}{}".format(rank, value, gene, marker))

    # Warn when the cut-off falls inside a tied group: the genes at the
    # boundary are interchangeable on this metric, and which of them lands
    # inside the top N is decided by sort order alone.
    if len(ranked) > args.top:
        boundary_degree = ranked[args.top - 1][1]
        if ranked[args.top][1] == boundary_degree:
            contenders = sorted(g for g, v, _, _ in ranked if v == boundary_degree)
            inside = sorted(g for g, _, _, _ in top if _degree_of(ranked, g) == boundary_degree)
            print()
            print("  WARNING: the top-{} cut falls inside a tie at degree {}.".format(
                args.top, boundary_degree))
            print("  {} genes share that degree: {}".format(len(contenders), ", ".join(contenders)))
            print("  {} of them fit inside the cut ({}); the rest are excluded on".format(
                len(inside), ", ".join(inside)))
            print("  sort order alone, not on any difference in centrality.")


def _degree_of(ranked, gene):
    for name, value, _, _ in ranked:
        if name == gene:
            return value
    return None


if __name__ == "__main__":
    main()
