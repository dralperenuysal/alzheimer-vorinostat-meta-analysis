#!/usr/bin/env python3
"""Draw Figures 3 and 4 (PPI networks and hub genes) as vector PDFs.

Everything comes from the committed STRING output in ``results/ppi``:

  Fig3   the whole network per direction, one node per connected protein
  Fig4   the subnetwork induced on the hub genes

Degrees are read from ``*_degree.tsv`` rather than recomputed, so the sizes and
colours here are the same numbers the hub gene ranking and the manuscript use.
Recomputing them on a pruned graph would silently disagree with both.

Isolated proteins (degree 0) carry no interaction and so never appear in the
edge table; the node counts drawn are therefore the connected counts the
manuscript reports, not the query-list length.

Colour encodes degree on a single-hue sequential ramp per direction, matching
the direction colours of Figure 1, and node area encodes it a second time so
the figure survives greyscale printing and colour-vision deficiency. Figure 3
carries no gene labels — at 146 and 402 nodes any label set collides — so the
identities are named in Figure 4 instead.

Usage:  python3 04_figures/make_network_figures.py
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
PPI = os.path.join(ROOT, "results", "ppi")
OUT = os.path.join(ROOT, "results", "figures")

MM = 1 / 25.4
FULL_W = 174 * MM               # Springer double-column width

SEED = 20240917                 # fixed so the layout is reproducible

# One hue per direction, light -> dark, anchored on the Figure 1 colours.
RAMPS = {
    "up": LinearSegmentedColormap.from_list("up", ["#FDE3D0", "#D55E00", "#7A3600"]),
    "down": LinearSegmentedColormap.from_list("down", ["#D6E8F5", "#0072B2", "#003F63"]),
}
TITLES = {"up": "Upregulated", "down": "Downregulated"}


def load(direction):
    stem = os.path.join(PPI, "intersected_genes_{}_alz".format(direction))
    with open(stem + "_edges.tsv", newline="") as fh:
        edges = [(r["preferredName_A"], r["preferredName_B"], float(r["score"]))
                 for r in csv.DictReader(fh, delimiter="\t")]
    with open(stem + "_degree.tsv", newline="") as fh:
        degree = {r["gene"]: int(r["degree"]) for r in csv.DictReader(fh, delimiter="\t")}
    with open(stem + "_hub.txt") as fh:
        hubs = [l.strip() for l in fh if l.strip() and l.strip() != "gene"]
    g = nx.Graph()
    for a, b, s in edges:
        g.add_edge(a, b, score=s)
    return g, degree, hubs


def ordered_subgraph(g, nodes):
    """Induced subgraph with a canonical node and edge order.

    ``Graph.subgraph`` iterates the node set internally, so the order of the
    result depends on Python's string hash seed. spring_layout seeds its
    starting positions by node order, so that alone makes the layout differ
    between runs even with a fixed seed. Rebuilding the graph from sorted
    sequences removes the dependency.
    """
    keep = set(nodes)
    sub = nx.Graph()
    sub.add_nodes_from(sorted(keep))
    sub.add_edges_from(sorted((a, b) if a <= b else (b, a)
                              for a, b in g.edges() if a in keep and b in keep))
    return sub


def packed_layout(g, fragment_scale=0.045, k=None):
    """Force-directed layout, laid out per component and packed.

    spring_layout on a disconnected graph pushes the components apart until
    they dominate the frame. The giant component is laid out on its own and
    kept large; the remaining fragments are parked in a strip underneath it,
    where they are still visible but no longer compete for the eye.
    """
    # connected_components yields sets, whose iteration order depends on
    # Python's string hash seed and therefore changes between runs. Sorting the
    # members and then the components themselves makes the layout reproducible
    # without having to pin PYTHONHASHSEED.
    comps = sorted((sorted(c) for c in nx.connected_components(g)),
                   key=lambda c: (-len(c), c[0]))
    pos = nx.spring_layout(ordered_subgraph(g, comps[0]), seed=SEED, k=k,
                           iterations=250)
    if len(comps) == 1:
        return pos

    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    width = max(xs) - min(xs)
    scale = width * fragment_scale             # fragment radius
    gap = scale * 2.9
    x = min(xs)
    y = min(ys) - width * 0.13                 # strip sits below the component
    for comp in comps[1:]:
        sub = nx.spring_layout(ordered_subgraph(g, comp), seed=SEED, iterations=60)
        for n, (nx_, ny_) in sub.items():
            pos[n] = (x + nx_ * scale, y + ny_ * scale)
        x += gap
        if x > max(xs):                        # wrap onto another row
            x = min(xs)
            y -= gap
    return pos


def draw(ax, g, degree, hubs, direction, label_hubs, node_scale, edge_width,
         fragment_scale):
    pos = packed_layout(g, fragment_scale=fragment_scale)
    degs = [degree.get(n, 0) for n in g.nodes()]
    norm = Normalize(vmin=min(degs), vmax=max(degs))
    cmap = RAMPS[direction]

    nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#9AA0A6",
                           width=edge_width, alpha=0.55)
    nx.draw_networkx_nodes(
        g, pos, ax=ax,
        node_size=[node_scale * (0.35 + norm(d)) for d in degs],
        node_color=[cmap(norm(d)) for d in degs],
        edgecolors="white", linewidths=0.4)

    if label_hubs:
        present = [h for h in hubs if h in g]
        # In the full networks only the few strongest are labelled; at 402
        # nodes a label per hub is a pile of overlapping text.
        if label_hubs != "all":
            present = sorted(present, key=lambda n: -degree.get(n, 0))[:label_hubs]
        # Hubs sit on top of each other in the middle of a force-directed
        # layout, so labels are pushed radially outward from the centroid
        # rather than drawn on the node.
        cx = sum(p[0] for p in pos.values()) / len(pos)
        cy = sum(p[1] for p in pos.values()) / len(pos)
        span = max(max(p[0] for p in pos.values()) - min(p[0] for p in pos.values()),
                   1e-9)
        for n in present:
            x, y = pos[n]
            dx, dy = x - cx, y - cy
            norm_d = (dx * dx + dy * dy) ** 0.5 or 1.0
            off = span * 0.055
            ax.annotate(
                n, xy=(x, y), xytext=(x + dx / norm_d * off, y + dy / norm_d * off),
                fontsize=4.6, ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="-", lw=0.3, color="#5F6368",
                                shrinkA=0, shrinkB=1))

    ax.set_axis_off()
    ax.margins(0.06)
    return norm, cmap


def figure(kind, path, node_scale, edge_width, label_hubs, fragment_scale):
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, FULL_W * 0.46))
    for ax, direction, tag in zip(axes, ("up", "down"), ("(a)", "(b)")):
        g, degree, hubs = load(direction)
        if kind == "hub":
            g = ordered_subgraph(g, [h for h in hubs if h in g])
        norm, cmap = draw(ax, g, degree, hubs, direction,
                          label_hubs=label_hubs, node_scale=node_scale,
                          edge_width=edge_width, fragment_scale=fragment_scale)
        ax.set_title("{}  {} — {} proteins, {} interactions".format(
            tag, TITLES[direction], g.number_of_nodes(), g.number_of_edges()),
            fontsize=6.5, loc="left")
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                          fraction=0.028, pad=0.01)
        cb.set_label("degree", fontsize=5.5)
        cb.ax.tick_params(labelsize=5, length=2)
        cb.outline.set_visible(False)
    fig.tight_layout(pad=0.4)
    fig.savefig(path, format="pdf")
    plt.close(fig)
    print("wrote {}".format(os.path.relpath(path, ROOT)))


def main():
    os.makedirs(OUT, exist_ok=True)
    # Fig3 carries no gene labels: at 146 and 402 nodes the hubs sit on top of
    # one another and any label set collides. Figure 4 names them instead.
    figure("full", os.path.join(OUT, "Fig3.pdf"),
           node_scale=30, edge_width=0.22, label_hubs=0, fragment_scale=0.045)
    figure("hub", os.path.join(OUT, "Fig4.pdf"),
           node_scale=340, edge_width=0.7, label_hubs="all", fragment_scale=0.30)


if __name__ == "__main__":
    main()
