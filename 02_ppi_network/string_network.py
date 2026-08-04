#!/usr/bin/env python3
"""Build a STRING protein-protein interaction network for a list of genes.

This replaces the manual step that was originally done in the STRING web
interface: paste the gene list, set the organism to Homo sapiens, restrict to
high-confidence interactions, export. The REST API returns the same network,
and pinning the API to a specific STRING release means the result does not
drift when STRING publishes a new version.

Outputs, written into --out and named after the input file:

  <name>_edges.tsv     every interaction STRING returned, all evidence
                       subscores kept, one row per pair
  <name>_network.sif   the same edges as a Cytoscape-importable SIF
  <name>_nodes.tsv     the mapped identifiers, so unmapped genes are visible
  <name>_network.svg   the STRING network image (only with --image)

Usage:
    python3 02_ppi_network/string_network.py \
        --genes results/meta_analysis/intersected_genes_up_alz.txt \
        --out results/ppi
"""

import argparse
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

# The manuscript reports STRING v12.0. The unversioned https://string-db.org
# endpoint always serves the current release, so it is pinned here instead:
# version-12-0.string-db.org will keep answering with v12.0 data after v13
# ships. Change this only if you intend to analyse a different release.
STRING_VERSION = "12.0"

# STRING asks callers to identify themselves so it can contact heavy users
# rather than silently rate-limiting them.
CALLER = "alzheimer_vorinostat_meta_analysis"

# STRING combined score, 0-1000. 700 is STRING's own "high confidence"
# threshold and is the setting reported in the manuscript.
DEFAULT_SCORE = 700

SPECIES_HUMAN = 9606


def version_host(version):
    """https://version-12-0.string-db.org for version '12.0'."""
    return "https://version-{}.string-db.org".format(version.replace(".", "-"))


def read_genes(path):
    """One gene symbol per line. A 'gene' header line is tolerated."""
    genes = []
    seen = set()
    with open(path) as handle:
        for line in handle:
            symbol = line.strip()
            if not symbol or symbol.lower() == "gene":
                continue
            if symbol not in seen:
                seen.add(symbol)
                genes.append(symbol)
    if not genes:
        sys.exit("No gene symbols found in {}".format(path))
    return genes


def string_post(endpoint, genes, score, version, extra=None, binary=False):
    """POST to a STRING API endpoint.

    POST rather than GET because a few hundred gene symbols overrun the URL
    length limit that STRING enforces on GET.
    """
    payload = {
        "identifiers": "\r".join(genes),
        "species": str(SPECIES_HUMAN),
        "required_score": str(score),
        "caller_identity": CALLER,
    }
    if extra:
        payload.update(extra)
    url = "{}/api/{}".format(version_host(version), endpoint)
    request = urllib.request.Request(url, data=urllib.parse.urlencode(payload).encode())
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            raw = response.read()
    except urllib.error.HTTPError as error:
        sys.exit("STRING returned HTTP {} for {}:\n{}".format(
            error.code, endpoint, error.read().decode("utf-8", "replace")[:500]))
    except urllib.error.URLError as error:
        sys.exit("Could not reach STRING ({}). This step needs network access.".format(error.reason))
    return raw if binary else raw.decode()


def parse_tsv(text):
    lines = [line for line in text.strip().split("\n") if line]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    return header, [line.split("\t") for line in lines[1:]]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", required=True,
                        help="file with one gene symbol per line")
    parser.add_argument("--out", default="results/ppi",
                        help="output directory (default: results/ppi)")
    parser.add_argument("--score", type=int, default=DEFAULT_SCORE,
                        help="minimum STRING combined score, 0-1000 "
                             "(default: %(default)s = high confidence)")
    parser.add_argument("--string-version", default=STRING_VERSION,
                        help="STRING release to query (default: %(default)s)")
    parser.add_argument("--image", action="store_true",
                        help="also download STRING's own network drawing as SVG")
    args = parser.parse_args()

    genes = read_genes(args.genes)
    name = os.path.splitext(os.path.basename(args.genes))[0]
    os.makedirs(args.out, exist_ok=True)

    print("{}: {} unique gene symbols".format(name, len(genes)))
    print("querying STRING v{} at score >= {}".format(args.string_version, args.score))

    # Identifier mapping first, so that genes STRING could not resolve are
    # recorded rather than silently dropped from the node count.
    header, rows = parse_tsv(string_post(
        "tsv/get_string_ids", genes, args.score, args.string_version,
        extra={"limit": "1", "echo_query": "1"}))
    nodes_path = os.path.join(args.out, name + "_nodes.tsv")
    with open(nodes_path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    mapped = len(rows)
    if mapped < len(genes):
        print("  note: {} of {} symbols mapped to a STRING protein; "
              "the rest are listed as absent from {}".format(mapped, len(genes), nodes_path))

    header, rows = parse_tsv(string_post(
        "tsv/network", genes, args.score, args.string_version))
    if not rows:
        sys.exit("STRING returned no interactions at score >= {}.".format(args.score))

    a_index = header.index("preferredName_A")
    b_index = header.index("preferredName_B")

    # STRING reports each interaction once per direction. Collapse to
    # undirected pairs so that degree is not counted twice.
    seen_pairs = set()
    unique_rows = []
    for row in rows:
        pair = tuple(sorted((row[a_index], row[b_index])))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        unique_rows.append(row)
    unique_rows.sort(key=lambda row: tuple(sorted((row[a_index], row[b_index]))))

    edges_path = os.path.join(args.out, name + "_edges.tsv")
    with open(edges_path, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in unique_rows:
            handle.write("\t".join(row) + "\n")

    sif_path = os.path.join(args.out, name + "_network.sif")
    with open(sif_path, "w") as handle:
        for row in unique_rows:
            handle.write("{}\tpp\t{}\n".format(row[a_index], row[b_index]))

    connected = {row[a_index] for row in unique_rows} | {row[b_index] for row in unique_rows}
    print("  {} interactions among {} connected proteins".format(len(unique_rows), len(connected)))
    print("  wrote {}".format(edges_path))
    print("  wrote {}  (File > Import > Network from File in Cytoscape)".format(sif_path))

    if args.image:
        svg = string_post("svg/network", genes, args.score, args.string_version, binary=True)
        image_path = os.path.join(args.out, name + "_network.svg")
        with open(image_path, "wb") as handle:
            handle.write(svg)
        print("  wrote {}".format(image_path))


if __name__ == "__main__":
    main()
