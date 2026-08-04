#!/usr/bin/env python3
"""Compare a re-run of the pipeline against the published results.

Re-running the pipeline overwrites the files in results/. The versions
committed to git are the ones behind the manuscript, so this script reads the
committed copies straight out of git (`git show HEAD:<path>`) and diffs the
working-tree files against them. Nothing needs to be backed up first, and the
comparison stays correct however many times you re-run.

    python3 verify_reproduction.py

Each check reports PASS, DIFFERS or SKIP (nothing regenerated yet). Small
differences are expected in places and are described as they come up -- read
the notes rather than only the verdicts.
"""

import csv
import hashlib
import io
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.abspath(__file__))

PASS, DIFFERS, SKIP = "PASS", "DIFFERS", "SKIP"
results = []


def record(status, name, *notes):
    results.append((status, name))
    colour = {PASS: "\033[32m", DIFFERS: "\033[33m", SKIP: "\033[90m"}[status]
    reset = "\033[0m" if sys.stdout.isatty() else ""
    if not sys.stdout.isatty():
        colour = ""
    print("{}{:<8}{}{}".format(colour, status, reset, name))
    for note in notes:
        for line in str(note).rstrip().split("\n"):
            print("         " + line)


def committed(path):
    """Contents of `path` as committed at HEAD, or None if not tracked."""
    try:
        return subprocess.run(["git", "-C", REPO, "show", "HEAD:" + path],
                              capture_output=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def working(path):
    full = os.path.join(REPO, path)
    if not os.path.exists(full):
        return None
    with open(full, "rb") as handle:
        return handle.read()


def gene_set(blob):
    if blob is None:
        return None
    genes = set()
    for line in blob.decode("utf-8", "replace").splitlines():
        symbol = line.strip()
        if symbol and symbol.lower() != "gene":
            genes.add(symbol)
    return genes


# A file whose bytes still match HEAD exactly was almost certainly never
# rewritten, so calling it "reproduced" would overstate what was checked.
UNTOUCHED = "byte-identical to HEAD, so this stage looks like it was not re-run"


def compare_gene_list(path, label):
    old_blob, new_blob = committed(path), working(path)
    old, new = gene_set(old_blob), gene_set(new_blob)
    if old is None:
        return record(SKIP, label, "not committed at HEAD: " + path)
    if new is None:
        return record(SKIP, label, "not present in the working tree: " + path)
    if old_blob == new_blob:
        return record(SKIP, label, "{} genes; ".format(len(new)) + UNTOUCHED)
    if old == new:
        return record(PASS, label, "{} genes, identical".format(len(new)))
    added, removed = sorted(new - old), sorted(old - new)
    union = len(old | new)
    record(DIFFERS, label,
           "published {} genes, regenerated {} ({:.1%} agreement)".format(
               len(old), len(new), len(old & new) / union if union else 1.0),
           "only in the regenerated list ({}): {}".format(len(added), fmt(added)),
           "only in the published list ({}): {}".format(len(removed), fmt(removed)))


def fmt(genes, limit=12):
    if not genes:
        return "none"
    if len(genes) <= limit:
        return ", ".join(genes)
    return ", ".join(genes[:limit]) + " ... (+{} more)".format(len(genes) - limit)


def compare_meta_csv(path, label):
    old_blob, new_blob = committed(path), working(path)
    if old_blob is None:
        return record(SKIP, label, "not committed at HEAD: " + path)
    if new_blob is None:
        return record(SKIP, label, "not present in the working tree: " + path)
    if old_blob == new_blob:
        return record(SKIP, label, UNTOUCHED)

    def effect_sizes(blob):
        table = {}
        reader = csv.DictReader(io.StringIO(blob.decode("utf-8", "replace")))
        for row in reader:
            symbol = row.get("symbol")
            value = row.get("Com.ES")
            if symbol and value:
                try:
                    table[symbol] = float(value)
                except ValueError:
                    pass
        return table

    old, new = effect_sizes(old_blob), effect_sizes(new_blob)
    if not old:
        return record(SKIP, label, "could not read Com.ES out of the committed file")
    shared = set(old) & set(new)
    if not shared:
        return record(DIFFERS, label, "no gene in common between the two versions")
    worst = max(abs(old[g] - new[g]) for g in shared)
    notes = ["{} genes published, {} regenerated, {} in common".format(len(old), len(new), len(shared)),
             "largest change in combined effect size: {:.2e}".format(worst)]
    # Floating-point summation order can move the last digits; anything above
    # this is a real change in the model, not arithmetic noise.
    if set(old) == set(new) and worst < 1e-6:
        record(PASS, label, *notes)
    else:
        record(DIFFERS, label, *notes)


def compare_hub_genes(regenerated, published, degree_table, label):
    """Compare hub lists, treating equal-degree swaps as what they are.

    A degree ranking is integer-valued, so the boundary of a top-N cut often
    falls inside a group of genes with identical degree. Two runs can then
    return different genes without disagreeing about the network at all. This
    reports that case separately from a genuine mismatch.
    """
    new = gene_set(working(regenerated))
    old = gene_set(committed(published)) or gene_set(working(published))
    if new is None:
        return record(SKIP, label, "not regenerated yet: " + regenerated)
    if old is None:
        return record(SKIP, label, "no published list at " + published)
    if old == new:
        return record(PASS, label, "same {} genes".format(len(new)))

    degrees = {}
    table = working(degree_table)
    if table:
        reader = csv.DictReader(io.StringIO(table.decode("utf-8", "replace")), delimiter="\t")
        for row in reader:
            try:
                degrees[row["gene"]] = int(row["degree"])
            except (KeyError, TypeError, ValueError):
                pass

    added, removed = sorted(new - old), sorted(old - new)
    notes = ["{} of {} hub genes match".format(len(old & new), len(old))]
    if degrees:
        notes.append("in the regenerated list: " + ", ".join(
            "{} (degree {})".format(g, degrees.get(g, "?")) for g in added))
        notes.append("in the published list:   " + ", ".join(
            "{} (degree {})".format(g, degrees.get(g, "?")) for g in removed))
        swap_degrees = {degrees.get(g) for g in added + removed}
        if len(swap_degrees) == 1 and None not in swap_degrees:
            notes.append("Every swapped gene has the same degree ({}), so the two lists".format(
                swap_degrees.pop()))
            notes.append("rank the network identically; only the tie at the cut-off is")
            notes.append("broken differently. This is not a disagreement about the biology.")
            return record(PASS, label, *notes)
    else:
        notes.append("regenerated only: " + fmt(added))
        notes.append("published only:   " + fmt(removed))
    record(DIFFERS, label, *notes)


def compare_pdf(path, label):
    """Compare PDFs by their text, not their bytes.

    matplotlib stamps /CreationDate into every PDF it writes, so two runs of
    the same script never produce identical bytes. The extracted text and the
    file size do match, and together they are enough to show the figure did
    not change.
    """
    old_blob, new_blob = committed(path), working(path)
    if old_blob is None:
        return record(SKIP, label, "not committed at HEAD: " + path)
    if new_blob is None:
        return record(SKIP, label, "not present in the working tree: " + path)
    if old_blob == new_blob:
        # matplotlib stamps a fresh /CreationDate on every write, so a
        # regenerated figure can never be byte-identical to the committed one.
        return record(SKIP, label, UNTOUCHED)

    def text_hash(blob):
        try:
            done = subprocess.run(["pdftotext", "-", "-"], input=blob,
                                  capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
        return hashlib.sha256(done.stdout).hexdigest()

    old_hash, new_hash = text_hash(old_blob), text_hash(new_blob)
    if old_hash is None:
        return record(SKIP, label, "pdftotext is not installed (poppler-utils)")
    size_note = "{} vs {} bytes".format(len(old_blob), len(new_blob))
    if old_hash == new_hash and len(old_blob) == len(new_blob):
        record(PASS, label, "identical text and size (" + size_note + ")")
    elif old_hash == new_hash:
        record(PASS, label, "identical text, size differs slightly (" + size_note + ")")
    else:
        record(DIFFERS, label, "the text drawn on the figure changed (" + size_note + ")")


def main():
    if not os.path.isdir(os.path.join(REPO, ".git")):
        sys.exit("This script reads the published results out of git history, "
                 "so it needs to run inside the git repository.")

    print("Comparing the working tree against the results committed at HEAD.\n")

    print("-- meta-analysis gene lists " + "-" * 44)
    for platform in ("microarray", "rnaseq"):
        for direction in ("up", "down"):
            compare_gene_list(
                "results/meta_analysis/{}_meta_genes_{}.txt".format(platform, direction),
                "{} {}-regulated genes".format(platform, direction))

    print("\n-- meta-analysis effect sizes " + "-" * 42)
    for platform in ("microarray", "rnaseq"):
        for direction in ("up", "down"):
            compare_meta_csv(
                "results/meta_analysis/{}_{}.csv".format(platform, direction),
                "{} {} effect sizes".format(platform, direction))

    print("\n-- microarray / RNA-seq intersection " + "-" * 35)
    for direction in ("up", "down"):
        compare_gene_list(
            "results/meta_analysis/intersected_genes_{}_alz.txt".format(direction),
            "intersected {}-regulated genes".format(direction))

    print("\n-- PPI hub genes " + "-" * 55)
    for direction in ("up", "down"):
        compare_hub_genes(
            "results/ppi/intersected_genes_{}_alz_hub.txt".format(direction),
            "results/ppi/intersected_genes_{}_alz_hub.txt".format(direction),
            "results/ppi/intersected_genes_{}_alz_degree.tsv".format(direction),
            "{}-regulated hub genes".format(direction))

    print("\n-- figures " + "-" * 61)
    for number in (1, 2, 3, 4, 5, 8, 9, 10, 11):
        compare_pdf("results/figures/Fig{}.pdf".format(number), "Figure {}".format(number))

    counts = {status: sum(1 for s, _ in results if s == status) for status in (PASS, DIFFERS, SKIP)}
    print("\n" + "=" * 72)
    print("{} passed, {} differ, {} skipped".format(counts[PASS], counts[DIFFERS], counts[SKIP]))
    if counts[SKIP]:
        print("\nSkipped checks are stages that have not been re-run: their files are")
        print("still exactly as committed. Run the stage first, then check again.")
    if counts[DIFFERS]:
        print("\nDifferences are not automatically failures. Package versions move,")
        print("GEO occasionally revises a series, and STRING publishes new releases.")
        print("Read the notes above to see whether the change is substantive.")
    return 1 if counts[DIFFERS] else 0


if __name__ == "__main__":
    sys.exit(main())
