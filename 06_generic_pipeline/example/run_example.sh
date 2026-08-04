#!/usr/bin/env bash
#
# Smoke test for the generic pipeline. Runs both assay branches over the toy
# data and checks that the planted signal comes back out.
#
# Run from the repository root:
#     bash 06_generic_pipeline/example/run_example.sh
#
# Takes about a minute and needs no network access. If this passes, your R
# environment is set up correctly and you can point the pipeline at real data.

set -euo pipefail

EXAMPLE="06_generic_pipeline/example"
OUT="${EXAMPLE}/output"

if [ ! -f "06_generic_pipeline/run_meta_analysis.R" ]; then
    echo "Run this from the repository root." >&2
    exit 1
fi

echo "== regenerating toy data =="
python3 "${EXAMPLE}/make_example_data.py"

echo
echo "== RNA-seq branch =="
Rscript 06_generic_pipeline/run_meta_analysis.R \
    --expression "${EXAMPLE}/expression" \
    --metadata "${EXAMPLE}/metadata" \
    --type rnaseq \
    --case case --control control \
    --no-id-harmonisation \
    --out "${OUT}/rnaseq"

echo
echo "== microarray branch =="
Rscript 06_generic_pipeline/run_meta_analysis.R \
    --expression "${EXAMPLE}/expression_microarray" \
    --metadata "${EXAMPLE}/metadata" \
    --type microarray \
    --case case --control control \
    --no-id-harmonisation \
    --out "${OUT}/microarray"

echo
echo "== recovery of the planted signal =="
python3 - "$OUT" "$EXAMPLE" <<'PYTHON'
import csv
import sys

out_dir, example_dir = sys.argv[1], sys.argv[2]

truth = {}
with open(example_dir + "/planted_genes.tsv") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        truth[row["gene"]] = row["direction"]

failures = 0
for branch in ("rnaseq", "microarray"):
    for direction in ("up", "down"):
        with open("{}/{}/meta_genes_{}.txt".format(out_dir, branch, direction)) as handle:
            called = [line.strip() for line in handle if line.strip()]
        planted = [g for g, d in truth.items() if d == direction]
        correct = [g for g in called if truth.get(g) == direction]
        wrong = sorted(set(called) - set(correct))
        recall = len(correct) / len(planted)
        precision = len(correct) / len(called) if called else 0.0
        print("  {:<10} {:<4}  recall {:>3}/{:<3} ({:.0%})   precision {:.0%}   {} false positive(s)".format(
            branch, direction, len(correct), len(planted), recall, precision, len(wrong)))
        # Loose thresholds on purpose: this checks that the pipeline runs and
        # finds signal, not that it hits some particular number.
        if recall < 0.5 or precision < 0.8:
            print("    UNEXPECTED: {} {} recovered poorly. False positives: {}".format(
                branch, direction, ", ".join(wrong) or "none"))
            failures += 1

if failures:
    print("\nFAILED: {} of 4 checks looked wrong.".format(failures))
    sys.exit(1)
print("\nOK: both branches ran and recovered the planted signal.")
PYTHON

echo
echo "Results are under ${OUT}/. Delete that directory to clean up."
