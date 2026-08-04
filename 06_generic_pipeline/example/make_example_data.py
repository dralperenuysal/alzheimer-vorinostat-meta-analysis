#!/usr/bin/env python3
"""Generate the toy dataset used to smoke-test the generic pipeline.

Three small studies with a planted signal: 40 genes are up in the case group
and 40 are down, in every study, on top of gene-specific baselines and
study-specific library sizes. Everything else is noise. Running
run_meta_analysis.R over these should recover most of the 80 planted genes
and few others, which is enough to tell a working install from a broken one.

Each study is written twice from the same underlying signal: as raw counts
(for --type rnaseq) and as log2 intensities with a scatter of missing values
(for --type microarray), so both branches can be exercised. The two share one
metadata directory.

The gene symbols are real so that identifier harmonisation can be exercised;
the values are not real data and mean nothing biologically.

    python3 06_generic_pipeline/example/make_example_data.py
"""

import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20260803

STUDIES = [("studyA", 8, 8), ("studyB", 6, 7), ("studyC", 10, 9)]
N_GENES = 600
N_CHANGED = 40
CASE, CONTROL = "case", "control"

GENES = """UBB ATP5F1A PRKACA PRKACB PSMD14 UQCRC2 CDC42 NDUFS4 SDHB NDUFB5
NDUFV1 SNAP25 PDHB COX4I1 GRB2 HDAC1 CREBBP NFKB1 KAT2B SP1 IKBKB ITGB5 ITGB8
YY1 HDAC4 GSN ITGA6 YAP1 FLNB EZR APP MAPT PSEN1 PSEN2 ACTB GAPDH TUBB SDHA
CALM1 CALM2 CALM3 ATP1A1 ATP1A2 ATP1A3 ATP2B1 ATP2B2 SYN1 SYN2 SYP SYT1
NEFL NEFM NEFH GFAP AIF1 CD68 TNF IL1B IL6 CXCL8 CCL2 NFKBIA RELA TP53 MDM2
CDKN1A BAX BCL2 CASP3 CASP9 PARP1 ATM ATR BRCA1 BRCA2 RAD51 XRCC1 ERCC1
MKI67 PCNA CCND1 CCNE1 CDK2 CDK4 CDK6 RB1 E2F1 MYC MAX JUN FOS EGR1 ATF3
STAT1 STAT3 STAT5A JAK1 JAK2 SOCS1 SOCS3 IRF1""".split()


def gene_names(n):
    """Real symbols first, then padded with unique synthetic ones."""
    names = list(GENES)
    i = 1
    while len(names) < n:
        names.append("SYNTH{:04d}".format(i))
        i += 1
    return names[:n]


def main():
    rng = random.Random(SEED)
    genes = gene_names(N_GENES)

    # The same genes move in the same direction in every study; that is the
    # signal a meta-analysis is supposed to find.
    shifted = rng.sample(genes, N_CHANGED * 2)
    up_genes = set(shifted[:N_CHANGED])
    down_genes = set(shifted[N_CHANGED:])

    baseline = {gene: rng.uniform(3.0, 9.0) for gene in genes}

    expression_dir = os.path.join(HERE, "expression")
    microarray_dir = os.path.join(HERE, "expression_microarray")
    metadata_dir = os.path.join(HERE, "metadata")
    for directory in (expression_dir, microarray_dir, metadata_dir):
        os.makedirs(directory, exist_ok=True)

    for study, n_case, n_control in STUDIES:
        samples = (["{}_case{:02d}".format(study, i + 1) for i in range(n_case)] +
                   ["{}_ctrl{:02d}".format(study, i + 1) for i in range(n_control)])
        labels = [CASE] * n_case + [CONTROL] * n_control
        depth = rng.uniform(0.7, 1.4)  # study-specific library size

        rows = []
        for gene in genes:
            counts = []
            for label in labels:
                shift = 0.0
                if label == CASE:
                    if gene in up_genes:
                        shift = 1.2
                    elif gene in down_genes:
                        shift = -1.2
                mean = 2.0 ** (baseline[gene] + shift) * depth
                # Negative-binomial-ish spread: gamma-mixed Poisson.
                spread = rng.gammavariate(6.0, 1 / 6.0)
                counts.append(max(0, int(rng.gauss(mean * spread, (mean * spread) ** 0.5))))
            rows.append((gene, counts))

        path = os.path.join(expression_dir, "{}_counts.tsv".format(study))
        with open(path, "w") as handle:
            handle.write("gene\t" + "\t".join(samples) + "\n")
            for gene, counts in rows:
                handle.write(gene + "\t" + "\t".join(str(c) for c in counts) + "\n")

        # Microarray view of the same signal: log2 intensities, with ~1% of
        # values knocked out so the kNN imputation branch actually runs.
        path = os.path.join(microarray_dir, "{}_expression.tsv".format(study))
        with open(path, "w") as handle:
            handle.write("GeneSymbol\t" + "\t".join(samples) + "\n")
            for gene, counts in rows:
                values = []
                for count in counts:
                    if rng.random() < 0.01:
                        values.append("NA")
                    else:
                        values.append("{:.4f}".format(math.log2(count + 1)))
                handle.write(gene + "\t" + "\t".join(values) + "\n")

        path = os.path.join(metadata_dir, "{}_metadata.tsv".format(study))
        with open(path, "w") as handle:
            handle.write("sample\tstate\tstudy\n")
            for sample, label in zip(samples, labels):
                handle.write("{}\t{}\t{}\n".format(sample, label, study))

    truth_path = os.path.join(HERE, "planted_genes.tsv")
    with open(truth_path, "w") as handle:
        handle.write("gene\tdirection\n")
        for gene in sorted(up_genes):
            handle.write(gene + "\tup\n")
        for gene in sorted(down_genes):
            handle.write(gene + "\tdown\n")

    print("Wrote {} studies as counts ({}/) and as log2 intensities ({}/)".format(
        len(STUDIES), os.path.basename(expression_dir), os.path.basename(microarray_dir)))
    print("Planted signal recorded in {}".format(truth_path))


if __name__ == "__main__":
    main()
