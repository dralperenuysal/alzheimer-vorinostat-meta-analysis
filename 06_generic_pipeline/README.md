# Meta-analysis pipeline for your own datasets

The rest of this repository reproduces one specific study. This directory is
the same machinery with the Alzheimer's specifics removed, so you can point it
at your own expression matrices.

You provide a directory of expression matrices and a directory of sample
metadata. It pairs them by file name, normalises each dataset according to its
assay type, runs the DExMA random-effects model, and writes ranked gene lists
that feed straight into the PPI step in `02_ppi_network/`.

## Where this pipeline starts, and what you have to do first

**It assumes your data is already prepared.** It downloads nothing and curates
nothing. You bring a gene × sample matrix per dataset and a metadata table per
dataset in which one column already says, for every sample, which group it
belongs to.

That assumption is deliberate, not a shortcut. Getting from a raw accession to
a matrix this pipeline can read is the part that cannot be automated: which
metadata column holds the diagnosis, and how a case sample is spelled, differ
in every series, and probe identifiers differ per platform, so each dataset
needs its own annotation table and its own probe-to-gene collapse. Somebody has
to open each one, read what the submitters wrote, and decide.

So the work before this pipeline is:

1. Get the expression data, whatever the source.
2. Map probes or transcript identifiers to gene symbols, and collapse
   duplicates. Use one identifier type across every dataset.
3. Build a metadata table per dataset with a sample ID column and a group
   column holding two consistent labels.
4. Check the sample IDs in the metadata match the column names in the matrix.

Step 4 is where most first runs fail, which is what `--dry-run` is for: it
reports the sample overlap per dataset without starting the analysis.

## Quick check that it works

```bash
bash 06_generic_pipeline/example/run_example.sh
```

This builds three small synthetic studies with a known planted signal, runs
both the RNA-seq and microarray branches over them, and reports how much of
the planted signal came back. It takes about a minute, needs no network
access, and touches nothing outside `06_generic_pipeline/example/output/`.

A healthy run looks roughly like this:

```
  rnaseq     up    recall  39/40  (98%)   precision 100%   0 false positive(s)
  rnaseq     down  recall  28/40  (70%)   precision 100%   0 false positive(s)
  microarray up    recall  39/40  (98%)   precision 100%   0 false positive(s)
  microarray down  recall  37/40  (92%)   precision 100%   0 false positive(s)
```

The exact numbers are not the point — recall in the 70–100% range with no
false positives is. If it recovers almost nothing, or produces long lists of
genes that were never planted, something in the R environment is wrong.

## Running it on real data

```bash
Rscript 06_generic_pipeline/run_meta_analysis.R \
    --expression my_study/counts \
    --metadata   my_study/metadata \
    --type       rnaseq \
    --case       tumour \
    --control    normal \
    --out        my_study/results
```

Always do a `--dry-run` first. It reads every file, pairs them up, checks your
group labels are actually present, and prints a table of what it found —
without starting the analysis:

```
 dataset genes samples case control other unmatched
  studyA   600      16    8       8     0         0
  studyB   600      13    6       7     0         0
  studyC   600      19   10       9     0         0
```

`unmatched` counts samples that appear in one file but not the other. A large
number there almost always means your sample IDs do not line up — GEO
accessions on one side and sample titles on the other, for instance.

`Rscript 06_generic_pipeline/run_meta_analysis.R --help` lists every option.

## Input format

**Expression** — one file per dataset, genes in rows, samples in columns. The
first column holds the gene identifier; every other column name is a sample
ID. Duplicate gene identifiers are averaged.

| | |
| --- | --- |
| `--type rnaseq` | raw integer counts |
| `--type microarray` | normalised intensities |

The distinction matters. The RNA-seq branch runs `filterByExpr`,
`calcNormFactors` and `voom` itself and cannot undo a normalisation you
already applied; the microarray branch runs kNN imputation over missing
values and leaves the scale alone. Feeding already-normalised RNA-seq data to
`--type rnaseq` produces a warning and results you should not trust.

**Metadata** — one file per dataset. One column holds the sample IDs (the
first column by default, or name it with `--id-column`) and another holds the
group labels (`--group-column`, default `state`). Any other columns are
carried through untouched.

**Pairing** — file names are matched on the part before a recognised suffix:

```
expression/studyA_counts.tsv        <->  metadata/studyA_metadata.tsv
expression/GSE12345_processed_expression.tsv
                                    <->  metadata/GSE12345_processed_sampleinfo.tsv
```

Recognised expression suffixes are `_processed_expression`,
`_expression_matrix`, `_expression`, `_counts_matrix`, `_raw_counts`,
`_counts`, `_matrix`, `_exprs`; metadata suffixes are
`_processed_sampleinfo`, `_sample_information`, `_sampleinfo`, `_phenotype`,
`_metadata`, `_samples`, `_pheno`, `_meta`, `_coldata`. Names with no
recognised suffix pair on the whole file name. Anything that fails to pair is
listed and skipped, never dropped silently.

`.tsv`, `.csv` and `.txt` are all read, and the separator is detected per
file, so you can mix them.

## Thresholds

The defaults are the ones used in the manuscript:

| Option | Default | Meaning |
| --- | --- | --- |
| `--prop-dataset` | 1 | the gene must be measured in every dataset |
| `--zval` | 2 | minimum \|Z\| of the combined effect |
| `--fdr` | 0.05 | maximum FDR |
| `--effect-size` | 0.5 | minimum \|combined effect size\| |

`--prop-dataset 1` is strict. It is reasonable when your datasets share a
platform or you have already harmonised identifiers, and punishing otherwise:
one dataset missing a gene removes it from the whole analysis. The script
tells you how many identifiers are shared by all datasets before it starts,
and refuses to run if that number is zero.

`--protein-coding` restricts the output to protein-coding genes, as the
manuscript does. It needs the `annotables` package, which is GitHub-only:
`remotes::install_github("stephenturner/annotables")`.

## Identifier harmonisation

By default DExMA's `allSameID` maps every dataset onto a common identifier
type (`--gene-id`, default `GeneSymbol`) for `--organism` (default
`Homo sapiens`). If your matrices already use one consistent identifier type,
`--no-id-harmonisation` skips that step, which is faster and avoids the
annotation packages it pulls in. The example uses the flag for exactly that
reason.

## Output

```
run_parameters.txt          what was run, with which filters, on what date
meta_analysis_all_genes.csv every gene tested, unfiltered
meta_genes_up.csv / .txt    genes passing the filters, positive effect
meta_genes_down.csv / .txt  genes passing the filters, negative effect
```

`run_parameters.txt` records the datasets used, the datasets that failed, the
thresholds, and the DExMA and R versions, so a result directory can still be
interpreted a year later.

The `.txt` files are plain gene symbols, which is what the PPI step expects:

```bash
python3 02_ppi_network/string_network.py \
    --genes my_study/results/meta_genes_up.txt --out my_study/ppi
python3 02_ppi_network/hub_genes.py \
    --edges my_study/ppi/meta_genes_up_edges.tsv --out my_study/ppi
```

## What this does not do

- **No batch correction.** Each dataset is modelled separately and combined by
  random effects; that is the point of a meta-analysis, but it is not a
  substitute for handling batch structure *within* a dataset.
- **No covariates.** The model compares two groups. If your groups differ in
  age, sex or post-mortem interval, this will not adjust for it.
- **Two groups only.** Multi-level designs need to be reduced to a pairwise
  contrast before they get here.
- **It does not check your labels.** `--case` and `--control` are taken at
  face value. Getting them backwards silently reverses every effect size, and
  nothing downstream will notice.
