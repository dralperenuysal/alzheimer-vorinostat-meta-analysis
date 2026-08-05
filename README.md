# Reversing Alzheimer's Disease Transcriptomic Signatures with the HDAC Inhibitor Vorinostat

Analysis code for the manuscript

> Uysal SA, Rostamlou A, Tezcanlı Kaymaz B. *Reversing Alzheimer's Disease
> Transcriptomic Signatures with the HDAC Inhibitor Vorinostat: An In-Vitro
> Validation Study.*

The study has two halves. A random-effects meta-analysis of 9 microarray and 11
bulk RNA-seq GEO series nominates genes consistently dysregulated in Alzheimer's
disease; protein–protein interaction network analysis prioritises hub genes among
them; and the prioritised genes are then measured by qPCR in differentiated
SH-SY5Y cells under hydrogen-peroxide stress with and without vorinostat.

The pipeline starts from **prepared expression matrices and sample metadata**.
Retrieving the series from GEO, annotating probes and labelling samples are
assumed done; the next section gives the exact format expected. From there,
`run_all.sh` runs the whole chain, and `verify_reproduction.py` diffs the
result against what is committed.

```bash
bash run_all.sh                 # list the stages
bash run_all.sh all             # run the chain end to end
python3 verify_reproduction.py  # diff the re-run against the published results

Rscript 06_generic_pipeline/run_meta_analysis.R --help   # your own data
```

Directories 01-04 are the pipeline stages, numbered in the order they run and
each corresponding to a subsection of the manuscript's Methods. 05 and 06 are
not stages: primer design was a one-off bench step, and the generic pipeline is
the meta-analysis with the Alzheimer's specifics removed.

## Layout

| Directory | Methods section | Contents |
| --- | --- | --- |
| `01_meta_analysis/` | 2.2 | DExMA random-effects meta-analysis and the microarray/RNA-seq intersection |
| `02_ppi_network/` | 2.3–2.4 | STRING network construction and hub gene ranking by degree centrality |
| `03_druggability/` | 2.5 Target Druggability Assessment | ChEMBL ligand retrieval for the hub genes, restricted to IC50 with a pChEMBL value |
| `04_figures/` | 2.11–2.12 | ΔCt / ΔΔCt relative expression and Manuscript Figures 1–5, 8–11 |
| `05_primer_design/` | 2.11 Real-Time PCR | Primer3 design, BLAST specificity check, hairpin screening. **Not a pipeline stage**, see below |
| `06_generic_pipeline/` | — | The meta-analysis, generalised to arbitrary input datasets |
| `results/` | — | Committed outputs (see below) |
| `env/` | — | Conda environment specifications |
| `run_all.sh` | — | Stage runner for the whole pipeline |
| `verify_reproduction.py` | — | Diffs a re-run against the published results |

**`05_primer_design/` is deliberately not a stage.** The chain this repository
reproduces runs from the expression matrices to the hub genes, and on to the
qPCR analysis and the figures. Primer design happened once, between the two
halves. The scripts are here as the record of how it was done, and they cover
all fourteen genes in the manuscript's primer table, but they are not re-run
and not verified. See `05_primer_design/README.md`.

## Input data

The two DExMA scripts in `01_meta_analysis/` read one directory per dataset
under `data/`, which is not committed — the twenty matrices come to roughly
300 MB, none of it original to this study, and it is all rebuildable from the
GEO accessions listed below. Nothing else in the pipeline needs it: every later
stage starts from files that are committed under `results/`.

Prepare `data/` in this layout:

```
data/
├── microarray_processed/
│   └── GSE110226/
│       ├── GSE110226_processed_expression.tsv
│       └── GSE110226_processed_sampleinfo.tsv
└── rnaseq_processed/
    └── GSE153873/
        ├── GSE153873_processed_expression.tsv
        └── GSE153873_processed_sampleinfo.tsv
```

The directory name is the dataset identifier and is repeated in both file
names. Which directories are read is driven by
`results/dataset_selection/{microarray,rnaseq}_check.tsv` — the rows marked
`suitable`. Use your own identifiers there if you are not reproducing this
study.

Both branches expect the same two shapes; only the values differ.

**Expression.** Tab-separated, genes in rows, samples in columns. The first
column is headed exactly `gene`; every other column name is a sample ID.

```
gene	GSM4649255	GSM4649256	GSM4649257
HDAC1	412	388	455
CREBBP	1204	1150	1331
NFKB1	87	102	95
```

For `rnaseq_processed/` these are **raw counts** — do not normalise, the script
runs `filterByExpr`, `calcNormFactors` and `voom` itself. For
`microarray_processed/` they are **normalised log2 intensities**; missing
values are imputed by kNN, so blanks are allowed.

**Sample metadata.** One row per sample. The first column is headed exactly
`samples` and holds the sample IDs; a `state` column holds exactly `AD` or
`control`. Any other columns are carried through and ignored.

```
samples	state
GSM4649255	AD
GSM4649256	AD
GSM4649257	control
```

Every sample named in the metadata must exist as a column in the expression
matrix. Samples the matrix has and the metadata does not are dropped.

**Gene identifiers are the join key**, so sort them out before this point. The
`gene` column may hold HGNC symbols or Ensembl gene IDs — `allSameID` maps every
dataset onto `GeneSymbol` — but it must hold **one identifier type per
dataset**, and probe IDs are not one of the two. The meta-analysis then keeps
only the symbols shared across all datasets, so unmapped probes, outdated
aliases and platform-specific identifiers never fail loudly; they quietly
shrink the gene space, and a single badly annotated dataset can shrink it a
long way. On every dataset: resolve probes to symbols or Ensembl IDs, collapse
duplicate rows to one row per identifier, and drop rows that map to nothing.

## Datasets

The 20 GEO series behind the meta-analysis. Selection was manual; the full
candidate list with the reason each series was kept or dropped is in
`results/dataset_selection/`.

**Microarray (9):** GSE110226, GSE29378, GSE61196, GSE122063, GSE15222,
GSE109887, GSE150696, GSE39420, GSE138260

**Bulk RNA-seq (11):** GSE148822, GSE125583, GSE125050, GSE153873, GSE226901,
GSE95587, GSE104704, GSE113524, GSE163877, GSE203206, GSE261817

## Running it

```bash
conda env create -f env/r.yaml
conda env create -f env/python.yaml
```

Then either run everything,

```bash
bash run_all.sh all
```

or pick up at a stage, which is what you want if you only care about part of
the chain:

```bash
bash run_all.sh from 02     # PPI step onwards
bash run_all.sh 02 04       # just the PPI step and the figures
```

`run_all.sh` reports each script as it runs and keeps going after a failure;
anything that failed is listed again at the end.

Scripts expect to be run **from the repository root**, not from their own
directory, whether through `run_all.sh` or by hand:

```bash
Rscript 01_meta_analysis/DExMA_meta_microarray.R
python3 04_figures/make_qpcr_figures.py
```

`data/` is deliberately not committed — none of it is original to this study,
and it runs to roughly a gigabyte. Everything from stage 01 onwards reads and
writes `results/`, which is committed, so the meta-analysis and everything
downstream can be inspected, and re-run, without rebuilding the input.

Two stages need network access: `02` queries STRING and `03` queries ChEMBL,
both public and both needing no credentials.

## Checking a re-run

Re-running overwrites `results/`. Rather than asking you to back it up first,
`verify_reproduction.py` reads the published versions straight out of git
history and diffs the working tree against them:

```bash
python3 verify_reproduction.py
```

It compares the per-platform gene lists, the combined effect sizes, the
intersection, the hub genes and the qPCR figures, and reports `PASS`,
`DIFFERS` or `SKIP` for each. `SKIP` means the file is still byte-identical to
what was committed — that stage has not actually been re-run yet, and calling
it reproduced would be overstating things.

**How much of it you can check without the input matrices.** Only the two
DExMA scripts in `01_meta_analysis/` read `data/`. Everything after them starts
from files that are committed, so without preparing any input at all you can
re-run the intersection, the STRING network and hub gene ranking, the ChEMBL
query, the ΔCt calculation and all four qPCR figures, and diff every one of
them against the published version. Half the checks below run for real on a
bare clone; the other half — the per-platform gene lists and effect sizes,
which are the direct output of the DExMA step — report `SKIP` until you supply
`data/` yourself. The published values for those are committed in
`results/meta_analysis/`, so they can be read and inspected, just not
regenerated from scratch.

Two comparisons are deliberately not byte-exact:

- **Figures.** matplotlib writes a fresh `/CreationDate` into every PDF, so
  two runs of the same script never produce identical bytes. The check
  compares the extracted text and the file size instead.
- **Hub genes.** Degree is an integer over a few hundred proteins, so a top-N
  cut can fall inside a tied group and two runs return different genes while
  ranking the network identically. The current cut, 14, was chosen so this
  does not happen — see below — but the check still verifies it rather than
  assuming it, and reports the degree of every swapped gene if a future rerun
  ever does land on a tie.

`DIFFERS` is not automatically a failure. Package versions move, GEO
occasionally revises a series, and STRING publishes new releases; the notes
under each result are there to let you judge which it is.

## What is in `results/`

| Path | What it is |
| --- | --- |
| `dataset_selection/` | Every candidate GEO series with its inclusion verdict (`suitable`, `blood`, `iPSC`). **Input**: the `suitable` rows drive which directories stage 01 reads |
| `meta_analysis/{microarray,rnaseq}_{up,down}.csv` | Full DExMA output — combined effect size, Z, p, FDR per gene. The volcano plots (Figure 1) are drawn from these |
| `meta_analysis/*_meta_genes_*.txt` | Gene symbols passing the thresholds, per platform and direction |
| `meta_analysis/intersected_genes_{up,down}_alz.txt` | Microarray ∩ RNA-seq — the input to the PPI network (Figure 2) |
| `ppi/*_edges.tsv` | Every STRING interaction at score > 0.7 among the intersected genes, with all evidence subscores |
| `ppi/*_network.sif` | The same network, ready to import into Cytoscape |
| `ppi/*_degree.tsv` | Every connected protein ranked by degree, with the size of its tie group |
| `ppi/*_nodes.tsv` | STRING's identifier mapping, so any gene it could not resolve is visible |
| `ppi/*_hub.txt` | The 14 hub genes per direction by degree centrality (Figure 4). Ten of these, five per direction, were carried into qPCR |
| `qpcr/supplementary_qpcr_data.csv` | One row per well with ΔCt and a detection flag; the manuscript's supplementary file. **Input**: every qPCR figure and log2FC is computed from it |
| `qpcr/*_primers*.tsv` | Designed primers before and after hairpin screening |
| `druggability/hub_gene_ligands.tsv` | One row per hub gene, ChEMBL target and ligand, with the pChEMBL potency and the retrieval date |
| `druggability/hub_gene_ligand_counts.tsv` | Distinct ligands per hub gene — the evidence behind the target selection |
| `figures/Fig{1,2}.pdf` | Volcano plots and platform overlap, drawn by `04_figures/{volcano_plots,venn_diagrams}.R` |
| `figures/Fig{3,4}.pdf` | PPI networks and hub gene subnetworks, drawn by `04_figures/make_network_figures.py` |
| `figures/Fig5.pdf` | Target selection: centrality against ligand count, and where vorinostat sits among HDAC1's ligands |
| `figures/Fig{8,9,10,11}.pdf` | qPCR panels, drawn by `04_figures/make_qpcr_figures.py`; byte-identical in content to the submitted figures |

`string_network.py --image` additionally downloads STRING's own drawing of the
network as SVG. Those files are not committed — STRING embeds every protein
bubble as a base64 image, which comes to about 35 MB for the two networks here.

## The PPI and hub gene step

For the manuscript, both steps were done inside Cytoscape v3.10.3: the network
was built with the STRING plug-in (stringApp), which queries the STRING v12.0
database, and the hub genes were ranked with the CytoHubba plugin. Neither step
left a script behind. `02_ppi_network/` reconstructs both without the GUI, and
goes to the same STRING backend the plug-in uses:

```bash
python3 02_ppi_network/string_network.py \
    --genes results/meta_analysis/intersected_genes_up_alz.txt --out results/ppi
python3 02_ppi_network/hub_genes.py \
    --edges results/ppi/intersected_genes_up_alz_edges.tsv --out results/ppi
```

`string_network.py` queries STRING's REST API, pinned to
`version-12-0.string-db.org` so the network does not drift when STRING
publishes a new release. `hub_genes.py` counts edges per node, which is exactly
what CytoHubba's Degree method does — no layout and no plugin are involved in
the ranking, only in the drawing. Python's standard library is all it needs.

**Which genes went into the qPCR, and why they are not the top five.** The ten
genes measured by qPCR were prioritised on an interim version of the
meta-analysis, before two microarray series were found to violate the inclusion
criteria and removed — GSE138261, a SuperSeries whose samples duplicate the
separately included GSE138260, and GSE117586, an iPS-cell model rather than
brain tissue. Both are marked as excluded in
`results/dataset_selection/microarray_check.tsv`. Removing them reordered the
hub genes without changing which genes are hubs: all ten remain in the hub gene
lists of their networks.

| Direction | Gene | Degree | Rank |
| --- | --- | ---: | ---: |
| up | HDAC1 | 17 | 1 |
| up | CREBBP | 11 | 2 (tied with 1) |
| up | NFKB1 | 11 | 2 (tied with 1) |
| up | KAT2B | 7 | 8 (tied with 1) |
| up | SP1 | 6 | 10 (tied with 4) |
| down | ATP5F1A | 30 | 1 |
| down | UBB | 28 | 2 |
| down | PSMD14 | 24 | 5 |
| down | PRKACA | 22 | 7 (tied with 3) |
| down | PRKACB | 22 | 7 (tied with 3) |

The networks are 185 interactions among 146 connected proteins (up) and 1,109
among 402 (down), at STRING combined score > 0.7.

**Ties are common, which is why the cut is 14 and not 10.** Degree is an
integer over a few hundred proteins, so a top-N boundary easily falls inside a
group of genes with the same degree — at 10 it would have split a five-way tie
in the upregulated network and a four-way tie in the downregulated one. At 14
both cuts land on a group boundary: degree drops 6 → 5 (up) and 21 → 19 (down)
between ranks 14 and 15, so no protein is admitted or excluded arbitrarily.
`hub_genes.py` breaks any remaining tie alphabetically, prints a warning
whenever a cut falls inside a tied group, and `results/ppi/*_degree.tsv` lists
the tie size for every protein, so none of this has to be taken on trust.

Figures 3 and 4 are drawn from these tables by
`04_figures/make_network_figures.py`, which reads the degrees rather than
recomputing them, so the sizes and colours are the same numbers as the ranking
above. The `.sif` files reload the identical network into Cytoscape if you want
to lay it out by hand instead.

## Steps that still have no code

- **Figure 7** — light microscopy images of differentiated and undifferentiated
  cells.

`results/figures/` therefore holds Figures 1–5 and the four qPCR figures. Figure 6 is drawn
by a script kept with the manuscript source, from the raw plate readings kept
alongside it, since it also embeds a panel that is not regenerable.

## Using the method on your own data

`06_generic_pipeline/` is the meta-analysis of `01_meta_analysis/` with the
Alzheimer's specifics removed. Give it a directory of expression matrices and a
directory of sample metadata, and it pairs them by file name, normalises each
dataset by assay type, and runs the same DExMA model under the same filters.

Like `01_meta_analysis/`, it assumes your data is already prepared: it
downloads nothing and curates nothing. Annotating probes, harmonising gene
symbols and deciding which sample belongs to which group are yours to do first.
The difference is that it takes flat directories of files and a pair of group
labels on the command line, rather than the fixed layout and `AD`/`control`
labels described under **Input data** above.

```bash
bash 06_generic_pipeline/example/run_example.sh      # one-minute smoke test

Rscript 06_generic_pipeline/run_meta_analysis.R \
    --expression my_study/counts --metadata my_study/metadata \
    --type rnaseq --case tumour --control normal --out my_study/results
```

Both assay branches are covered: `--type rnaseq` for raw counts, which runs
`filterByExpr`/`calcNormFactors`/`voom`, and `--type microarray` for normalised
intensities, which runs kNN imputation. Output gene lists feed straight into
`02_ppi_network/`. `06_generic_pipeline/README.md` has the full description,
including what the pipeline deliberately does not handle.

## Software

R 4.3.3, Python 3.12. Key packages: DExMA 1.10.7 (meta-analysis), limma 3.58.1
(differential expression), edgeR 4.0.16 (RNA-seq normalisation), primer3-py and
Primer3 (primer design), chembl_webresource_client (druggability), matplotlib
(figures). Full lists in `env/`.

`02_ppi_network/` and `verify_reproduction.py` use only the Python standard
library, so they run without either conda environment. `verify_reproduction.py`
uses `pdftotext` (poppler-utils) for the figure comparison and skips that check
if it is not installed.

`05_primer_design/` needs `rentrez`, `reticulate` and `primer3-py`, which are
not in `env/` because the pipeline does not run it. Its own README lists them.

## License

MIT — see `LICENSE`. The GEO datasets analysed here remain under their original
terms; accessions are listed above.
