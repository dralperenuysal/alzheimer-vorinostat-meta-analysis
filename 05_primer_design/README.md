# Primer design

**This directory is not part of the pipeline.** `run_all.sh` does not run it, and
`verify_reproduction.py` does not check it. Primer design was a one-off step between the
computational half of the study and the bench half: the hub genes came out of
`02_ppi_network/`, primers were designed for the top five of each direction plus the
reference and marker genes, and the qPCR was run with them. The
scripts here cover all fourteen genes in the manuscript's primer table, and the template
accession for every gene is pinned in `TEMPLATES`, so the sequence each pair was designed
on is documented exactly.

## What is here

| File | What it does |
| --- | --- |
| `design_primer.R` | The design route: Entrez retrieval by pinned accession, then Primer3 through its Python bindings |
| `primer_specificity_blastn.py` | BLAST screen of the candidates against the human genome |
| `hairpin_check.py` | Hairpin and self-dimer screen over the survivors |
| `pcr_functions.py` | Helpers shared by the Python scripts |
| `design_primer.py` | An earlier, incomplete pass; kept for the record, never used (see its docstring) |

## Running it

Nothing here is needed to reproduce the analysis. If you want to design primers for your
own genes, `env/python.yaml` already carries `primer3-py` and `biopython`; the two R
packages are not in `env/r.yaml`, because the pipeline does not need them:

```r
install.packages(c("rentrez", "reticulate"))
```

`reticulate` has to find the interpreter that has `primer3-py` — set `RETICULATE_PYTHON`
to it if it picks the wrong one.

```bash
export ENTREZ_EMAIL="you@example.org"   # NCBI asks for a contact address
Rscript 05_primer_design/design_primer.R
```

`DESIGN_PARAMS` sets a 150–250 bp amplicon, 18–24 nt primers around an optimum of 20,
Tm 57–63 °C around an optimum of 60, and 40–60% GC around an optimum of 50.

A gene with no pinned accession falls back to the old top-hit search and says so loudly.
That ranking is not stable and is not guaranteed to return the gene you asked for, so
check what it retrieved before ordering anything designed on it.

The script writes every candidate pair Primer3 returns. Choosing among them is the
specificity and hairpin screens' job, and then yours.
