#!/usr/bin/env Rscript
#
# Design qPCR primer pairs for the genes assayed in this study.
#
#     Rscript 05_primer_design/design_primer.R
#
# For each gene: fetch its template sequence from NCBI's nucleotide database by
# accession, and hand it to Primer3 through its Python bindings. Up to five
# pairs per gene are written out, and the specificity and hairpin screens that
# follow in this directory then filter them.
#
# This is not a pipeline stage. It covers all fourteen genes in the manuscript's
# primer table, but it is a reconstruction of how they were obtained, not a
# rerun of it: see README.md in this directory for exactly what it does and does
# not establish.
#
# Outputs, into results/qpcr/:
#     upreg_primers.tsv     the five up-regulated hub genes
#     downreg_primers.tsv   the five down-regulated hub genes
#     marker_primers.tsv    the reference gene and the three marker genes
#
# Requirements
# ------------
# R:      rentrez, reticulate
# Python: primer3-py, reachable from reticulate. env/python.yaml installs it;
#         point RETICULATE_PYTHON at that interpreter if reticulate picks the
#         wrong one.
# NCBI asks for a contact address with every Entrez request: set ENTREZ_EMAIL.

suppressPackageStartupMessages({
  library(rentrez)
  library(reticulate)
})

# ------------------------------------------------------------------ settings

ENTREZ_EMAIL <- Sys.getenv("ENTREZ_EMAIL", "your.name@example.com")
ORGANISM <- "Homo sapiens"

UP_GENES_FILE <- "results/ppi/intersected_genes_up_alz_hub.txt"
DOWN_GENES_FILE <- "results/ppi/intersected_genes_down_alz_hub.txt"
UP_OUTPUT <- "results/qpcr/upreg_primers.tsv"
DOWN_OUTPUT <- "results/qpcr/downreg_primers.tsv"
MARKER_OUTPUT <- "results/qpcr/marker_primers.tsv"

# The top five of each hub list are the genes carried into qPCR.
N_HUB <- 5

# ACTB is the reference gene; APP, PSEN1 and PSEN2 are the Alzheimer's marker
# genes. They are not hub genes and do not come out of stage 03b, so they are
# named here.
MARKER_GENES <- c("ACTB", "APP", "PSEN1", "PSEN2")

# Primer3 settings. Amplicon 150-250 bp, primers 18-24 nt around an optimum of
# 20, melting temperature 57-63 degrees around an optimum of 60, GC content
# 40-60% around an optimum of 50. The published primers span 153-239 bp, are
# 20-23 nt long and are 43-60% GC, so they sit inside this envelope.
#
# The nine hub genes designed through this script were originally run with the
# narrower 150-200 bp setting; the marker genes, designed in NCBI's web
# interface with otherwise identical parameters, reach 239 bp. 150-250 is the
# range that covers the table as published, and is what the manuscript reports.
DESIGN_PARAMS <- list(
  PRIMER_PRODUCT_SIZE_RANGE = as.integer(c(150, 250)),
  PRIMER_OPT_SIZE = 20L,
  PRIMER_MIN_SIZE = 18L,
  PRIMER_MAX_SIZE = 24L,
  PRIMER_OPT_TM = 60.0,
  PRIMER_MIN_TM = 57.0,
  PRIMER_MAX_TM = 63.0,
  PRIMER_MIN_GC = 40.0,
  PRIMER_OPT_GC_PERCENT = 50.0,
  PRIMER_MAX_GC = 60.0
)

# ------------------------------------------------------------------ templates
#
# One pinned accession per gene. Pinning matters: the original run asked Entrez
# for `<gene>[Gene] AND Homo sapiens[Organism]` and took whatever came back
# first, and that ranking is not stable, not guaranteed to be a transcript, and
# not guaranteed to be the right gene. Each accession below was established by
# locating the primer pair for that gene inside the record rather than by
# assuming, so the template every pair was designed on is documented exactly.
#
TEMPLATES <- c(
  # up-regulated hub genes
  HDAC1   = "NG_047018.1",     # RefSeqGene, genomic
  CREBBP  = "NG_009873.2",     # RefSeqGene, genomic
  NFKB1   = "NM_001319226.2",  # transcript variant 3, mRNA
  KAT2B   = "NM_003884.5",     # mRNA
  SP1     = "NM_138473.3",     # transcript variant 1, mRNA

  # down-regulated hub genes
  UBB     = "NM_001281720.2",  # transcript variant 6, mRNA
  ATP5F1A = "NG_041769.2",     # RefSeqGene, genomic
  PRKACA  = "NM_001304349.2",  # transcript variant 3, mRNA
  PRKACB  = "NM_001242860.3",  # transcript variant 7, mRNA
  PSMD14  = "NM_005805.6",     # mRNA; the genomic record is a whole chromosome

  # reference gene and Alzheimer's marker genes
  ACTB    = "NM_001101.5",     # mRNA
  APP     = "NM_000484.4",     # transcript variant 1, mRNA
  PSEN1   = "NM_000021.4",     # transcript variant 1, mRNA
  PSEN2   = "NM_000447.3"      # transcript variant 1, mRNA
)

# Exon-exon junctions, as positions in the pinned transcript, for genes whose
# primers are required to straddle one. A pair that spans a junction cannot
# amplify genomic DNA carry-over, so the assay measures transcript and not
# contaminating template.
#
# Only genes listed here get the constraint. It is not applied everywhere
# because three of the templates above (HDAC1, CREBBP, ATP5F1A) are genomic
# records, which have no junctions to straddle.
#
# The positions are derived from the transcript, not typed in by hand:
# junctions_from_genomic() below aligns the transcript against the gene's
# RefSeqGene record and returns the block boundaries.
JUNCTION_TEMPLATES <- list(
  SP1 = list(refseqgene = "NG_030361.1")
)

#' Recover exon-exon junction positions for a transcript.
#'
#' Walks the transcript, anchoring successive 30-mers in the genomic record and
#' extending each match as far as it goes. Every place the walk has to jump
#' forward in the genomic sequence is an intron, and the transcript position
#' where it jumped is a junction.
junctions_from_genomic <- function(transcript, genomic, anchor_length = 30L) {
  boundaries <- integer(0)
  position <- 1L
  search_from <- 1L
  n <- nchar(transcript)
  while (position <= n - anchor_length) {
    anchor <- substr(transcript, position, position + anchor_length - 1L)
    hit <- regexpr(anchor, substr(genomic, search_from, nchar(genomic)), fixed = TRUE)
    if (hit < 0) break
    g <- search_from + hit - 1L
    i <- 0L
    while (position + i <= n && g + i <= nchar(genomic) &&
           substr(transcript, position + i, position + i) ==
           substr(genomic, g + i, g + i)) i <- i + 1L
    position <- position + i
    search_from <- g + i
    if (position <= n) boundaries <- c(boundaries, position - 1L)
  }
  boundaries
}

# ------------------------------------------------------- sequence retrieval

#' Fetch one accession from NCBI's nucleotide database.
fetch_accession <- function(gene, accession, email = ENTREZ_EMAIL) {
  fasta <- tryCatch(
    entrez_fetch(db = "nucleotide", id = accession, rettype = "fasta", email = email),
    error = function(e) {
      message(sprintf("  %s (%s): fetch failed - %s", gene, accession, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fasta)) return(NULL)

  lines <- unlist(strsplit(fasta, "\n", fixed = TRUE))
  sequence <- paste(lines[-1], collapse = "")
  message(sprintf("  %s: %s (%d bp)", gene, sub("^>", "", lines[1]), nchar(sequence)))
  sequence
}

#' Fall back to a search for a gene with no pinned accession.
#'
#' Only reached for genes that are not in TEMPLATES, which means genes this
#' study did not assay. The result depends on how Entrez ranks records on the
#' day, so it is reported loudly rather than silently.
fetch_top_hit <- function(gene, organism = ORGANISM, email = ENTREZ_EMAIL) {
  term <- sprintf("%s[Gene] AND %s[Organism]", gene, organism)
  search <- tryCatch(
    entrez_search(db = "nucleotide", term = term, retmax = 1, email = email),
    error = function(e) {
      message(sprintf("  %s: search failed - %s", gene, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(search) || length(search$ids) == 0) {
    message(sprintf("  %s: no nucleotide record found", gene))
    return(NULL)
  }
  message(sprintf("  %s: NOT PINNED - using the current Entrez top hit. Check it", gene))
  message("      is the right gene before ordering anything designed on it.")
  fetch_accession(gene, search$ids[1], email = email)
}

fetch_all_sequences <- function(genes) {
  message(sprintf("Retrieving %d sequence(s) from NCBI Entrez ...", length(genes)))
  sequences <- list()
  for (gene in genes) {
    sequence <- if (gene %in% names(TEMPLATES)) {
      fetch_accession(gene, TEMPLATES[[gene]])
    } else {
      fetch_top_hit(gene)
    }
    if (!is.null(sequence)) sequences[[gene]] <- sequence
  }
  sequences
}

# ---------------------------------------------------------- primer design

primer3 <- NULL

#' Design primer pairs for a named list of sequences.
#'
#' Primer3 is reached through its Python bindings rather than the command-line
#' binary, which is why reticulate is a dependency.
design_primers <- function(sequences, params = DESIGN_PARAMS) {
  if (is.null(primer3)) {
    primer3 <<- tryCatch(
      import("primer3"),
      error = function(e) {
        stop("Could not import primer3 from Python. Install primer3-py (it is in\n",
             "env/python.yaml) and, if reticulate is finding the wrong interpreter,\n",
             "set RETICULATE_PYTHON to the one that has it.\n  ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }

  designs <- list()
  for (gene in names(sequences)) {
    sequence <- sequences[[gene]]
    length_bp <- nchar(sequence)
    if (length_bp < params$PRIMER_PRODUCT_SIZE_RANGE[1]) {
      message(sprintf("  %s: sequence shorter than the minimum amplicon, skipped", gene))
      next
    }
    seq_args <- list(
      SEQUENCE_ID = gene,
      SEQUENCE_TEMPLATE = sequence,
      SEQUENCE_INCLUDED_REGION = as.integer(c(0, length_bp - 1))
    )

    # Genes listed in JUNCTION_TEMPLATES are constrained to straddle an
    # exon-exon junction, so their pairs cannot amplify genomic DNA.
    gene_params <- params
    if (gene %in% names(JUNCTION_TEMPLATES)) {
      genomic <- fetch_accession(
        paste0(gene, " (genomic, for junctions)"),
        JUNCTION_TEMPLATES[[gene]]$refseqgene
      )
      junctions <- if (is.null(genomic)) integer(0) else
        junctions_from_genomic(sequence, genomic)
      if (length(junctions) == 0) {
        message(sprintf("  %s: no junctions recovered; designing without the constraint", gene))
      } else {
        message(sprintf("  %s: %d exon-exon junction(s) at %s", gene,
                        length(junctions), paste(junctions, collapse = ", ")))
        seq_args$SEQUENCE_OVERLAP_JUNCTION_LIST <- as.integer(junctions)
        gene_params$PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION <- 7L
        gene_params$PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION <- 4L
      }
    }

    designs[[gene]] <- tryCatch(
      primer3$bindings$design_primers(seq_args, gene_params),
      error = function(e) {
        message(sprintf("  %s: Primer3 failed - %s", gene, conditionMessage(e)))
        NULL
      }
    )
  }
  designs
}

#' Flatten Primer3's nested output into one row per primer pair.
collect_primers <- function(designs) {
  rows <- list()
  for (gene in names(designs)) {
    primers <- designs[[gene]]
    if (is.null(primers) || is.null(primers$PRIMER_PAIR) || length(primers$PRIMER_PAIR) == 0) {
      message(sprintf("  %s: no valid primer pairs", gene))
      next
    }
    for (i in seq_along(primers$PRIMER_PAIR)) {
      left <- primers$PRIMER_LEFT[[i]]
      right <- primers$PRIMER_RIGHT[[i]]
      rows[[length(rows) + 1]] <- data.frame(
        Gene = gene,
        Pair = i,
        LeftPrimer = left$SEQUENCE,
        RightPrimer = right$SEQUENCE,
        ProductSize = primers$PRIMER_PAIR[[i]]$PRODUCT_SIZE,
        LeftTm = left$TM,
        RightTm = right$TM,
        LeftGC = left$GC_PERCENT,
        RightGC = right$GC_PERCENT,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(data.frame(
      Gene = character(), Pair = integer(), LeftPrimer = character(),
      RightPrimer = character(), ProductSize = integer(), LeftTm = numeric(),
      RightTm = numeric(), LeftGC = numeric(), RightGC = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

# -------------------------------------------------------------------- main

read_hub_genes <- function(path, n) {
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run stage 03b first.", call. = FALSE)
  }
  genes <- read.delim(path, header = TRUE)$gene
  head(genes, n)
}

run_group <- function(label, genes, output) {
  message("\n== ", label, " ==")
  message("Genes: ", paste(genes, collapse = ", "))

  sequences <- fetch_all_sequences(genes)
  if (length(sequences) == 0) {
    message("No sequences retrieved; nothing written for ", label, ".")
    return(invisible(NULL))
  }

  message("Designing primers ...")
  primers <- collect_primers(design_primers(sequences))
  if (nrow(primers) == 0) {
    message("No primer pairs designed; nothing written for ", label, ".")
    return(invisible(NULL))
  }

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  write.table(primers, file = output, sep = "\t")
  message(sprintf("Wrote %s (%d pairs across %d genes)",
                  output, nrow(primers), length(unique(primers$Gene))))
}

if (identical(ENTREZ_EMAIL, "your.name@example.com")) {
  message("Note: ENTREZ_EMAIL is unset, so requests carry a placeholder address.\n",
          "NCBI asks for a real one and may throttle or block anonymous traffic.\n")
}

run_group("Up-regulated hub genes",
          read_hub_genes(UP_GENES_FILE, N_HUB), UP_OUTPUT)
run_group("Down-regulated hub genes",
          read_hub_genes(DOWN_GENES_FILE, N_HUB), DOWN_OUTPUT)
run_group("Reference and marker genes",
          MARKER_GENES, MARKER_OUTPUT)

message("\nNext: primer_specificity_blastn.py, then hairpin_check.py")
