#!/usr/bin/env Rscript
#
# Random-effects transcriptomic meta-analysis over your own datasets.
#
# This is the meta-analysis of 01_meta_analysis/ with the Alzheimer's specifics
# taken out. You supply a directory of expression matrices and a directory of
# sample metadata; it pairs them up by file name, normalises each dataset
# according to its assay type, and runs the same DExMA random-effects model
# under the same significance filters.
#
#   Rscript 06_generic_pipeline/run_meta_analysis.R \
#     --expression my_study/expression \
#     --metadata   my_study/metadata \
#     --type       rnaseq \
#     --case       tumour \
#     --control    normal \
#     --out        my_study/results
#
# Add --dry-run to check that your files pair up and your group labels are
# present before committing to a long run. It reads every file and reports
# what it found, but stops before the meta-analysis.
#
# Input format
# ------------
# Expression: one file per dataset, genes in rows, samples in columns. First
# column holds the gene identifier, remaining column names are sample IDs.
# Microarray files should hold normalised intensities; RNA-seq files should
# hold raw integer counts, because this script does its own edgeR/voom
# normalisation and cannot undo a previous one.
#
# Metadata: one file per dataset. One column holds the sample IDs (by default
# the first column) and another holds the group label (--group-column, default
# "state"). Extra columns are carried through untouched.
#
# Pairing: file names are matched on the part before a recognised suffix, so
# GSE12345_processed_expression.tsv pairs with GSE12345_processed_sampleinfo.tsv,
# and study_a_counts.csv pairs with study_a_metadata.csv. Anything unmatched is
# reported and skipped, never silently dropped.
#
# .tsv, .csv and .txt are all accepted; the separator is detected per file.

suppressPackageStartupMessages({
  library(Biobase)
  library(DExMA)
})

# ---------------------------------------------------------------- arguments

DEFAULTS <- list(
  expression = NULL,
  metadata = NULL,
  type = NULL,
  out = "meta_analysis_results",
  `group-column` = "state",
  `id-column` = NULL,
  case = NULL,
  control = NULL,
  `prop-dataset` = "1",
  zval = "2",
  fdr = "0.05",
  `effect-size` = "0.5",
  organism = "Homo sapiens",
  `gene-id` = "GeneSymbol",
  datasets = NULL,
  `impute-k` = "7"
)

FLAGS <- c("dry-run", "no-id-harmonisation", "protein-coding", "help")

usage <- function() {
  # The header comment above is the help text. Read it back off disk rather
  # than duplicating it in a string, stopping at the first line of code.
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  if (!is.na(script) && file.exists(script)) {
    lines <- readLines(script, warn = FALSE)
    lines <- lines[-1]                                   # drop the shebang
    lines <- lines[seq_len(which(!startsWith(lines, "#"))[1] - 1)]
    cat(sub("^# ?", "", lines), sep = "\n")
  }
  cat("\nOptions:\n")
  cat("  --expression DIR      directory of expression matrices        [required]\n")
  cat("  --metadata DIR        directory of sample metadata files      [required]\n")
  cat("  --type TYPE           'rnaseq' (raw counts) or 'microarray'   [required]\n")
  cat("  --case LABEL          group label for the condition of interest [required]\n")
  cat("  --control LABEL       group label for the reference group     [required]\n")
  cat("  --out DIR             output directory                        [meta_analysis_results]\n")
  cat("  --group-column NAME   metadata column holding the labels      [state]\n")
  cat("  --id-column NAME      metadata column holding sample IDs      [first column]\n")
  cat("  --datasets A,B,C      restrict to these dataset IDs           [all]\n")
  cat("\n  --prop-dataset F      min fraction of datasets a gene appears in [1]\n")
  cat("  --zval F              min |Z| of the combined effect          [2]\n")
  cat("  --fdr F               max FDR                                 [0.05]\n")
  cat("  --effect-size F       min |combined effect size|              [0.5]\n")
  cat("\n  --organism NAME       for identifier harmonisation           [Homo sapiens]\n")
  cat("  --gene-id ID          target identifier type                  [GeneSymbol]\n")
  cat("  --no-id-harmonisation skip DExMA's allSameID step\n")
  cat("  --protein-coding      keep protein-coding genes only (needs annotables)\n")
  cat("  --impute-k N          k for microarray kNN imputation         [7]\n")
  cat("  --dry-run             validate inputs and stop\n")
  cat("  --help                show this message\n\n")
  quit(status = 0)
}

parse_args <- function(argv) {
  opts <- DEFAULTS
  for (flag in FLAGS) opts[[flag]] <- FALSE
  i <- 1
  while (i <= length(argv)) {
    token <- argv[i]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    key <- substring(token, 3)
    if (key %in% FLAGS) {
      opts[[key]] <- TRUE
      i <- i + 1
    } else if (key %in% names(DEFAULTS)) {
      if (i == length(argv)) stop("--", key, " needs a value", call. = FALSE)
      opts[[key]] <- argv[i + 1]
      i <- i + 2
    } else {
      stop("Unknown option --", key, ". Run with --help for the list.", call. = FALSE)
    }
  }
  opts
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0) usage()
opts <- parse_args(argv)
if (isTRUE(opts$help)) usage()

required <- c("expression", "metadata", "type", "case", "control")
missing_args <- required[vapply(required, function(k) is.null(opts[[k]]), logical(1))]
if (length(missing_args)) {
  stop("Missing required option(s): --", paste(missing_args, collapse = ", --"),
       "\nRun with --help for usage.", call. = FALSE)
}
if (!opts$type %in% c("rnaseq", "microarray")) {
  stop("--type must be 'rnaseq' or 'microarray', not '", opts$type, "'", call. = FALSE)
}
for (dir_opt in c("expression", "metadata")) {
  if (!dir.exists(opts[[dir_opt]])) {
    stop("--", dir_opt, " directory does not exist: ", opts[[dir_opt]], call. = FALSE)
  }
}

as_number <- function(key) {
  value <- suppressWarnings(as.numeric(opts[[key]]))
  if (is.na(value)) stop("--", key, " must be a number, got '", opts[[key]], "'", call. = FALSE)
  value
}
prop_dataset <- as_number("prop-dataset")
min_zval <- as_number("zval")
max_fdr <- as_number("fdr")
min_effect <- as_number("effect-size")
impute_k <- as_number("impute-k")

# ------------------------------------------------------------ file pairing

# Suffixes stripped when deriving a dataset ID from a file name. Longest
# first, so that "_processed_expression" wins over "_expression".
EXPRESSION_SUFFIXES <- c("_processed_expression", "_expression_matrix", "_expression",
                         "_counts_matrix", "_raw_counts", "_counts", "_matrix", "_exprs")
METADATA_SUFFIXES <- c("_processed_sampleinfo", "_sample_information", "_sampleinfo",
                       "_phenotype", "_metadata", "_samples", "_pheno", "_meta", "_coldata")

dataset_id <- function(path, suffixes) {
  stem <- tools::file_path_sans_ext(basename(path))
  for (suffix in suffixes) {
    if (endsWith(tolower(stem), tolower(suffix))) {
      return(substring(stem, 1, nchar(stem) - nchar(suffix)))
    }
  }
  stem
}

index_directory <- function(dir, suffixes) {
  files <- list.files(dir, pattern = "\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)
  index <- list()
  for (file in files) {
    id <- dataset_id(file, suffixes)
    if (!is.null(index[[id]])) {
      stop("Two files in ", dir, " both map to dataset '", id, "':\n  ",
           basename(index[[id]]), "\n  ", basename(file),
           "\nRename one of them.", call. = FALSE)
    }
    index[[id]] <- file
  }
  index
}

expression_files <- index_directory(opts$expression, EXPRESSION_SUFFIXES)
metadata_files <- index_directory(opts$metadata, METADATA_SUFFIXES)

if (!length(expression_files)) {
  stop("No .tsv/.csv/.txt files in ", opts$expression, call. = FALSE)
}

paired <- intersect(names(expression_files), names(metadata_files))
if (!is.null(opts$datasets)) {
  wanted <- trimws(strsplit(opts$datasets, ",")[[1]])
  unknown <- setdiff(wanted, paired)
  if (length(unknown)) {
    stop("--datasets names ", paste(unknown, collapse = ", "),
         ", which did not pair up. Available: ", paste(paired, collapse = ", "), call. = FALSE)
  }
  paired <- wanted
}
paired <- sort(paired)

orphan_expression <- setdiff(names(expression_files), names(metadata_files))
orphan_metadata <- setdiff(names(metadata_files), names(expression_files))
if (length(orphan_expression)) {
  message("Skipping ", length(orphan_expression), " expression file(s) with no matching metadata: ",
          paste(orphan_expression, collapse = ", "))
}
if (length(orphan_metadata)) {
  message("Skipping ", length(orphan_metadata), " metadata file(s) with no matching expression: ",
          paste(orphan_metadata, collapse = ", "))
}
if (length(paired) < 2) {
  stop("Need at least 2 paired datasets to meta-analyse; found ", length(paired),
       ". A meta-analysis of one dataset is just that dataset.", call. = FALSE)
}

# ------------------------------------------------------------ file reading

read_table_auto <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  separator <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
  read.delim(path, sep = separator, check.names = FALSE, stringsAsFactors = FALSE)
}

load_dataset <- function(id) {
  metadata <- read_table_auto(metadata_files[[id]])
  id_column <- opts$`id-column`
  if (is.null(id_column)) {
    id_column <- names(metadata)[1]
  } else if (!id_column %in% names(metadata)) {
    stop("Dataset ", id, ": --id-column '", id_column, "' is not a column in ",
         basename(metadata_files[[id]]), ". Columns are: ",
         paste(names(metadata), collapse = ", "), call. = FALSE)
  }
  if (!opts$`group-column` %in% names(metadata)) {
    stop("Dataset ", id, ": --group-column '", opts$`group-column`, "' is not a column in ",
         basename(metadata_files[[id]]), ". Columns are: ",
         paste(names(metadata), collapse = ", "), call. = FALSE)
  }

  samples <- as.character(metadata[[id_column]])
  metadata <- metadata[, setdiff(names(metadata), id_column), drop = FALSE]
  rownames(metadata) <- samples

  expression <- read_table_auto(expression_files[[id]])
  gene_column <- names(expression)[1]
  genes <- as.character(expression[[gene_column]])
  expression <- expression[, -1, drop = FALSE]

  # Collapse duplicate gene identifiers by their mean.
  matrix_data <- as.matrix(expression)
  storage.mode(matrix_data) <- "numeric"
  if (anyDuplicated(genes)) {
    collapsed <- rowsum(matrix_data, group = genes, reorder = TRUE, na.rm = TRUE)
    counts <- as.vector(table(genes)[rownames(collapsed)])
    matrix_data <- collapsed / counts
  } else {
    rownames(matrix_data) <- genes
  }

  shared <- intersect(colnames(matrix_data), rownames(metadata))
  if (!length(shared)) {
    stop("Dataset ", id, ": no sample name is present in both files.\n",
         "  expression columns start: ", paste(head(colnames(matrix_data), 3), collapse = ", "),
         "\n  metadata IDs start:       ", paste(head(rownames(metadata), 3), collapse = ", "),
         call. = FALSE)
  }
  dropped <- length(union(colnames(matrix_data), rownames(metadata))) - length(shared)
  matrix_data <- matrix_data[, shared, drop = FALSE]
  metadata <- metadata[shared, , drop = FALSE]

  labels <- as.character(metadata[[opts$`group-column`]])
  found <- intersect(c(opts$case, opts$control), unique(labels))
  if (length(found) < 2) {
    stop("Dataset ", id, ": column '", opts$`group-column`, "' does not contain both '",
         opts$case, "' and '", opts$control, "'. It contains: ",
         paste(unique(labels), collapse = ", "), call. = FALSE)
  }

  list(id = id, expression = matrix_data, metadata = metadata,
       n_case = sum(labels == opts$case), n_control = sum(labels == opts$control),
       n_other = sum(!labels %in% c(opts$case, opts$control)), dropped = dropped)
}

message("\nPairing ", length(paired), " dataset(s) from\n  ", opts$expression, "\n  ", opts$metadata, "\n")
datasets <- list()
for (id in paired) datasets[[id]] <- load_dataset(id)

summary_table <- data.frame(
  dataset = vapply(datasets, `[[`, character(1), "id"),
  genes = vapply(datasets, function(d) nrow(d$expression), integer(1)),
  samples = vapply(datasets, function(d) ncol(d$expression), integer(1)),
  case = vapply(datasets, `[[`, integer(1), "n_case"),
  control = vapply(datasets, `[[`, integer(1), "n_control"),
  other = vapply(datasets, `[[`, integer(1), "n_other"),
  unmatched = vapply(datasets, `[[`, numeric(1), "dropped"),
  row.names = NULL
)
print(summary_table, row.names = FALSE)

shared_genes <- Reduce(intersect, lapply(datasets, function(d) rownames(d$expression)))
message("\n", length(shared_genes), " gene identifiers are present in all ",
        length(datasets), " datasets.")
if (prop_dataset == 1 && length(shared_genes) == 0) {
  stop("No gene is shared by every dataset, but --prop-dataset is 1, so nothing ",
       "could pass. Either harmonise your gene identifiers or lower --prop-dataset.",
       call. = FALSE)
}
if (sum(summary_table$other) > 0) {
  message("Note: ", sum(summary_table$other), " sample(s) carry a label other than '",
          opts$case, "' or '", opts$control, "'. DExMA ignores them.")
}

if (isTRUE(opts$`dry-run`)) {
  message("\n--dry-run: inputs are valid. Remove the flag to run the meta-analysis.")
  quit(status = 0)
}

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------- per-assay handling

# The two branches below are the normalisation used in
# 01_meta_analysis/DExMA_meta_microarray.R and DExMA_meta_rnaseq.R
# respectively, unchanged apart from being parameterised.
build_eset <- function(dataset) {
  expression <- dataset$expression
  metadata <- dataset$metadata

  if (opts$type == "rnaseq") {
    suppressPackageStartupMessages({ library(edgeR); library(limma) })
    if (any(expression %% 1 != 0, na.rm = TRUE)) {
      warning("Dataset ", dataset$id, " holds non-integer values. --type rnaseq expects ",
              "raw counts; if these are already normalised, use --type microarray instead.",
              call. = FALSE, immediate. = TRUE)
    }
    dge <- DGEList(counts = expression)
    groups <- factor(metadata[[opts$`group-column`]])
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    expression <- voom(dge, design, plot = FALSE)$E
  } else {
    suppressPackageStartupMessages(library(impute))
    if (anyNA(expression)) {
      expression <- impute.knn(expression, k = impute_k)$data
    }
  }

  ExpressionSet(assayData = expression,
                phenoData = AnnotatedDataFrame(metadata))
}

message("\nNormalising (", opts$type, ") ...")
meta_object <- list()
failed <- character(0)
for (id in names(datasets)) {
  result <- tryCatch({
    eset <- build_eset(datasets[[id]])
    elementObjectMA(expressionMatrix = eset,
                    groupPheno = opts$`group-column`,
                    expGroup = opts$case,
                    refGroup = opts$control)
  }, error = function(e) {
    message("  ", id, ": FAILED - ", conditionMessage(e))
    NULL
  })
  if (is.null(result)) {
    failed <- c(failed, id)
  } else {
    meta_object[[id]] <- result
    message("  ", id, ": ok")
  }
}
if (length(failed)) {
  message("\n", length(failed), " dataset(s) failed and are excluded: ",
          paste(failed, collapse = ", "))
}
if (length(meta_object) < 2) {
  stop("Only ", length(meta_object), " dataset(s) survived normalisation; ",
       "nothing to meta-analyse.", call. = FALSE)
}

if (!isTRUE(opts$`no-id-harmonisation`)) {
  message("\nHarmonising identifiers to ", opts$`gene-id`, " (", opts$organism, ") ...")
  meta_object <- allSameID(meta_object, finalID = opts$`gene-id`, organism = opts$organism)
}
meta_object <- dataLog(meta_object)

message("\nHeterogeneity:")
print(heterogeneityTest(meta_object))

# ------------------------------------------------------------ meta-analysis

message("\nRandom-effects meta-analysis over ", length(meta_object), " datasets ...")
results <- as.data.frame(metaAnalysisDE(meta_object, typeMethod = "REM"))
results$symbol <- rownames(results)

write.csv(results, file.path(opts$out, "meta_analysis_all_genes.csv"), row.names = FALSE)

filtered <- results[
  results$propDataset >= prop_dataset &
    abs(results$Zval) > min_zval &
    results$FDR < max_fdr &
    abs(results$Com.ES) > min_effect, , drop = FALSE]

if (isTRUE(opts$`protein-coding`)) {
  if (!requireNamespace("annotables", quietly = TRUE)) {
    stop("--protein-coding needs the annotables package:\n",
         "  remotes::install_github(\"stephenturner/annotables\")", call. = FALSE)
  }
  coding <- unique(annotables::grch38$symbol[annotables::grch38$biotype == "protein_coding"])
  before <- nrow(filtered)
  filtered <- filtered[filtered$symbol %in% coding, , drop = FALSE]
  message("Protein-coding filter: ", before, " -> ", nrow(filtered), " genes")
}

filtered <- filtered[!duplicated(filtered$symbol), , drop = FALSE]
up <- filtered[filtered$Com.ES > min_effect, , drop = FALSE]
down <- filtered[filtered$Com.ES < -min_effect, , drop = FALSE]
up <- up[order(-up$Com.ES), , drop = FALSE]
down <- down[order(down$Com.ES), , drop = FALSE]

write.csv(up, file.path(opts$out, "meta_genes_up.csv"), row.names = FALSE)
write.csv(down, file.path(opts$out, "meta_genes_down.csv"), row.names = FALSE)
write.table(up$symbol, file.path(opts$out, "meta_genes_up.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(down$symbol, file.path(opts$out, "meta_genes_down.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# A record of exactly what produced these files, so a result directory can be
# interpreted months later without the command history.
writeLines(c(
  paste("run date:        ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("assay type:      ", opts$type),
  paste("datasets:        ", paste(names(meta_object), collapse = ", ")),
  paste("datasets failed: ", if (length(failed)) paste(failed, collapse = ", ") else "none"),
  paste("case / control:  ", opts$case, "/", opts$control),
  paste("group column:    ", opts$`group-column`),
  paste("filters:         propDataset >=", prop_dataset, "| |Z| >", min_zval,
        "| FDR <", max_fdr, "| |Com.ES| >", min_effect),
  paste("protein-coding:  ", isTRUE(opts$`protein-coding`)),
  paste("genes tested:    ", nrow(results)),
  paste("genes passing:   ", nrow(filtered), "(", nrow(up), "up,", nrow(down), "down )"),
  paste("DExMA version:   ", as.character(packageVersion("DExMA"))),
  paste("R version:       ", R.version.string)
), file.path(opts$out, "run_parameters.txt"))

message("\nDone. ", nrow(results), " genes tested; ", nrow(up), " up, ", nrow(down),
        " down passed the filters.")
message("Results in ", normalizePath(opts$out))
message("\nNext step: feed either gene list to the PPI/hub script, e.g.")
message("  python3 02_ppi_network/string_network.py --genes ",
        file.path(opts$out, "meta_genes_up.txt"), " --out ", file.path(opts$out, "ppi"))
