library(DExMA)
library(Biobase)
library(glue)
library(limma)
library(impute)
library(dplyr)
library(annotables)
library(tibble)
library(ggplot2)

suitable <- read.table("results/dataset_selection/microarray_check.tsv", header = T) %>% filter(status == "suitable")
series <- suitable$accession

# Initialize an empty list to store ExpressionSet objects
esets <- list()
pheno.list <- list()

# Loop through each series in the 'series' vector
for (serie in series) {
  expression_file <- glue("data/microarray_processed/{serie}/{serie}_processed_expression.tsv")
  pheno_file <- glue("data/microarray_processed/{serie}/{serie}_processed_sampleinfo.tsv")
  
  # Check if both files exist
  if (file.exists(expression_file) && file.exists(pheno_file)) {
    # Use tryCatch to handle potential errors
    tryCatch({
      # Read expression data. The gene identifiers live in a column named
      # gene, the same as the RNA-seq branch expects. They may be symbols or
      # Ensembl IDs; allSameID() maps them onto GeneSymbol below.
      exprs <- read.delim(expression_file, sep = '\t', check.names = F) %>% distinct(.keep_all = T)
      if (!"gene" %in% colnames(exprs)) {
        stop(glue("{expression_file}: no gene column. The first column must be headed gene."))
      }
      exprs <- column_to_rownames(exprs, var = "gene")
      exprs <- as.matrix(exprs)  # Convert to matrix

      # Read phenotype data
      phenodata <- read.delim(pheno_file, sep = '\t', check.names = F, row.names = "samples")
      pheno.list[[serie]] <- phenodata
      phenodata.adf <- AnnotatedDataFrame(phenodata)  # Convert to AnnotatedDataFrame
      
      # Impute the data
      exprs.imp <- impute.knn(exprs, k = 7)
      exprs <- exprs.imp$data
      
      # Create ExpressionSet and add to list
      esets[[serie]] <- ExpressionSet(assayData = exprs,
                                      phenoData = phenodata.adf)
    }, error = function(e) {
      # Handle the error: print a message and continue
      message(glue("Error processing series {serie}: {e$message}"))
    })
  } else {
    warning(glue("Files for series {serie} do not exist: {expression_file}, {pheno_file}"))
  }
}

intersection <- intersect(names(pheno.list), names(esets))
pheno.list <- pheno.list[intersection]
meta.obj <- list()

# Loop through each eset in esets using names
for (name in names(esets)) {
  eset <- esets[[name]]
  
  meta.element <- elementObjectMA(expressionMatrix = eset,
                                  groupPheno = "state",
                                  expGroup = "AD",
                                  refGroup = "control")
  #imputation <- missGenesImput(meta.element, k=7)
  #meta.element.imput <- imputation$meta.element
  meta.obj[[name]] <- meta.element
}

meta.obj <- allSameID(meta.obj, finalID="GeneSymbol", organism = "Homo sapiens")
meta.obj <- dataLog(meta.obj)

heterogeneityTest(meta.obj)

resultsMA <- metaAnalysisDE(meta.obj, typeMethod="REM")
results.df <- as.data.frame(resultsMA)
results.df.filtered <- results.df %>%
  filter(propDataset == 1) %>%
  filter(abs(Zval) > 2) %>%
  filter(FDR < 0.05)

# Filter the protein_coding genes
annot <- grch38 %>% dplyr::select(c(3,4,5,6,8,9)) %>%
  filter(biotype == "protein_coding")

results.df.filtered$symbol <- rownames(results.df.filtered)
results.df.filtered <- results.df.filtered %>%
  filter(symbol %in% annot$symbol)
results.merged <- merge(results.df.filtered, annot)
results.merged <- results.merged %>% distinct(symbol, .keep_all = T)
results.merged.up <- results.merged %>%
  filter(Com.ES > 0) %>%
  arrange(desc(Com.ES))
results.merged.down <- results.merged %>%
  filter(Com.ES < 0) %>%
  arrange(Com.ES)

write.csv(results.merged.up, "results/meta_analysis/microarray_up.csv")
write.csv(results.merged.down, "results/meta_analysis/microarray_down.csv")

write.table(results.merged.up$symbol, file = "results/meta_analysis/microarray_meta_genes_up.txt", quote = F, col.names = F, row.names = F)
write.table(results.merged.down$symbol, file = "results/meta_analysis/microarray_meta_genes_down.txt", quote = F, col.names = F, row.names = F)

makeHeatmap(objectMA = meta.obj, resMA = results.df.filtered, scaling = "zscor",
            regulation = "up", numSig=30,
            fdrSig = 0.05, logFCSig = 0.5, show_rownames = TRUE)

makeHeatmap(objectMA = meta.obj, resMA = results.df, scaling = "zscor",
            regulation = "down", numSig=10,
            fdrSig = 0.05, logFCSig = 0.5, show_rownames = TRUE)

