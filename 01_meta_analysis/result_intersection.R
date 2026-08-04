library(dplyr)

# upregulation
rnaseq <- read.delim("results/meta_analysis/rnaseq_meta_genes_up.txt", header = F)
microarray <- read.delim("results/meta_analysis/microarray_meta_genes_up.txt", header = F)

intersection <- intersect(rnaseq$V1, microarray$V1)

write.table(intersection, file = "results/meta_analysis/intersected_genes_up_alz.txt", row.names = F, col.names = F, quote = F)

# downregulation
rnaseq_down <- read.delim("results/meta_analysis/rnaseq_meta_genes_down.txt", header = F)
micro_down <- read.delim("results/meta_analysis/microarray_meta_genes_down.txt", header = F)

intersection.alz.down <- intersect(rnaseq_down$V1, micro_down$V1)
write.table(intersection.alz.down, file = "results/meta_analysis/intersected_genes_down_alz.txt", row.names = F, col.names = F, quote = F)

