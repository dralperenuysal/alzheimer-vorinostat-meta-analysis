# Figure 2 - overlap between the microarray and RNA-seq DEG lists.
#
# Two two-set diagrams, one per direction. The intersections drawn here are the
# same ones result_intersection.R writes and the PPI step consumes.

library(tidyverse)
library(patchwork)

OUT <- "results/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MM <- 1 / 25.4
FULL_W <- 174 * MM

read_genes <- function(path) unique(trimws(readLines(path)))

circle <- function(cx, r = 1, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(t), y = r * sin(t))
}

venn <- function(rna_file, micro_file, title) {
  rna <- read_genes(rna_file)
  micro <- read_genes(micro_file)
  both <- intersect(rna, micro)

  shapes <- bind_rows(
    circle(-0.55) %>% mutate(set = "RNA-seq"),
    circle(0.55) %>% mutate(set = "Microarray")
  )
  labels <- tibble(
    x = c(-1.05, 0, 1.05), y = 0,
    label = c(length(rna) - length(both), length(both),
              length(micro) - length(both))
  )
  headers <- tibble(
    x = c(-0.85, 0.85), y = 1.35,
    label = c("RNA-seq", "Microarray"),
    set = c("RNA-seq", "Microarray")
  )

  ggplot() +
    geom_polygon(data = shapes, aes(x, y, fill = set, colour = set),
                 alpha = 0.28, linewidth = 0.3) +
    geom_text(data = labels, aes(x, y, label = label), size = 2.8) +
    geom_text(data = headers, aes(x, y, label = label, colour = set),
              size = 2.6, fontface = "bold") +
    scale_fill_manual(values = c("RNA-seq" = "#C2185B", "Microarray" = "#3949AB")) +
    scale_colour_manual(values = c("RNA-seq" = "#C2185B", "Microarray" = "#3949AB")) +
    coord_equal(clip = "off") +
    labs(title = title) +
    theme_void(base_size = 8) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4))
}

a <- venn("results/meta_analysis/rnaseq_meta_genes_up.txt",
          "results/meta_analysis/microarray_meta_genes_up.txt", "Upregulated")
b <- venn("results/meta_analysis/rnaseq_meta_genes_down.txt",
          "results/meta_analysis/microarray_meta_genes_down.txt", "Downregulated")

fig <- a + b +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(file.path(OUT, "Fig2.pdf"), fig,
       width = FULL_W, height = FULL_W * 0.38, units = "in", device = cairo_pdf)

cat(sprintf("wrote %s/Fig2.pdf\n", OUT))
