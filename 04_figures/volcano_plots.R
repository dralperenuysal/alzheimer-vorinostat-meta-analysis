# Figure 1 - volcano plots of the two meta-analyses.
#
# Both input files are already the filtered DEG sets (detected in every
# dataset, |Z| > 2, FDR < 0.05), so every point plotted is significant; the
# colour carries the direction, given by the sign of the combined effect size.
# The dashed vertical line at zero is that boundary, not a threshold.

library(tidyverse)
library(patchwork)

OUT <- "results/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MM <- 1 / 25.4
FULL_W <- 174 * MM          # Springer double-column width

read_pair <- function(up, down) {
  bind_rows(read.csv(up), read.csv(down)) %>%
    as_tibble() %>%
    mutate(minus_log10FDR = -log10(FDR),
           direction = if_else(Com.ES > 0, "Up", "Down"))
}

volcano <- function(df, title) {
  ggplot(df, aes(x = Com.ES, y = minus_log10FDR)) +
    geom_point(aes(color = direction, size = minus_log10FDR),
               alpha = 0.55, stroke = 0) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    scale_color_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
    scale_size_continuous(range = c(0.15, 1.5)) +
    labs(x = "Meta effect size (Com.ES)",
         y = expression(-log[10]("FDR")),
         color = NULL, size = expression(-log[10]("FDR")),
         title = title) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "right",
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(size = 8))
}

micro <- read_pair("results/meta_analysis/microarray_up.csv",
                   "results/meta_analysis/microarray_down.csv")
rna <- read_pair("results/meta_analysis/rnaseq_up.csv",
                 "results/meta_analysis/rnaseq_down.csv")

fig <- volcano(micro, "Microarray") + volcano(rna, "RNA-seq") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(file.path(OUT, "Fig1.pdf"), fig,
       width = FULL_W, height = FULL_W * 0.42, units = "in", device = cairo_pdf)

cat(sprintf("wrote %s/Fig1.pdf  (microarray %d, RNA-seq %d points)\n",
            OUT, nrow(micro), nrow(rna)))
