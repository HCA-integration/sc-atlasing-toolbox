library(data.table)
library(ggplot2)
library(ComplexUpset)

min_size <- snakemake@params$min_size

dt <- fread(snakemake@input$tsv)
dt[, value := TRUE]
dt_cast <- dcast(dt, index ~ sample, value.var = 'value', fill=FALSE)[, .SD, .SDcols = !'index']

n_samples <- uniqueN(dt$sample)
min_size <- min(min_size, uniqueN(dt$index))

value_counts <- dt[, uniqueN(sample), by = index]
n_cols = length(value_counts[V1 >= min_size, V1])

upset(
    dt_cast,
    colnames(dt_cast),
    name = 'samples',
    height_ratio = 5,
    width_ratio = 0.2,
    min_size = min_size,
    base_annotations=list('Size'=(intersection_size(counts=FALSE)))
)
ggsave(snakemake@output$png, width = 5 + 0.09 * n_cols, height = 2.5 + 0.2 * n_samples, dpi = 300)
