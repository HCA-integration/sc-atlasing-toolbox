input_file <- snakemake@input$tsv
output_tsv <- snakemake@output$tsv
output_pdf <- snakemake@output$pdf

id_vars <- snakemake@params$id_vars
value_var <- snakemake@params$value_var
variable_var <- snakemake@params$variable_var
weight_batch <- snakemake@params$weight_batch
n_top <- snakemake@params$n_top

tryCatch({
 library(funkyheatmap)
  }, error = function(e) {
    install.packages('funkyheatmap', repos = 'http://cran.us.r-project.org')
  }
)
suppressPackageStartupMessages({
  library(funkyheatmap)
  library(data.table)
  library(tidyverse)
  library(dynutils)
})

dt <- fread(input_file)
#print(head(dt))

# define groups of interest to be plotted in funkyheatmap
integration_setup <- id_vars
bio_metrics <- unique(dt[metric_type == 'bio_conservation', metric])  # c('ari', 'clisi_y', 'asw_label', 'asw_label_y', 'isolated_label_asw', 'isolated_label_f1', 'clisi', 'nmi', 'cell_cycle')
batch_metrics <- unique(dt[metric_type == 'batch_correction', metric])  # c('graph_connectivity', 'asw_batch', 'asw_batch_y', 'ilisi', 'ilisi_y', 'pcr')
metrics <- c(bio_metrics, batch_metrics)

# subset data.table to columns of interest & transform
all_columns <- c(id_vars, variable_var, value_var, 'metric_type')
id_formula <- as.formula(paste(paste(id_vars, collapse='+'), variable_var, sep = '~'))
metrics_tab <- dcast(
  subset(dt, select = all_columns),
  id_formula,
   value.var = value_var
)

# scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
scaled_metrics_tab <- as.matrix(metrics_tab[, metrics, with=FALSE])
scaled_metrics_tab <- as.data.table(apply(scaled_metrics_tab, 2, function(x) scale_minmax(x)))

# calculate average score by group and 
score_group_batch <- rowMeans(scaled_metrics_tab[, batch_metrics, with=FALSE], na.rm = T)
score_group_bio <- rowMeans(scaled_metrics_tab[, bio_metrics, with=FALSE], na.rm = T)

# weighted overall score
score_all <- (weight_batch * score_group_batch + (1 - weight_batch) * score_group_bio)
metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "dataset")
metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .after = "Overall Score")
metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .before = bio_metrics[1])

# order methods by the overall score
metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]

# write output summary.csv file before plotting ? Michaela to check the output dir? 
fwrite(metrics_tab, output_tsv, sep='\t')

#add column info metadata for plotting using funkyheatmap
dt1 <- data.table(id=integration_setup, group="Integration setup", geom='text', palette=NA, hjust=0.5)
dt2 <- data.table(id="Overall Score", geom="bar", palette="Greens")
dt3 <- data.table(id="Bio conservation", group="Bio conservation", geom='bar', palette='Blues')
dt4 <- data.table(id=bio_metrics, group="Bio conservation", geom='funkyrect', palette='Blues')
dt5 <- data.table(id="Batch Correction", group="Batch correction", geom='bar', palette='Reds')
dt6 <- data.table(id=batch_metrics, group="Batch correction", geom='circle', palette='Reds')
column_info <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, fill=TRUE)
print(column_info)

n_top <- min(n_top, nrow(metrics_tab))
g <- funky_heatmap(
  metrics_tab[1:n_top],
  column_info = column_info,
  col_annot_offset = 3.5,
  add_abc = TRUE,
  scale_column = FALSE
)
ggsave(output_pdf, g, device = cairo_pdf, width = g$width, height = g$height, dpi=300)