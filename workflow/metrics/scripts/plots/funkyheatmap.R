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
    repos = snakemake@params$cran_url
    if(is.null(repos)){
      repos = 'https://cloud.r-project.org'
    }
    install.packages('funkyheatmap', repos = repos)
  }
)
suppressPackageStartupMessages({
  library(funkyheatmap)
  library(data.table)
  library(tidyverse)
  library(dynutils)
})

# read files
dt <- fread(input_file)
dt[, (id_vars) := lapply(.SD, as.character), .SDcols = id_vars]
print(head(dt))

extra_columns <- readLines(snakemake@input$extra_columns)
print(extra_columns)
id_vars <- unique(c(id_vars, extra_columns))

# remove unintegrated output types without corresponding method
if ('unintegrated' %in% dt$method & uniqueN(dt$method) > 1) {
  ot_count <- dt[method != 'unintegrated', .(output_type, method)] %>%
    unique %>%
    .[, output_type] %>%
    table
  keep_ot <- ot_count[ot_count > 1]
  df <- dt[, output_type %in% keep_ot]
}

# define groups of interest to be plotted in funkyheatmap
integration_setup <- id_vars
bio_metrics <- unique(dt[metric_type == 'bio_conservation' & !is.na(score), metric]) 
batch_metrics <- unique(dt[metric_type == 'batch_correction' & !is.na(score), metric])
metrics <- c(bio_metrics, batch_metrics)

# subset data.table to columns of interest & transform
all_columns <- c(id_vars, variable_var, value_var, 'metric_type')
id_formula <- paste(paste(id_vars, collapse='+'), variable_var, sep = '~')
cat(paste(id_formula, '\n'))
metrics_tab <- dcast(
  subset(dt, select = all_columns),
  as.formula(id_formula),
  value.var = value_var,
  fun.aggregate = mean
)
print(metrics_tab)

# scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
scaled_metrics_tab <- as.matrix(metrics_tab[, metrics, with=FALSE])
if (nrow(scaled_metrics_tab) > 1) {
  scaled_metrics_tab <- as.data.table(apply(scaled_metrics_tab, 2, function(x) scale_minmax(x)))
} else {
  scaled_metrics_tab <- as.data.table(scaled_metrics_tab)
}
cat('Scaled metrics table: \n')
print(scaled_metrics_tab)

# calculate average score by group and 
score_group_batch <- rowMeans(scaled_metrics_tab[, batch_metrics, with=FALSE], na.rm = TRUE)
score_group_bio <- rowMeans(scaled_metrics_tab[, bio_metrics, with=FALSE], na.rm = TRUE)

# weighted overall score
score_all <- (weight_batch * score_group_batch + (1 - weight_batch) * score_group_bio)
if(length(score_group_batch) > 0){
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .before = batch_metrics[1])
}
if(length(score_group_bio) > 0){
  metrics_tab <- add_column(metrics_tab, "Bio Conservation" = score_group_bio, .before = bio_metrics[1])
}
if(length(score_all) > 0){
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .before = "Batch Correction")
  metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
}
# order methods by the overall score

# write output summary.csv file before plotting ? Michaela to check the output dir?
cat('Writing output file: ', output_tsv, '\n')
fwrite(metrics_tab, output_tsv, sep='\t')

#add column info metadata for plotting using funkyheatmap
dt1 <- data.table(id=integration_setup, group="Integration setup", geom='text', palette=NA, hjust=0.5)
dt2 <- data.table(id="Overall Score", geom="bar", palette="Greens")
dt3 <- data.table(id="Bio conservation", group="Bio conservation", geom='bar', palette='Blues')
dt4 <- data.table(id=bio_metrics, group="Bio conservation", geom='funkyrect', palette='Blues')
dt5 <- data.table(id="Batch Correction", group="Batch correction", geom='bar', palette='Reds')
dt6 <- data.table(id=batch_metrics, group="Batch correction", geom='circle', palette='Reds')
column_info <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, fill=TRUE)
column_info = column_info[column_info$id %in% colnames(metrics_tab)]
print("column_info")
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
