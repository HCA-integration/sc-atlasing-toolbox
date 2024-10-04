options(show.error.locations = TRUE)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(GloScope)
  library(reticulate)
  library(anndata)
})
io <- import('utils.io')
ad <- import('anndata')

input_file <- snakemake@input$zarr
output_file <- snakemake@output$zarr

use_rep <- snakemake@params$use_rep
sample_key <- snakemake@params$sample_key
k <- snakemake@params$k
seed <- snakemake@params$seed
threads <- snakemake@threads


message("Read data...")
adata <- io$read_anndata(
    input_file,
    obs='obs',
    obsm='obsm',
)

embedding <- adata$obsm[[use_rep]]
rownames(embedding) <- adata$obs_names
sample_ids <- adata$obs[[sample_key]] # TOOD: is this correct?
obs_dt <- data.table(adata$obs)
sample_dt <- obs_dt[, lapply(.SD, first), by = sample_key]

message(paste(c("Data dimensions:", paste(dim(embedding), collapse=', '))))
message(paste("Number of IDs:", uniqueN(sample_ids)))

message("Calculating distance matrix...")
dist_matrix <- gloscope(
    embedding,
    sample_ids,
    dens = "GMM",
    dist_mat = "KL",
    k = k,
    BPPARAM = BiocParallel::MulticoreParam(workers = 4, RNGseed = seed)
)
message(paste(c("Distance matrix dimensions:", paste(dim(dist_matrix), collapse=', '))))

# create new anndata object
message("Creating new anndata object...")
adata <- ad$AnnData(
    X = dist_matrix,
    obs = sample_dt,
    var = sample_dt,
)
print(adata)

message(paste("Writing to", output_file, "..."))
io$write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=c('X', 'layers' ,'obs', 'var', 'obsm', 'varm', 'obsp', 'varp', 'uns'),
)