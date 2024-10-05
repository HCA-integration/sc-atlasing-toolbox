options(show.error.locations = TRUE)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(GloScope)
  library(reticulate)
  library(anndata)
})
io <- import('utils.io')
sc <- import('scanpy')

input_file <- snakemake@input$zarr
prepare_file <- snakemake@input$prepare
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
adata <- sc$AnnData(
    obs = sample_dt,
    obsm = list(distances = dist_matrix)
)

# compute kNN graph
sc$pp$neighbors(adata, use_rep = 'distances', metric = 'precomputed')

message(paste("Writing to", output_file, "..."))
print(adata)
io$write_zarr_linked(
    adata,
    in_dir=prepare_file,
    out_dir=output_file,
    files_to_keep=c('obsm', 'obsp', 'uns'),
)