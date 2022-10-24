# Load data

Given a TSV file of datasets with URLs and further information per dataset, do the following:

1. Download data from URL & add metadata
2. Merge datasets to single h5ad per organ

For example configuration files, checkout out `test/config.yaml` and `test/datasets.tsv`.
Complete dataset configurations should be available at the top-level pipeline (git root of this repository).

## Download data

AnnData object must contain:

+ `.X` raw counts
+ `.uns['dataset']` name of task/dataset
+ `.obs['dataset']`
+ `.obs['organ']`
+ `.uns['organ']`
+ `.uns['meta']` metadata related to data download
  + `.uns['meta']['organ']`
  + `.uns['meta']['dataset_name']`
  + `.uns['meta']['dataset_id']`
  + `.uns['meta']['url']`
  + `.uns['meta']['modalities']`
  + `.uns['meta']['n_samples']`
  + `.uns['meta']['n_cells']`
  + `.uns['meta']['health_status']`
  + `.uns['meta']['demographic']`
  + `.uns['meta']['10x_assay']`
  + `.uns['meta']['sample_column']`
  + `.uns['meta']['donor_column']`


## Merged data

+ `.X` raw counts
+ `.uns['dataset']` name of organ
+ `.obs['organ']`
+ `.obs['dataset']`
+ `.obs['donor']`
+ `.obs['sample']`
+ `.obs['label']`


## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n
bash test/run_test.sh -c
```
