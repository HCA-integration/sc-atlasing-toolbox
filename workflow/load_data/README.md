# Load data

Given a TSV file of datasets with URLs and further information per dataset, do the following:

1. Download data from URL
2. Add metadata
3. Merge if multiple datasets are available per study
4. Filter cell per study
5. Merge data per study to single file per organ or additionally specified subset

For example configuration files, checkout out `test/config.yaml` and `test/datasets.tsv`.
Complete dataset configurations should be available at the top-level pipeline (git root of this repository).

## Testing

### Prepare test data
For the from file loader feature, you must download a test dataset first before you test the complete pipeline.
The following script downloads the SchulteSchrepping dataset and then copies it to the location that is defined in `dataset.tsv`.

```
bash download_test_data.sh -c1
```

This needs to be done only once.

### Run pipeline on test configuration

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n
bash test/run_test.sh -c
```

## Download data

Data is downloaded by the `download` rule.
URLs are either retrieved from the by the CxG collection and dataset IDs or from the input TSV directly.

## Add Metadata

The `metadata` rule adds additional dataset-level information that is included from the input TSV file.
This steps expects the data to follow
the [CELLxGENE schema 3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md)
and is extended as described below.
Other steps include:

+ adding external annotations if annotation file and columns are available in the TSV
+ saving donor IDs under `.obs['donor']`
+ inferring sample ID from input TSV and saving it under `.obs['sample']`

AnnData object must contain:

+ `.X` raw counts, sparse format
+ `.uns['meta']` metadata from TSV file
+ `.obs` columns from CELLxGENE schema 3.0.0 + a subset of information in `.uns['meta']` from `EXTRA_COLUMNS`
  from `scripts/utils.py`
+ `.obs['dataset']`, `.uns['dataset']` name of task/dataset
+ `.obs['organ']`, `.uns['organ']` organ
+ `.obs['donor']` donor ID
+ `.obs['sample']` sample ID (inferred from input TSV)
+ `.obs['barcode']` cell barcodes as declared in index
+ `.obs['author_annotation']` author annotation under the `author_annotation` column of the input TSV
+ `.var` gene information as specified in CELLxGENE schema 3.0.0
+ `.obs.index` unique cell identifiers e.g. dataset + numerical index

The anndata is saved as a zarr file for a better speed to compression tradeoff compared ot gzipped h5ad files.

## Merge Data

This operation is applied by the following rules

+ `merge_study`: merge all datasets of a study if multiple datasets are available, else create a symlink
+ `merge_organ`: merge all studies that belong to an organ
+ `merge_organ_filter`: merge all cells removed by filtering per organ
+ `merge_subset`: merge datasets by subset defined in the input TSV under `subset` (overlapping subsets allowed)

The anndata must contain all the slots described in [Add Metadata](#add-metadata) apart from the `.uns` slot.

## Filter

Filter cells per study depending on the `config.yaml` specification.
Two keys are available for controlling the filtering behaviour, `filter_per_organ` specifies global filter paramters
for all datasets per organ and `filter_per_study` allows for study specific filter options.
An example of an organ-level filter specification is shown below:

```yaml
filter_per_organ:
  blood:
    cells_per_sample:
      min: 50
      max: 10000
    mito_pct: 30
    remove_by_colum:
      dataset:
        - Lee2020_2
```

All organ-level filtering decisions are applied per study.
The `remove_by_column` key can include any columns that are available in the anndata objects per study.

An example of per study filters shows that only the studies that require further filtering need to be overwritten.
The filter options are the same as for the organ-level filters.

```yaml
filter_per_study:
  SchulteSchrepping2020:
    remove_by_colum:
      sample:
        - Schulte-Schrepping_C2P01H_d0
        - Schulte-Schrepping_C2P05F_d0
        - Schulte-Schrepping_C2P07H_d0
        - Schulte-Schrepping_C2P10H_d0
        - Schulte-Schrepping_C2P13F_d0
        - Schulte-Schrepping_C2P15H_d0
        - Schulte-Schrepping_C2P16H_d0
        - Schulte-Schrepping_C2P19H_d0
      donor:
        - C19-CB-0008
      disease:
        - influenza
```

Both keys can be empty or missing from the config file.
In that case, no filtering is applied.
