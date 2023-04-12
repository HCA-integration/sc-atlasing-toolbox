# Label Harmonisation

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

## Input

## Config file

```yaml
DATASETS:
  dataset_name:  # can replace with a representative name
    input:
      label_harmonization: anndata_file_path.h5ad  # can be output from another module
    label_harmonization:
      mapping:  # mapping setup
        file:  test/input/mapping.tsv  # TSV file with cell label mapping
        order:  # order of which label column should be matched to which
          - cell_type  # first entry MUST be an existing column in the anndata object
          - harmonized_label  # first remapping level (mapping 'cell_type' to 'harmonized_label')
          - lineage  # second remapping level (mapping 'harmonized_label' to 'lineage')
```

### Anndata
AnnData file h5ad or zarr file with the following:

+ `.obs` name of task/dataset, must contain the first entry of `config['DATASETS']['dataset_name']['label_harmonization']['mapping']['order']
