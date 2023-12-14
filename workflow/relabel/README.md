# Label Harmonisation

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

## Input

### Config file

```yaml
DATASETS:
  dataset_name:  # can replace with a representative name
    input:
      relabel:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    relabel:
      new_columns:  # mapping setup for new columns
        file:  test/input/mapping.tsv  # TSV file with cell label mapping
        order:  # order of which label column should be matched to which
          - cell_type  # first entry MUST be an existing column in the anndata object
          - harmonized_label  # first remapping level (mapping 'cell_type' to 'harmonized_label')
          - lineage  # second remapping level (mapping 'harmonized_label' to 'lineage')
      merge_columns:  # mapping setup for merging existing columns
        file:  test/input/merge_test.tsv  # TSV file with cell label mapping
        sep: '-'  # separator for merging columns
```

### Mapping files

#### New columns

The new column mapping file must be in TSV format with at least one column in the anndata object.
The columns must include all columns that are specified in the `order` list.

Example for `test/input/pbmc68k.h5ad`:

```
bulk_labels	lineage
CD14+ Monocyte	myeloid
Dendritic	lymphoid
CD56+ NK	lymphoid
CD4+/CD25 T Reg	lymphoid
```

#### Merge existing columns

The merge column mapping file must be in TSV format with the following columns:

* `file_id`: input file id that is used in the pipeline to match which file to add the new column to
* `column_name`: new column name
* `columns`: columns to merge, separated by `,` (no whitespaces)

Example for `test/input/pbmc68k.h5ad`:

```
file_id	column_name	columns
file_1	cc_joint	S_score,G2M_score
file_1	ct_joint	bulk_labels,louvain
file_2	counts_joint	n_genes,n_counts
```

### Anndata
AnnData file h5ad or zarr file with the following:

+ `.obs` with all columns specified in `config['DATASETS'][dataset]['relabel']['new_columns']['order']`