# Split Data

This module is used to split cells into multiple files depending on a categorical value of a column in ``.obs``.

## Input

Example config:

```yaml
output_dir: test/out
images: test/images

DATASETS:
  test:
    input:
      split_data:
        pbmc: test/input/pbmc68k.h5ad
    split_data:
      key: bulk_labels
      values:
        - CD4+_CD45RA+_CD25-_Naive_T
        - Dendritic
        - CD14+_Monocyte
        - CD19+_B
```

Example output:

```shell
test/out/split_data
├── dataset~test
│   └── file_id~pbmc
│       └── key~bulk_labels
│           ├── value~CD14+_Monocyte.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD14+_Monocyte.zarr
│           ├── value~CD19+_B.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD19+_B.zarr
│           ├── value~CD4+_CD45RA+_CD25-_Naive_T.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~CD4+_CD45RA+_CD25-_Naive_T.zarr
│           └── value~Dendritic.zarr -> ../../../splits/dataset~test/file_id~pbmc/key~bulk_labels/value~Dendritic.zarr
├── input_files.tsv
└── splits
    └── dataset~test
        └── file_id~pbmc
            └── key~bulk_labels
```