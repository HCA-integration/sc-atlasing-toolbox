# Workflow

Given a TSV file and a schema mapping, the pipeline does the following:

1. [Load data from CELLxGENE, DCP, URL or file](#load-data)
2. [Harmonize metadata](#harmonize-metadata)
3. [Merge datasets per study, organ or custom subset](#merge-data)
4. [Filter cells](#filter)

For example configuration files, check out `test/config.yaml` and `test/datasets.tsv`.
Complete dataset configurations should be available under `configs` at the top-level pipeline (git root of this repository).


## Load data

Files can either be read from a specified input file or downloaded from an URL, CELLxGENE or DCP directly.
Which way a dataset is loaded, depends on the [dataset mapping](#dataset-file-mapping).

## Harmonize Metadata

The `metadata` rule adds additional dataset-level information that is included from the input TSV file.
This steps expects the data to follow
the [CELLxGENE schema 3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md)
and is extended as described below.
Other steps include:

+ adding external annotations if annotation file and columns are available in the TSV
+ saving donor IDs under `.obs['donor']`
+ inferring sample ID from input TSV and saving it under `.obs['sample']`

The input AnnData objects must contain:

+ `.X` raw counts, sparse format
+ `.obs` columns from the schema defined in the schema 
+ `.var` gene information as specified in CELLxGENE schema 3.0.0
+ `.obs.index` cell barcode

The output AnnData objects will contain:
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

The `AnnData` is saved as a zarr file for a better speed to compression tradeoff compared ot gzipped h5ad files.

## Merge Data

This operation is applied by the following rules

+ `merge_study`: merge all datasets of a study if multiple datasets are available, else create a symlink
+ `merge_organ`: merge all studies that belong to an organ
+ `merge_organ_filter`: merge all cells removed by filtering per organ
+ `merge_subset`: merge datasets by subset defined in the input TSV under `subset` (overlapping subsets allowed)

The `AnnData` must contain all the slots described in [Harmonize Metadata](#harmonize-metadata) apart from the `.uns` slot.

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

# Preparing the input data

In order to specify your input, you need to define the following files with the file locations and dataset-level metadata. Additionally, you need to prepare your input files to contain certain metadata.

* Dataset file mapping (`configs/datasets.tsv`)
* Schema mapping (e.g. `configs/schema_mapping.tsv`)
* Configuration file (e.g. `configs/imorted/config.yaml`)
* DCP metadata (optional)

## Dataset file mapping

The dataset mapping is the main configuration for which datasets you want to include together with dataset-level configuration.
This is the place to define which datasets (`AnnData` objects) you want to define in your atlas building workflow.

| Column            | Description                                                                                                                                                                                                         |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| dataset           | Name of the dataset, multiple datasets make up a study                                                                                                                                                              |
| study             | Name of the study, used for aggregating datasets to study level                                                                                                                                                     |
| organ             | Name of the organ, can be tissue or any other name for aggregating the atlas                                                                                                                                        |
| donor_column      | Column with donor IDs                                                                                                                                                                                               |
| sample_column     | Column with sample IDs                                                                                                                                                                                              |
| author_annotation | Column with author annotations (needed for annotation quality assessment and label harmonisation)                                                                                                                   |
| cell_type         | Column with cell ontology labels (needed to for different versions of the CELLxGENE schema, TODO: deprecate)                                                                                                        |
| schema            | Name of schema to be mapped to `cellxgene`. Naming must match the columns in schema mapping.                                                                                                                        |
| url               | URL or file name of the `h5ad` file. For data from CELLxGENE or DCP, the `cellxgene` or `dcp` are sufficient, if `collection_id` and `dataset_id` are defined for `cellxgene` or `project_uuid` is defined for `dcp |
| collection_id     | Only for `url=cellxgene`, if data should be downloaded directly from CELLxGENE.                                                                                                                                     |
| dataset_id        | Only for `url=cellxgene`, if data should be downloaded directly from CELLxGENE.                                                                                                                                     |
| project_uuid      | Only for `url=dcp`, if data should be downloaded directly from HCA DCP data portal                                                                                                                                  |
| annotation_file   | Optional. Any additional annotations that are not in the `AnnData` object. Needs to have a matching barcode column.                                                                                                 |
| barcode_column    | Optional. Column in `AnnData.obs` for merging external annotations                                                                                                                                                  |

All other columns are optional and will be added to `AnnData.uns['meta']`.

## Schema Mapping

The data loader ensures that the data adheres to the [CELLxGENE schema 3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md) specifications.
For datasets that do not adhere to that schema, the schema mapping file allows to provide a mapping of custom `AnnData.obs` columns to the ones defined in CELLxGENE.
Below is an example of a schema mapping for the schemas `custom` and `dcp` to `cellxgene`.
The column names are used to map the dataset to its corresponding schema.

```
cellxgene                custom                       dcp                                                            
study                    study_PI                     project.contributors.name                                      
sample                   sample_ID                    specimen_from_organism.biomaterial_core.biomaterial_id         
donor_id                 subject_ID                   donor_organism.biomaterial_core.biomaterial_id                 
development_stage        subject_developmental_state  donor_organism.development_stage.text                          
sex                      sex                          donor_organism.sex                                             
self_reported_ethnicity  ethnicity_free_text          donor_organism.human_specific.ethnicity.text                   
suspension_type          biological_unit              library_preparation_protocol.nucleic_acid_source               
assay                    library_platform             library_preparation_protocol.library_construction_method.text  
organism                 species                      donor_organism.genus_species.text                              
cell_type                cell_type                                                                                   
...                      ...                          ...                                                            
```


## DCP metadata (optional)

This file is optional and used for datasets for which users want to map additional DCP annotations.
The mapping should contain a `study` and a `filename` column, where `filename` is a TSV file that follows the [DCP metadata schema](https://data.humancellatlas.org/metadata).

The mapping does not have to include all the studies that the user wants to include.
If a study that is defined in the dataset mapping, but not in the DCP mapping, this just affects rules that require the DCP mapping.

TODO: extend to other metadata input.

## Configuration file

This file is the heart of the pipeline.
Snakemake can use multiple configuration files in one run and will concatinate files on order of definition.
For data loading, you just need to define the location of the files defined above.
By default, the `configs/config.yaml` file should already contain the correct paths, but if you configured your datasets with different files, you need to update them accordingly.

```yaml
dataset_meta: configs/datasets.tsv
schema_file: data/input/schema_mapping.tsv
dcp_metadata: configs/dcp_metadata.tsv

filter_per_organ:
...

filter_per_study:
...
```

# Testing

## Prepare test data
For the from file loader feature, you must download a test dataset first before you test the complete pipeline.
The following script downloads the SchulteSchrepping dataset and then copies it to the location that is defined in `dataset.tsv`.

```
bash download_test_data.sh -c1
```

This needs to be done only once.

## Run pipeline on test configuration

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n
bash test/run_test.sh -c
```
