import requests

# CELLxGENE columns of schema 3.0.0
# https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#general-requirements

CELLxGENE_OBS = [
    'assay',
    'assay_ontology_term_id',
    'cell_type',
    'cell_type_ontology_term_id',
    'development_stage',
    'development_stage_ontology_term_id',
    'disease',
    'disease_ontology_term_id',
    'donor_id',
    'is_primary_data',
    'organism',
    'organism_ontology_term_id',
    'self_reported_ethnicity',
    'self_reported_ethnicity_ontology_term_id',
    'sex',
    'sex_ontology_term_id',
    'suspension_type',
    'tissue',
    'tissue_ontology_term_id',
]

CELLxGENE_VARS = [
    'feature_name',
    'feature_reference',
    'feature_biotype'
]

EXTRA_COLUMNS = [
    'barcode',
    'organ',
    'donor',
    'sample',
    'author_annotation',
    'reference',
    'study',
    'dataset',
    'modalities',
    'pipeline_version',
    'institution',
]


def get_union(*args: list):
    """

    :param args: lists of values
    :return: list containing union of all list values
    """
    for arg in args:
        assert isinstance(arg, list)

    sets = [set(arg) for arg in args]
    if len(sets) <= 1:
        return list(sets[0])

    union_set = sets[0]
    for _set in sets[1:]:
        union_set = union_set.union(_set)

    return list(union_set)


def get_from_dataset(dataset_df, key, value, column=None, debug=False):
    """
    Retrieve values in column from dataset_df by subset (key, value)

    :param dataset_df:
    :param key: column in dataset_df to subset to
    :param value: value to subset to
    :param column: column in dataset_df
    :param debug:
    :return:
    """
    if column is None:
        column = dataset_df.columns
    sub = dataset_df.query(f'{key} == @value')
    if debug:
        print(f'##DEBUG## {key, value, column}')
        print(sub)
        print(sub[column])
    return sub[column]


def get_url(dataset_df, wildcards):
    url = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'url').iloc[0]
    if not isinstance(url, str):
        # infer CxG URL from collection and dataset IDs
        collection_id = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'collection_id').iloc[0]
        dataset_id = get_from_dataset(dataset_df, 'dataset', wildcards.dataset, 'dataset_id').iloc[0]
        url = f'https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}/datasets/{dataset_id}/assets/'
        assets = requests.get(url=url).json()
        asset = [a for a in assets if a["filetype"] == "H5AD"][0]
        url = asset["presigned_url"]
    return url


def unlist_dict(dictionary):
    return {
        k: v[0] if isinstance(v, list) and len(v) == 1 else v
        for k, v in dictionary.items()
    }
