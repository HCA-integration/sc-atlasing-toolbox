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
