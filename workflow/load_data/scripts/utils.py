CELLxGENE_COLUMNS = [
    'organism',
    #'donor_id',
    'development_stage',
    'sex',
    'ethnicity',
    'disease',
    'tissue',
    'cell_type',
    'assay',
    #'suspension_type'
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
