import scanpy as sc
import pandas as pd


def write_metrics(filename, scores, output_types, **kwargs):
    """
    Write metrics for output type specific scores
    :param filename: file to write to
    :param scores: list of scores, order must match that of output_types
    :param output_types: list of output types, order must match that of scores
    :param kwargs: additional information to add to output
    """
    meta_names = list(kwargs.keys())
    meta_values = [kwargs[col] for col in meta_names]

    records = [
        (*meta_values, output_type, score) 
        for score, output_type in zip(scores, output_types)
    ]
    df = pd.DataFrame.from_records(records, columns=meta_names + ['output_type', 'score'])
    df.to_csv(filename, sep='\t', index=False)


def rename_categories(adata, obs_col):
    s = adata.obs[obs_col]
    s = s.cat.rename_categories({i for i, _ in enumerate(s.cat.categories)})
    return s.to_numpy()


def select_neighbors(adata, output_type):
    # neighbors_key = f'neighbors_{output_type}'
    neighbors_key = 'neighbors'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    
    connectivities_key = adata.uns[neighbors_key]['connectivities_key']
    assert connectivities_key in adata.obsp, f'Connectivities key "{connectivities_key}" missing from adata.obsp {adata}'
    adata.obsp['connectivities'] = adata.obsp[connectivities_key]
    
    distances_key = adata.uns[neighbors_key]['distances_key']
    assert distances_key in adata.obsp, f'Distances key "{distances_key}" missing from adata.obsp {adata}'
    adata.obsp['distances'] = adata.obsp[distances_key]
    
    return adata
