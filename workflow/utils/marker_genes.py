from collections.abc import Mapping

from .ModuleConfig import ModuleConfig


# TODO: rename to get_gene_sets
def get_marker_gene_set(
    mcfg: ModuleConfig,
    wildcards: Mapping,
    config_key: str = 'marker_genes',
    gene_map_key: str = 'MARKER_GENES',
    flatten: bool = False,
    **kwargs
) -> dict:
    """
    Assumes that the marker genes are stored in the config file under the MARKER_GENES key and that the key is configured as a parameter
    """
    gene_sets = mcfg.get_from_parameters(
        wildcards,
        config_key,
        default='',
        **kwargs
    )
    
    genes = dict()
    
    if isinstance(gene_sets, dict):
        genes = gene_sets
    elif isinstance(gene_sets, str):
        # split by comma for multiple gene sets
        gene_sets.split(',')
        
        # strip trailing whitespace
        gene_sets = [g.strip() for g in gene_sets]
        
        for key in gene_sets:
            assert isinstance(key, str), f'Expected string, got {key}'
            genes[key] = mcfg.config.get(gene_map_key, {}).get(key, {})
    
    if flatten:
        # ignore supersets of gene set maps to single flattened gene set map
        # Note: this will overwrite gene sets with the same key across different supersets
        for key, gene_set in genes.items():
            if isinstance(gene_set, str):
                gene_set = {key: [gene_set]}
            elif isinstance(gene_set, list):
                continue
            genes |= gene_set
    
    return genes