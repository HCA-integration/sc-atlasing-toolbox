from collections.abc import Mapping

from .ModuleConfig import ModuleConfig


def get_marker_gene_set(mcfg: ModuleConfig, wildcards: Mapping):
    """
    Assumes that the marker genes are stored in the config file under the MARKER_GENES key and that the key is configured as a parameter
    """
    gene_set_keys =  mcfg.get_from_parameters(
        wildcards,
        'marker_genes',
        default=''
    ).split(',')  # split by comma for multiple gene sets
    
    # strip trailing whitespace
    gene_set_keys = [g.strip() for g in gene_set_keys]
    
    marker_genes = {}
    for key in gene_set_keys:
        assert isinstance(key, str), f'Expected string, got {key}'
        marker_genes[key] = mcfg.config.get('MARKER_GENES', {}).get(key, {})
    return marker_genes