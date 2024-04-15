from collections.abc import Mapping

from .ModuleConfig import ModuleConfig


def get_marker_gene_set(mcfg: ModuleConfig, wildcards: Mapping):
    """
    Assumes that the marker genes are stored in the config file under the MARKER_GENES key and that the key is configured as a parameter
    """
    gene_set_key = mcfg.get_from_parameters(wildcards, 'marker_genes', default='')
    assert isinstance(gene_set_key, str), f'Expected string, got {gene_set_key}'
    return mcfg.config.get('MARKER_GENES', {}).get(gene_set_key, {})