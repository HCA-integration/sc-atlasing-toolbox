def get_neighbors_file(wildcards):
    if mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.label_harmonization_neighbors.output.zarr
    return mcfg.get_input_file(**wildcards)


def get_umap_file(wildcards):
    wildcards = {k: wildcards[k] for k in ('dataset', 'file_id')}
    if mcfg.get_from_parameters(wildcards, 'recompute_umap', default=False) \
    or mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False):
        return rules.label_harmonization_umap.output.zarr
    return get_neighbors_file(wildcards)


use rule neighbors from preprocessing as label_harmonization_neighbors with:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'neighbors.zarr'),
    params:
        args=lambda wildcards: {
            'use_rep': mcfg.get_from_parameters(wildcards, 'cellhint').get('use_rep', 'X_pca'),
        },
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),


use rule umap from preprocessing as label_harmonization_umap with:
    input:
        anndata=get_neighbors_file,
        rep=get_neighbors_file,
    output:
        zarr=directory(out_dir / paramspace.wildcard_pattern / 'umap.zarr'),
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='gpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
