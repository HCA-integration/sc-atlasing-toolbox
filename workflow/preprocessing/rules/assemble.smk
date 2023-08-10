"""
Assemble anndata with different Preprocessing outputs
"""

def collect_files(wildcards):
    file_dict = {
        'counts': get_for_dataset(config, wildcards.dataset, ['input', module_name]),
        'normalize': rules.normalize.output.zarr,
        'highly_variable_genes': rules.highly_variable_genes.output.zarr,
        'pca': rules.pca.output.zarr,
        'neighbors': rules.neighbors.output.zarr,
        # 'umap': rules.umap.output.zarr,
    }
    assembly_config = get_for_dataset(config, wildcards.dataset, [module_name, 'assemble'])
    if assembly_config is None:
        return file_dict
    return {k: v for k, v in file_dict.items() if k in assembly_config}


rule assemble:
    input:
        unpack(collect_files)
    output:
        zarr=directory(out_dir / '{dataset}' / 'preprocessed.zarr')
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/assemble.py'