"""
Assemble h5ad with different Preprocessing outputs
"""

def collect_files(wildcards):
    file_dict = {
        'counts': get_for_dataset(config, wildcards.dataset, ['input', module_name]),
        'normalize': rules.normalize.output.h5ad,
        'highly_variable_genes': rules.highly_variable_genes.output.h5ad,
        'pca': rules.pca.output.h5ad,
        'neighbors': rules.neighbors.output.h5ad,
        # 'umap': rules.umap.output.h5ad,
    }
    assembly_config = get_for_dataset(config, wildcards.dataset, [module_name, 'assemble'])
    if assembly_config is None:
        return file_dict
    return {k: v for k, v in file_dict.items() if k in assembly_config}


rule assemble:
    input:
        unpack(collect_files)
    output:
        h5ad=out_dir / '{dataset}' / 'preprocessed.h5ad'
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/assemble.py'