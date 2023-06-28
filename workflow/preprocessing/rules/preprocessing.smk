"""
Preprocessing steps
Only the unique outputs per step are saved for storage efficiency. For assembled h5ad files, see `assemble.smk`.
"""

rule normalize:
    input:
        lambda w: get_for_dataset(config, w.dataset, ['input', module_name])
    output:
        h5ad=out_dir / '{dataset}' / 'normalized.h5ad'
    params:
        raw_counts=lambda w: get_for_dataset(config, w.dataset, [module_name, 'raw_counts']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/normalize.py'


rule highly_variable_genes:
    input:
        h5ad=rules.normalize.output.h5ad
    output:
        h5ad=out_dir / '{dataset}' / 'highly_variable_genes.h5ad'
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'highly_variable_genes']),
        batch=lambda w: get_for_dataset(config, w.dataset, [module_name, 'batch']),
        lineage=lambda w: get_for_dataset(config, w.dataset, [module_name, 'lineage']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/highly_variable_genes.py'


rule pca:
    input:
        h5ad=rules.highly_variable_genes.output.h5ad,
        counts=rules.normalize.output.h5ad,
    output:
        h5ad=out_dir / '{dataset}' / 'pca.h5ad'
    params:
        scale=lambda w: get_for_dataset(config, w.dataset, [module_name, 'scale'])
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda: '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/pca.py'


rule neighbors:
    input:
        h5ad=rules.pca.output.h5ad
    output:
        h5ad=out_dir / '{dataset}' / 'neighbors.h5ad'
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'neighbors']),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
    conda:
        ifelse(
            'os' not in config.keys() or config['os'] == 'm1',
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/neighbors.py'


rule umap:
    input:
        h5ad=rules.neighbors.output.h5ad
    output:
        h5ad=out_dir / '{dataset}' / 'umap.h5ad'
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
    conda:
        ifelse(
            'os' not in config.keys() or config['os'] == 'm1',
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/umap.py'
