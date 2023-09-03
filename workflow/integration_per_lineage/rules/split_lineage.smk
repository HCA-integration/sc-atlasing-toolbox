checkpoint split_lineage:
    message:
        """
        Split lineages: Split dataset "{wildcards.dataset}" by lineage key "{params.lineage_key}"
        input: {input}
        output: {output}
        wildcards: {wildcards}
        resources: {resources.mem_mb}MB
        """
    input:
        lambda wildcards: get_for_dataset(config, wildcards.dataset, query=['input', module_name])
    output:
        directory(out_dir / 'per_lineage' / 'split_lineage' / 'dataset~{dataset}')
    params:
        batch=lambda wildcards: get_for_dataset(config, wildcards.dataset, [module_name, 'batch']),
        label=lambda wildcards: get_for_dataset(config, wildcards.dataset, [module_name, 'label']),
        lineage_key=lambda wildcards: get_for_dataset(config, wildcards.dataset, [module_name, 'lineage']),
        # norm_counts=lambda wildcards: get_params(wildcards,parameters,'norm_counts'),
        hvg_args=lambda wildcards: get_for_dataset(config, wildcards.dataset, ['preprocessing', 'highly_variable_genes']),
    conda:
        get_env(config, 'scanpy', gpu_env='scanpy_rapids')
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
    script:
        '../scripts/split_anndata.py'
