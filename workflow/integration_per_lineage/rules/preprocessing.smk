pp_per_lineage_dir = out_dir / 'preprocessing' / 'dataset~{dataset}' / 'file_id~{file_id}' / 'lineage~{lineage}'


use rule normalize from preprocessing as preprocessing_per_lineage_normalize with:
    input:
        zarr=lambda w: get_checkpoint_output(checkpoints.split_lineage, **w) / f'lineage~{w.lineage}.zarr',
    output:
        zarr=directory(pp_per_lineage_dir / 'normalize.zarr'),
    params:
        raw_counts=lambda w: get_for_dataset(config, w.dataset, [module_name, 'raw_counts']),
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


use rule highly_variable_genes from preprocessing as preprocessing_per_lineage_highly_variable_genes with:
    input:
        zarr=rules.preprocessing_per_lineage_normalize.output.zarr,
    output:
        zarr=directory(pp_per_lineage_dir / 'highly_variable_genes.zarr'),
    params:
        args=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'highly_variable_genes']),
        batch=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'batch']),
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


use rule pca from preprocessing as preprocessing_per_lineage_pca with:
    input:
        zarr=rules.preprocessing_per_lineage_highly_variable_genes.output.zarr,
        counts=rules.preprocessing_per_lineage_highly_variable_genes.input.zarr,
    output:
        zarr=directory(pp_per_lineage_dir / 'pca.zarr'),
    params:
        scale=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'scale'])
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


use rule neighbors from preprocessing as preprocessing_per_lineage_neighbors with:
    input:
        zarr=rules.preprocessing_per_lineage_pca.output.zarr
    output:
        zarr=directory(pp_per_lineage_dir / 'neighbors.zarr'),
    params:
        args=lambda w: get_for_dataset(config, w.dataset, ['preprocessing', 'neighbors'], default={}),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_resource(config,profile='gpu',resource_key='mem_mb', attempt=attempt),


use rule assemble from preprocessing as preprocessing_per_lineage_assemble with:
    input:
        unpack(
            lambda wildcards: dict(
                counts=get_checkpoint_output(checkpoints.split_lineage, **wildcards) / f'lineage~{wildcards.lineage}.zarr',
                normalize=rules.preprocessing_per_lineage_normalize.output.zarr,
                highly_variable_genes=rules.preprocessing_per_lineage_highly_variable_genes.output.zarr,
                pca=rules.preprocessing_per_lineage_pca.output.zarr,
                neighbors=rules.preprocessing_per_lineage_neighbors.output.zarr,
            )
        )
    output:
        zarr=directory(pp_per_lineage_dir / 'assembled.zarr'),
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),

