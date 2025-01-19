rule prepare:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'prepare.zarr'),
    params:
        label_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label'),
        neighbor_args=lambda wildcards: mcfg.get_for_dataset(wildcards.dataset, ['preprocessing', 'neighbors'], default={}),
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        corrected_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'corrected', default='X'),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile='gpu'),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile='gpu'),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile='gpu', attempt=attempt),
        time="1-00:00:00",
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output)
    localrule: True


# Clustering for cluster-based metrics

use rule cluster from clustering as metrics_cluster with:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'cluster_resolutions' / '{algorithm}--{resolution}--{level}.zarr'),
    conda:
        get_env(config, 'scanpy') # force clustering on GPU
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),


use rule merge from clustering as metrics_cluster_collect with:
    input:
        zarr=rules.prepare.output.zarr,
        cluster_anno=expand(
            rules.metrics_cluster.output.zarr,
            resolution=[2 * (x + 1) / 10 for x in range(10)],
            algorithm='leiden',
            level=1,
            allow_missing=True,
        ),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'clustered.zarr'),


rule cluster_all:
    input:
        mcfg.get_output_files(rules.metrics_cluster_collect.output)
    localrule: True


# Gene scoring for gene-based metrics

rule score_genes:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'gene_scores' / '{gene_set}.zarr'),
    params:
        gene_set=lambda wildcards: mcfg.get_gene_sets(wildcards.dataset).get(wildcards.gene_set, []),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
    script:
        '../scripts/score_genes.py'


def get_score_genes(mcfg, datasets=None):
    if datasets is None:
        datasets = mcfg.get_datasets()
    targets = []
    for dataset in datasets:
        gene_set = mcfg.get_gene_sets(dataset)
        files = mcfg.get_output_files(
            rules.score_genes.output,
            extra_wildcards={'gene_set': gene_set}
        )
        targets.extend(files)
    return targets


rule score_genes_all:
    input: *get_score_genes(mcfg)
    localrule: True