rule marker_genes:
    input:
        zarr=rules.preprocessing_per_lineage_normalize.output.zarr,
        clusters=rules.integration_per_lineage_cluster.output.tsv
    output:
        tsv=out_dir / 'marker_genes' / paramspace.wildcard_pattern / 'lineage~{lineage}' / '{output_type}--{resolution}' / 'marker_genes.tsv',
        uns=out_dir / 'marker_genes' / paramspace.wildcard_pattern / 'lineage~{lineage}' / '{output_type}--{resolution}' / 'marker_genes.pkl',
        rankplot=image_dir / 'marker_genes' / paramspace.wildcard_pattern / 'lineage~{lineage}' / '{output_type}--{resolution}--rankplot.png',
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=get_resource(config,profile='cpu',resource_key='partition'),
        qos=get_resource(config,profile='cpu',resource_key='qos'),
        mem_mb=lambda w, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
        gpu=get_resource(config,profile='cpu',resource_key='gpu'),
    script:
        '../scripts/marker_genes.py'


rule marker_genes_collect:
    input: 
        lambda wildcards: expand(
            rules.marker_genes.output.tsv,
            resolution=get_for_dataset(
                config,
                wildcards.dataset,
                query=['clustering', 'resolutions'],
                default=[0.4, 0.8, 1.0]
            ),
            output_type=get_params(wildcards, parameters, 'output_type'),
            allow_missing=True
        )
    output:
        touch(out_dir / 'marker_genes' / paramspace.wildcard_pattern / 'lineage~{lineage}' / 'marker_genes.done')


rule marker_genes_per_lineage:
    input:
        unpack(lambda w: collect_lineages(w, rules.marker_genes_collect.output))
    output:
        touch(out_dir / 'marker_genes' / paramspace.wildcard_pattern / 'marker_genes.done')


rule marker_genes_all:
    input:
        expand(rules.marker_genes_per_lineage.output, zip, **parameters[wildcard_names].to_dict('list'))