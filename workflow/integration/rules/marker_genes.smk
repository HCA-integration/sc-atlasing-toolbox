rule marker_genes:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id),
        clusters=rules.integration_cluster.output.tsv
    output:
        tsv=out_dir / 'marker_genes' / paramspace.wildcard_pattern / '{resolution}' / 'marker_genes.tsv',
        uns=out_dir / 'marker_genes' / paramspace.wildcard_pattern / '{resolution}' / 'marker_genes.pkl',
        rankplot=image_dir / 'marker_genes' / paramspace.wildcard_pattern / '{resolution}--rankplot.png',
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell')
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb', attempt=attempt),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
    script:
        '../scripts/marker_genes.py'


rule marker_genes_collect:
    input: 
        lambda wildcards: expand(
            rules.marker_genes.output.tsv,
            resolution=mcfg.get_for_dataset(
                wildcards.dataset,
                query=['clustering', 'resolutions'],
                default=[0.4, 0.8, 1.0]
            ),
            allow_missing=True
        )
    output:
        touch(out_dir / 'marker_genes' / paramspace.wildcard_pattern / 'marker_genes.done')


rule marker_genes_all:
    input:
        mcfg.get_output_files(rules.marker_genes_collect.output)