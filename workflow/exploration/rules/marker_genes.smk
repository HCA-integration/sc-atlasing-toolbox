rule marker_genes:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        png=mcfg.image_dir / 'marker_genes' / f'{params.wildcard_pattern}.png',
    params:
        dataset=lambda wildcards: wildcards.file_id,
        markers=lambda wildcards: get_marker_gene_set(mcfg, wildcards),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/marker_genes.py'


rule marker_genes_all:
    input:
        mcfg.get_output_files(rules.marker_genes.output)
    localrule: True
