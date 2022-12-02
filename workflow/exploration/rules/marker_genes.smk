rule marker_genes:
    input:
        h5ad=rules.load_data_filter.output.h5ad
    output:
        png=out_dir / 'marker_genes' / '{study}.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        markers=lambda wildcards: config['ORGANS']['blood']['marker_genes']  # TODO: organ per dataset
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/marker_genes.py'


rule marker_genes_all:
    input:
        expand(rules.marker_genes.output,study=dataset_df['study'].unique())
