def get_markers(wildcards):
    organ = dataset_df[dataset_df['study'] ==  wildcards.study]['organ'].tolist()[0]
    return config['ORGANS'][organ]['marker_genes']


rule marker_genes:
    input:
        zarr=rules.load_data_filter.output.zarr
    output:
        png=out_dir / 'marker_genes' / '{study}.png',
    params:
        dataset=lambda wildcards: wildcards.study,
        markers=get_markers
    conda:
        '../envs/scanpy.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/marker_genes.py'


rule marker_genes_all:
    input:
        expand(rules.marker_genes.output,study=dataset_df['study'].unique())
