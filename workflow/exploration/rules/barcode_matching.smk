rule get_barcodes:
    input:
        h5ad=rules.load_data_metadata.output.h5ad
    output:
        tsv=out_dir / 'barcode_matching' / 'barcodes' / '{dataset}.tsv',
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/get_barcodes.py'


rule barcode_matching:
    input:
        tsv=rules.get_barcodes.output.tsv
    output:
        png=out_dir / 'barcode_matching' / '{dataset}.png',
    params:
        min_size=5
    conda:
        '../envs/plots.yaml'
    resources:
        mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/barcode_matching.R'


rule barcode_matching_all:
    input:
        expand(rules.barcode_matching.output,dataset=dataset_df['dataset'])