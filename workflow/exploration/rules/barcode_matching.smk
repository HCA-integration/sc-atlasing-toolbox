# rule get_barcodes:
#     input:
#         zarr=rules.load_data_metadata.output.zarr
#     output:
#         tsv=out_dir / 'barcode_matching' / 'barcodes' / '{dataset}.tsv',
#     conda:
#         '../envs/plots.yaml'
#     script:
#         '../scripts/get_barcodes.py'


# rule barcode_matching:
#     input:
#         tsv=rules.get_barcodes.output.tsv
#     output:
#         png=images_dir / 'barcode_matching' / '{dataset}.png',
#     params:
#         min_size=5
#     conda:
#         '../envs/plots.yaml'
#     resources:
#         mem_mb=get_resource(config,profile='cpu',resource_key='mem_mb')
#     script:
#         '../scripts/barcode_matching.R'


rule barcode_matching:
    input:
        zarr=rules.load_data_merge_study.output.zarr
    output:
        png=images_dir / 'barcode_matching' / '{study}.png',
    conda:
        '../envs/plots.yaml'
    script:
        '../scripts/barcode_matching.py'


rule barcode_matching_all:
    input:
        expand(rules.barcode_matching.output,study=dataset_df['study'])