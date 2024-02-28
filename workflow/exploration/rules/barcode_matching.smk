from utils.environments import get_env

# rule get_barcodes:
#     input:
#         zarr=rules.metadata.output.zarr
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
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        png=mcfg.image_dir / 'barcode_matching' / f'{params.wildcard_pattern}.png',
    conda:
        get_env(config, 'plots')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/barcode_matching.py'


rule barcode_matching_all:
    input:
        mcfg.get_output_files(rules.barcode_matching.output)