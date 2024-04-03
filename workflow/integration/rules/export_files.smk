rule export_files:
    input:
        zarr=rules.integration_compute_umap.output.zarr,
        clusters=rules.integration_cluster_merge.output.tsv,
    output:
        zarr=directory(out_dir / 'export' / paramspace.wildcard_pattern / 'adata.zarr'),
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    script:
        '../scripts/export_files.py'


rule export_files_all:
    input:
        mcfg.get_output_files(rules.export_files.output)
    localrule: True
