rule plot:
    input:
        zarr=rules.run_method.output.zarr
    output:
        histplot_path=mcfg.image_dir / paramspace.wildcard_pattern / 'distances_histplot.png'
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu', resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/distances_plot.py'