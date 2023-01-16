rule zarr_to_h5ad:
    input:
        zarr='{file}.zarr'
    output:
        h5ad='{file}.h5ad'
    conda:
        '../envs/scanpy.yaml'
    script:
        '../scripts/convert_zarr_h5ad.py'
