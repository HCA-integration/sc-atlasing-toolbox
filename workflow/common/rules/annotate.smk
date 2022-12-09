rule add_obs:
    input:
        h5ad='test/data/Lee2020.h5ad',
        tsv='test/data/Lee2020.h5ad'
    output:
        h5ad='test/out/annotated.h5ad'
    params:
        h5ad_column='barcode',
        tsv_column='barcode_DCP',
    script:
        '../scripts/add_obs.py'
