import scanpy as sc

organ = snakemake.wildcards.organ
files = snakemake.input
out_file = snakemake.output.h5ad


def read_adata(file):
    ad = sc.read(file)
    del ad.uns
    del ad.layers
    del ad.raw
    return ad


adatas = [read_adata(file) for file in files]
print(adatas)

print('Concatenate...')
adata = sc.concat(adatas)

print(f'Add metadata')
adata.uns['dataset'] = organ
adata.uns['organ'] = organ
adata.obs['organ'] = organ

print('Write...')
adata.write(out_file, compression='gzip')