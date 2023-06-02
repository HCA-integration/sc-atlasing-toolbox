import subprocess
import scanpy as sc

from utils import get_url

output_file = snakemake.output.h5ad
dataset_df = snakemake.params.dataset_df
tmpdir = snakemake.resources.tmpdir
wildcards = snakemake.wildcards

file_type = 'h5ad'
url = get_url(dataset_df, wildcards)
if isinstance(url, tuple):
    url, file_type = url

adata_file = f'{tmpdir}/{wildcards.dataset}.{file_type}'
print(adata_file)

subprocess.run(["wget", "-O", adata_file, url], check=True)

# read file and write to h5ad
if file_type == 'h5ad':
    subprocess.run(["cp", adata_file, output_file])
    exit(0)
elif file_type == 'loom':
    print('read loom...')
    adata = sc.read_loom(adata_file, sparse=True)
else:
    raise TypeError(f'Unknown file type {file_type}')

print('write to h5ad...')
adata.write(output_file, compression='lzf')
