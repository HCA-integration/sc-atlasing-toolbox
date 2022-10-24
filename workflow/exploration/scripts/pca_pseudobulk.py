from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc
import sc_toolbox as sct


sc.set_figure_params(dpi=100, dpi_save=150, frameon=False)
plt.rcParams['figure.figsize'] = 12, 8
input_h5ad = snakemake.input.h5ad
output_png = snakemake.output.png
bulk_by = snakemake.params['bulk_by']
color = snakemake.params['color']

adata = sc.read(input_h5ad)

adata.X = adata.X.todense()
pbulks_df = sct.tools.generate_pseudobulk(adata, group_key=bulk_by)
pbulks_df = pbulks_df.set_index('Genes')

obs = pd.DataFrame(pbulks_df.columns, columns=[bulk_by])
obs = obs.merge(adata.obs[[bulk_by, color]].drop_duplicates(), on='sample')

adata_bulk = anndata.AnnData(
    pbulks_df.transpose().reset_index(drop=True),
    obs=obs,
    dtype='float32'
)
sc.pp.log1p(adata_bulk)
sc.pp.highly_variable_genes(adata_bulk, n_top_genes=200)
sc.pp.pca(adata_bulk)

print(output_png)
sc.pl.pca(
    adata_bulk,
    color=color,
    title='Pseudobulk per sample, colored by donor',
    show=False,
)
plt.tight_layout()
plt.savefig(output_png)