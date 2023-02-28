import scanpy as sc
import muon as mu
from metrics.utils import compute_neighbors, get_from_adata
from utils.io import read_anndata_or_mudata


input_adata = snakemake.input.h5ad
output_file = snakemake.output.h5mu
lineage_key = snakemake.wildcards.lineage_key

adata = read_anndata_or_mudata(input_adata)
meta = get_from_adata(adata)

if isinstance(adata, mu.MuData):
    mudata = adata
elif lineage_key not in adata.obs.columns:
    mudata = mu.MuData({lineage_key: adata})
else:
    mudata = mu.MuData(
        {
            str(lineage): adata[adata.obs[lineage_key] == lineage]
             for lineage in adata.obs[lineage_key].unique()
        }
    )
mudata.uns = adata.uns

for lineage in mudata.mod:
    ad = mudata[lineage]
    for output_type in meta['output_types']:
        # compute neighbors by output type
        ad_knn = compute_neighbors(ad, output_type)
        ad.obsp['connectivities_' + output_type] = ad_knn.obsp['connectivities']
        ad.obsp['distances_' + output_type] = ad_knn.obsp['distances']

print(mudata)
mudata.write(output_file)