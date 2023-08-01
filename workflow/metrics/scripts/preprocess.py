import scanpy as sc
import muon as mu
import logging
logging.basicConfig(level=logging.INFO)

from metrics.utils import compute_neighbors, get_from_adata
from utils.io import read_anndata_or_mudata


input_adata = snakemake.input.h5ad
output_file = snakemake.output.h5mu
lineage_key = snakemake.wildcards.lineage_key

logging.info('Read file...')
adata = read_anndata_or_mudata(input_adata)
meta = get_from_adata(adata)

if isinstance(adata, mu.MuData):
    logging.info('Data is already a MuData object')
    mudata = adata
elif lineage_key not in adata.obs.columns:
    logging.info('Data is global AnnData object, use generic lineage name.')
    mudata = mu.MuData({lineage_key: adata})
else:
    logging.info('Data is AnnData object, split by lineage.')
    mudata = mu.MuData(
        {
            str(lineage): adata[adata.obs[lineage_key] == lineage]
             for lineage in adata.obs[lineage_key].unique()
        }
    )
mudata.uns = adata.uns

for lineage in mudata.mod:
    logging.info(f'Processing lineage {lineage}...')
    ad = mudata[lineage]
    for output_type in meta['output_types']:
        logging.info(f'Computing neighbors for output type {output_type}...')
        # compute neighbors by output type
        ad_knn = compute_neighbors(ad, output_type)
        ad.obsp[f'connectivities_{output_type}'] = ad_knn.obsp['connectivities']
        ad.obsp[f'distances_{output_type}'] = ad_knn.obsp['distances']

logging.info(mudata.__str__())
mudata.write(output_file)