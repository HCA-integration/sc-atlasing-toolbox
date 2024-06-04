"""
Highly variable gene selection
- HVG by group -> take union of HVGs from each group
- allow including user-specified genes
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, write_zarr_linked, csr_matrix_int64_indptr
from utils.misc import dask_compute
from utils.processing import filter_genes, sc, USE_GPU


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
extra_hvg_args = snakemake.params.get('extra_hvgs', {})
overwrite_args = extra_hvg_args.get('overwrite_args', {})
union_over = extra_hvg_args.get('union_over')
extra_genes = extra_hvg_args.get('extra_genes', [])
remove_genes = extra_hvg_args.get('remove_genes', [])

if args is None:
    args = {}
elif isinstance(args, dict):
    args.pop('subset', None) # don't support subsetting
if overwrite_args:
    args |= overwrite_args

logging.info(str(args))

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    backed=True,
    dask=True,
)
logging.info(adata.__str__())

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}
adata.uns['preprocessing']['extra_hvgs'] = args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.var['extra_hvgs'] = True
    adata.write_zarr(output_file)
    exit(0)

if args == False:
    logging.info('No highly variable gene parameters provided, including all genes...')
    adata.var['extra_hvgs'] = True
else:
    adata.var['extra_hvgs'] = False
    # union over groups
    if union_over is not None:
        if isinstance(union_over, str):
            union_over = [union_over]
        adata.obs['union_over'] = adata.obs[union_over].astype(str).apply(lambda x: '--'.join(x), axis=1)
        for group in adata.obs['union_over'].unique():
            logging.info(f'Subset to group={group}...')
            _ad = dask_compute(adata[adata.obs['union_over'] == group].copy())
            logging.info(f'{_ad.n_obs} cells, {_ad.n_vars} genes')
            
            logging.info('Filter genes...')
            _ad = filter_genes(
                _ad,
                min_cells=1,
                batch_key=args.get('batch_key'),
            )
            
            if _ad.n_obs < 2:
                logging.info(f'Group={group} has less than 2 cells, skipping...')
                continue
            
            if USE_GPU:
                sc.get.anndata_to_GPU(_ad)

            logging.info(f'Select features for group={group} with arguments: {args}...')
            sc.pp.highly_variable_genes(_ad, **args)
            
            # get union of gene sets
            adata.var['extra_hvgs'] = adata.var['extra_hvgs'] | _ad.var['highly_variable']
            del _ad
        del adata.obs['union_over']
    else:
        # default gene selection
        logging.info(f'Select features for all cells with arguments: {args}...')
        adata = dask_compute(adata)
        if USE_GPU:
            sc.get.anndata_to_GPU(adata)
        sc.pp.highly_variable_genes(adata, **args)
        adata.var['extra_hvgs'] = adata.var['highly_variable']

    # workaround for CxG datasets
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    # add user-provided genes
    if extra_genes is not None:
        n_genes = len(extra_genes)
        extra_genes = [gene for gene in extra_genes if gene in adata.var_names]
        
        if len(extra_genes) < n_genes:
            logging.warning(f'Only {len(extra_genes)} of {n_genes} user-provided genes found in data...')
        if extra_genes == 0:
            logging.info('No extra user genes added...')
        else:
            adata.var.loc[extra_genes, 'extra_hvgs'] = True

    # remove user-specified genes
    if remove_genes is not None:
        remove_genes = [gene for gene in remove_genes if gene in adata.var_names]
        logging.info(f'Remove {len(remove_genes)} genes...')
        adata.var.loc[remove_genes, 'extra_hvgs'] = False

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns', 'var']
)