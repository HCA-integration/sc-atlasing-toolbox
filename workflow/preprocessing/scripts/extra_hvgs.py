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
import anndata as ad

from utils.io import read_anndata, write_zarr_linked, csr_matrix_int64_indptr
from utils.misc import dask_compute
from utils.processing import filter_genes, sc, USE_GPU


def match_genes(var_df, gene_list, column=None):
    import urllib

    def is_url(url):
        try:
            result = urllib.parse.urlparse(url)
            return all([result.scheme, result.netloc])
        except ValueError:
            return False

    genes_from_path = dict()
    for gene in gene_list:
        if Path(gene).exists():
            with open(gene, 'r') as f:
                genes_from_path[gene] = f.read().splitlines()
        elif is_url(gene):
            try:
                with urllib.request.urlopen(gene) as f:
                    genes_from_path[gene] = f.read().decode('utf-8').splitlines()
            except Exception as e:
                logging.error(f'Error reading gene list from URL {gene}...')
                raise e

    for path, genes in genes_from_path.items():
        logging.info(f'Gene list from {path}: {len(genes)} genes')
        gene_list.extend(genes)
        gene_list.remove(path)

    try:
        genes = var_df.index if column is None else var_df[column]
        pattern = '|'.join(gene_list)
        return genes[genes.astype(str).str.contains(pattern, regex=True)].index
    except Exception as e:
        logging.error(f'Error: {e}')
        logging.error(f'Gene list: {gene_list}')
        logging.error(f'Pattern: {pattern}')
        logging.error(f'Gene names: {var_df.index}')
        raise e


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

logging.info(f'args: {args}')

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
    var = adata.var.copy()
    
    # workaround for CxG datasets
    feature_col = 'feature_name' if 'feature_name' in var.columns else None

    # remove user-specified genes
    if remove_genes:
        remove_genes = match_genes(var, remove_genes, column=feature_col)
        logging.info(f'Remove {len(remove_genes)} genes (subset data)...')
        adata = adata[:, ~adata.var_names.isin(remove_genes)].copy()

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
            
            min_cells = 10
            if _ad.n_obs < min_cells:
                logging.info(f'Group={group} has fewer than {min_cells} cells, skipping...')
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

    var['extra_hvgs'] = False
    var.loc[adata.var_names, 'extra_hvgs'] = adata.var['extra_hvgs']
    adata = ad.AnnData(var=var, uns=adata.uns)

    # add user-provided genes
    if extra_genes:
        n_genes = len(extra_genes)
        extra_genes = match_genes(var, extra_genes, column=feature_col)
        
        if len(extra_genes) < n_genes:
            logging.warning(f'Only {len(extra_genes)} of {n_genes} user-provided genes found in data...')
        if len(extra_genes) == 0:
            logging.info('No extra user genes added...')
        else:
            logging.info(f'Add {len(extra_genes)} user-provided genes...')
            var.loc[extra_genes, 'extra_hvgs'] = True

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns', 'var']
)