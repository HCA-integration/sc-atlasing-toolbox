from pathlib import Path
from pandas.api.types import is_numeric_dtype
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
# from matplotlib.axes import Axes
import seaborn as sns
import scanpy as sc
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
sns.set_theme(style='white')
sc.set_figure_params(frameon=False, fontsize=10, dpi_save=300, vector_friendly=True)

from utils.io import read_anndata
from qc_utils import parse_parameters, get_thresholds, plot_qc_joint


input_zarr = snakemake.input.zarr
output_joint = snakemake.output.joint

output_joint = Path(output_joint)
output_joint.mkdir(parents=True, exist_ok=True)

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')

# get parameters
file_id = snakemake.wildcards.file_id
sample, dataset, hues = parse_parameters(adata, snakemake.params)
hues = hues+['percent_mito']

# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('No data, skip plotting...')
    # plt.savefig(output_density)
    # plt.savefig(output_density_log)
    exit()

thresholds = get_thresholds(
    threshold_keys=['n_counts', 'n_genes', 'percent_mito'],
    autoqc_thresholds=adata.uns['scautoqc_ranges'],
    user_thresholds=snakemake.params.get('thresholds'),
    user_key=file_id,
)
logging.info(f'\n{pformat(thresholds)}')

# subset to max of 100k cells due to high computational cost
density_data = adata.obs.sample(n=int(min(1e5, adata.n_obs)), random_state=42)

coordinates = [
    ('n_counts', 'n_genes', 10, 10),
    ('n_genes', 'percent_mito', 2, 1),
]

scatter_plot_kwargs = dict(
    s=4,
    alpha=.5,
    linewidth=0,
)

kde_plot_kwargs = dict(
    fill=True,
    cmap='plasma',
    alpha=.8,
)

for x, y, log_x, log_y in coordinates:
    logging.info(f'Joint QC plots per {x} vs {y}...')
    
    density_png = output_joint / f'{x}_vs_{y}_density.png'
    density_log_png = output_joint / f'log_{x}_vs_log_{y}_density.png'
    
    plot_qc_joint(
        density_data,
        x=x,
        y=y,
        # main_plot_function=Axes.hexbin,
        main_plot_function=sns.kdeplot,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title='',
        **kde_plot_kwargs,
    )
    plt.tight_layout()
    plt.savefig(density_png)
    
    plot_qc_joint(
        density_data,
        x=x,
        y=y,
        main_plot_function=sns.kdeplot,
        log_x=log_x,
        log_y=log_y,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title='',
        **kde_plot_kwargs,
    )
    plt.tight_layout()
    plt.savefig(density_log_png)

    for hue in hues:
        logging.info(f'Joint QC plots for hue={hue}...')
        joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
        
        # determine plotting parameters
        if is_numeric_dtype(adata.obs[hue]):
            palette = 'plasma'
            legend = 'brief'
        else:
            palette = None # if adata.obs[hue].nunique() > 100 else 'plasma'
            legend = adata.obs[hue].nunique() <= 20
        scatter_plot_kwargs |= dict(palette=palette, legend=legend)
        
        # plot joint QC on regular scale
        png_file = output_joint / f'{x}_vs_{y}_hue={hue}.png'
        plot_qc_joint(
            adata.obs,
            x=x,
            y=y,
            hue=hue,
            marginal_hue=hue,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **scatter_plot_kwargs,
        )
        plt.tight_layout()
        plt.savefig(png_file)
        
        # assemble figure
        f, axes = plt.subplots(1, 2, figsize=(20, 10))
        axes[0].imshow(mpimg.imread(png_file))
        axes[1].imshow(mpimg.imread(density_png))
        for ax in axes.ravel():
            ax.set_axis_off()
        plt.suptitle(joint_title, fontsize=16)
        plt.tight_layout()
        plt.savefig(png_file)

        # plot in log scale 
        png_file = output_joint / f'log_{x}_vs_log_{y}_hue={hue}.png'
        g_scatter = plot_qc_joint(
            adata.obs,
            x=x,
            y=y,
            log_x=log_x,
            log_y=log_y,
            hue=hue,
            marginal_hue=hue,
            x_threshold=thresholds[x],
            y_threshold=thresholds[y],
            title='',
            **scatter_plot_kwargs,
        )
        plt.tight_layout()
        plt.savefig(png_file)
        
        # assemble figure
        f, axes = plt.subplots(1, 2, figsize=(20, 10))
        axes[0].imshow(mpimg.imread(png_file))
        axes[1].imshow(mpimg.imread(density_log_png))
        for ax in axes.ravel():
            ax.set_axis_off()
        plt.suptitle(joint_title, fontsize=16)
        plt.tight_layout()
        plt.savefig(png_file)
    
    # remove redundant plots
    density_png.unlink()
    density_log_png.unlink()
    plt.close('all')