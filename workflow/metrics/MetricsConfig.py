import pandas as pd
from pprint import pformat

from utils.ModuleConfig import ModuleConfig
from utils.marker_genes import get_marker_gene_set
from utils.misc import unique_dataframe


class MetricsConfig(ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        super().__init__(**kwargs, write_output_files=False)
        
        self.set_gene_sets()
        self.set_covariates()
        
        self.update_parameters(
            wildcards_df=unique_dataframe(self.parameters.wildcards_df),
            wildcard_names=self.parameters.wildcard_names,
            **self.parameters.paramspace_kwargs
        )


    def set_gene_sets(self, config_key='gene_sets'):
        gene_set_map = dict()
        wildcards_df = self.parameters.wildcards_df
        
        for dataset in self.datasets:
            gene_sets = get_marker_gene_set(
                self,
                wildcards={'dataset': dataset},
                config_key=config_key,
                flatten=True
            )
            gene_set_map[dataset] = gene_sets
            
            # convert into format for exploding
            # gene_sets = [{key: value} for key, value in gene_sets.items()]
            gene_sets = list(gene_sets.items())
            
            # set gene sets only for metrics that use it
            dataset_mask = wildcards_df['dataset'] == dataset
            use_config_mask = wildcards_df[f'use_{config_key}'] & dataset_mask

            wildcards_df.loc[dataset_mask, config_key] = None  # Default to None
            wildcards_df.loc[use_config_mask, config_key] = pd.Series(
                [gene_sets] * use_config_mask.sum(),
                dtype=object
            )

        # explode gene sets so that each metric setting has a single gene set
        self.wildcards_df = wildcards_df.explode(config_key)
        self.gene_set_map = gene_set_map
    
    
    def get_gene_sets(self, dataset, **kwargs):
        return self.gene_set_map[dataset]


    def set_covariates(self, config_key='covariates'):
        wildcards_df = self.wildcards_df
        wildcards_df.loc[~wildcards_df[f'use_{config_key}'], config_key] = None
        self.wildcards_df = wildcards_df.explode(config_key)