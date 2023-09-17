import numpy as np
import pandas as pd
from snakemake.utils import Paramspace
from snakemake.rules import Rule

from utils.ModuleConfig import ModuleConfig
from utils.config import _get_or_default_from_config
from utils.misc import expand_dict_and_serialize, unique_dataframe


class IntegrationConfig(ModuleConfig):
    
    def __init__(self, module_name: str, config: dict, parameters: pd.DataFrame = None, default_output: [str, Rule] = None):
        super().__init__(module_name, config, parameters, default_output)
        self.wildcard_names = ['dataset', 'file_id', 'batch', 'label', 'method', 'hyperparams']
        
        # remove redundant label wildcards
        if 'use_cell_type' in self.wildcards_df.columns:
            self.wildcards_df['label'] = np.where(
                self.wildcards_df['use_cell_type'],
                self.wildcards_df['label'],
                'None'
            )
            self.wildcards_df = unique_dataframe(self.wildcards_df)

        self.paramspace = Paramspace(
            self.wildcards_df[self.wildcard_names],
            filename_params=['method', 'hyperparams'],
            filename_sep='--',
        )
    
    
    def set_wildcards(
        self,
        config_params: list = None,
        wildcard_names: list = None,
        explode_by: list = None,
        config_keys: list = None,
        warn: bool = False,
    ):
        # parameters = pd.read_table(workflow.source_path('params.tsv'))
        # parameters['output_type'] = parameters['output_type'].str.split(',')
        self.set_hyperparams()
        super().set_wildcards(
            config_params=['methods', 'label', 'batch', 'norm_counts', 'raw_counts'],
            wildcard_names=['method', 'label', 'batch', 'norm_counts', 'raw_counts'],
            explode_by=['method', 'batch'],
        )#.merge(parameters,on='method')
        self.wildcards_df = self.wildcards_df.merge(
            self.hyperparams_df,
            on=['dataset', 'method'],
            how='left'
        )


    def set_hyperparams(self):
        """
        Get hyperparameters specific to each method of a module for all datasets

        :param config: config containing dataset specific information
        :param module: name of module, key must be present for each dataset entry
        """
        records = []
        for dataset, dataset_dict in self.datasets.items():
            methods_config = _get_or_default_from_config(
                config=dataset_dict,
                defaults=self.config['defaults'][self.module_name],
                key=self.module_name,
                value='methods',
                update=False,
            )
            if methods_config is None:
                continue
            for method, hyperparams_dict in methods_config.items():
                if isinstance(hyperparams_dict, dict):
                    records.extend(
                        (dataset, method, *rec)
                        for rec in expand_dict_and_serialize(hyperparams_dict)
                    )
                else:
                    records.append((dataset, method, str(hyperparams_dict), hyperparams_dict))
        self.hyperparams_df = pd.DataFrame(
            records,
            columns=['dataset', 'method', 'hyperparams', 'hyperparams_dict']
        )
        
        # write to file 
        unique_dataframe(
            self.hyperparams_df[['method', 'hyperparams', 'hyperparams_dict']]
        ).to_csv(self.out_dir / 'hyperparams.tsv', sep='\t', index=False)



    def get_hyperparams(self):
        return self.hyperparams_df
