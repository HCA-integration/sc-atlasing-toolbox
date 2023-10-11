from pathlib import Path
import numpy as np
import pandas as pd
from snakemake.rules import Rule

from utils.ModuleConfig import ModuleConfig
from utils.WildcardParameters import WildcardParameters
from utils.config import _get_or_default_from_config
from utils.misc import expand_dict_and_serialize, unique_dataframe


class IntegrationConfig(ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        # set default paramspace arguments
        if kwargs.get('paramspace_kwargs') is None:
            kwargs['paramspace_kwargs'] = dict(
                filename_params=['method', 'hyperparams'],
                filename_sep='--',
            )
        
        # parse output type info
        parameters = kwargs.get('parameters')
        if isinstance(parameters, str):
            parameters = pd.read_table(parameters)
        if isinstance(parameters, pd.DataFrame):
            parameters['output_type'] = parameters['output_type'].str.split(',')
            kwargs['parameters'] = parameters
        
        super().__init__(**kwargs)
        
        # set hyperparameters
        self.set_hyperparams()
        
        # remove redundant label wildcards
        wildcards_df = self.parameters.wildcards_df
        if 'use_cell_type' in wildcards_df.columns:
            wildcards_df['label'] = np.where(
                wildcards_df['use_cell_type'],
                wildcards_df['label'],
                'None'
            )
        self.parameters.update(
            wildcards_df=unique_dataframe(wildcards_df),
            wildcard_names=self.parameters.wildcard_names + ['hyperparams'],
        )


    def set_hyperparams(self):
        """
        Set hyperparameters specific to each method of a module for all datasets

        :param config: config containing dataset specific information
        :param module: name of module, key must be present for each dataset entry
        """
        defaults = self.parameters.default_config.get(self.module_name)
        records = []
        for dataset, dataset_dict in self.parameters.dataset_config.items():
            methods_config = _get_or_default_from_config(
                config=dataset_dict,
                defaults=defaults,
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
        
        # merge to wildcards
        self.parameters.wildcards_df = self.parameters.wildcards_df.merge(
            self.hyperparams_df,
            on=['dataset', 'method'],
            how='left'
        )


    def get_hyperparams(self):
        return self.hyperparams_df
