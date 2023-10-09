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
        module_name: str,
        config: dict,
        parameters: pd.DataFrame = None,
        default_output: [str, Rule] = None,
        wildcard_names: list = None,
    ):
        super().__init__(
            module_name=module_name,
            config=config,
            default_output=default_output,
        )
        
        self.parameters = IntegrationWildcardParameters(
            module_name=module_name,
            parameters=parameters,
            input_file_wildcards=self.input_files.get_wildcards(),
            dataset_config=self.datasets,
            default_config=self.config.get('defaults'),
            output_directory=self.out_dir,
            wildcard_names=wildcard_names,
        )
        
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
            filename_params=['method', 'hyperparams'],
            filename_sep='--',
        )


    def get_hyperparams(self):
        return self.paramters.get_hyperparams()



class IntegrationWildcardParameters(WildcardParameters):
    
    def __init__(
        self,
        module_name: str,
        parameters: pd.DataFrame,
        input_file_wildcards: pd.DataFrame,
        dataset_config: dict,
        default_config: dict,
        output_directory: [str, Path],
        wildcard_names: list = None,
    ):
        self.out_dir = output_directory
        if wildcard_names is None:
            wildcard_names = ['dataset', 'file_id', 'batch', 'label', 'method', 'hyperparams']
        super().__init__(
            module_name=module_name,
            parameters=parameters,
            input_file_wildcards=input_file_wildcards,
            dataset_config=dataset_config,
            default_config=default_config,
            wildcard_names=wildcard_names,
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
        Set hyperparameters specific to each method of a module for all datasets

        :param config: config containing dataset specific information
        :param module: name of module, key must be present for each dataset entry
        """
        defaults = self.default_config.get(self.module_name)
        records = []
        for dataset, dataset_dict in self.dataset_config.items():
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


    def get_hyperparams(self):
        return self.hyperparams_df
