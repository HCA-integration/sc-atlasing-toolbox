from typing import Union
from pathlib import Path
import pandas as pd
from snakemake.io import expand
from snakemake.rules import Rule

from .WildcardParameters import WildcardParameters
from .InputFiles import InputFiles
from .config import get_from_config, _get_or_default_from_config


class ModuleConfig:

    # module_name = ''
    # config = {}
    # datasets = {}
    # wildcard_names = []
    # parameters = WildcardParameters()
    # input_files = InputFiles()
    # default_output = None
    # out_dir = Path()
    # image_dir = Path()


    def __init__(
        self,
        module_name: str,
        config: dict,
        parameters: pd.DataFrame = None,
        default_output: [str, Rule] = None,
        wildcard_names: list = None,
    ):
        self.module_name = module_name
        self.default_output = default_output
        self.config = config
        self.set_defaults()
        
        self.out_dir = Path(self.config['output_dir']) / self.module_name
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.image_dir = Path(self.config['images']) / self.module_name

        # get subset of datasets for this module
        all_datasets = self.config['DATASETS']
        default_datasets = self.config['defaults'].get('datasets', all_datasets.keys())
        self.datasets = {
            dataset: entry for dataset, entry in all_datasets.items()
            if dataset in default_datasets
            and module_name in entry.get('input', {}).keys()
        }
        for dataset in self.datasets:
            self.set_defaults_per_dataset(dataset)

        self.input_files = InputFiles(
            module_name=module_name,
            dataset_config=self.datasets,
            output_directory=self.out_dir,
        )
        
        if wildcard_names is None:
            wildcard_names = ['dataset', 'file_id']
        self.parameters = WildcardParameters(
            module_name=module_name,
            parameters=parameters,
            input_file_wildcards=self.input_files.get_wildcards(),
            dataset_config=self.datasets,
            default_config=self.config.get('defaults'),
            wildcard_names=wildcard_names,
        )


    def set_defaults(self, warn: bool = False):
        self.config['DATASETS'] = self.config.get('DATASETS', {})
        self.config['defaults'] = self.config.get('defaults', {})
        self.config['defaults'][self.module_name] = self.config['defaults'].get(self.module_name, {})
        datasets = list(self.config['DATASETS'].keys())
        self.config['defaults']['datasets'] = self.config['defaults'].get('datasets', datasets)


    def set_defaults_per_dataset(self, dataset: str, warn: bool = False):
        # update entries for each dataset
        entry = _get_or_default_from_config(
            config=self.config['DATASETS'],
            defaults=self.config['defaults'],
            key=dataset,
            value=self.module_name,
            return_missing={},
            warn=warn,
            update=True,
        )
        
        # for TSV input make sure integration methods have the proper types TODO: deprecate
        if self.module_name == 'integration' and isinstance(entry, list):
            entry = {k: self.config['defaults'][self.module_name][k] for k in entry}
        
        # set entry in config
        self.config['DATASETS'][dataset][self.module_name] = entry
        self.datasets[dataset][self.module_name] = entry


    def get_for_dataset(
        self,
        dataset: str,
        query: list,
        default: Union[str,bool,float,int,dict,list, None] = None,
        warn: bool = False,
    ) -> Union[str,bool,float,int,dict,list, None]:
        """Get any key from the config via query

        Args:
            dataset (str): dataset key in config['DATASETS']
            query (list): list of keys to walk down the config
            default (Union[str,bool,float,int,dict,list, None], optional): default value if key not found. Defaults to None.

        Returns:
            Union[str,bool,float,int,dict,list, None]: value of query in config
        """
        return get_from_config(
            self.config['DATASETS'][dataset],
            query=query,
            default=default,
            warn=warn
        )


    def get_config(self) -> dict:
        """
        Get complete config with parsed input files
        """
        return self.config


    def get_datasets(self) -> dict:
        """
        Get config for datasets that use the module
        """
        return self.datasets


    def get_wildcards(self, **kwargs) -> [dict, pd.DataFrame]:
        """
        Retrieve wildcard instances as dictionary
        
        :param **kwargs: arguments passed to WildcardParameters.get_wildcards
        :return: dictionary of wildcards that can be applied directly for expanding target files
        """
        return self.parameters.get_wildcards(**kwargs)


    def get_wildcard_names(self):
        return self.parameters.wildcard_names


    def get_output_files(self, pattern: [str, Rule] = None, exclude: list = None) -> list:
        """
        Get output file based on wildcards
        """
        if pattern is None:
            pattern = self.default_output
        return expand(pattern, zip, **self.get_wildcards(exclude=exclude))


    def get_input_file_wildcards(self):
        return self.input_files.get_wildcards()


    def get_input_files(self):
        return self.input_files.get_files()


    def get_input_files_per_dataset(self, dataset):
        return self.input_files.get_files_per_dataset(dataset)


    def get_input_file(self, dataset, file_id):
        return self.input_files.get_file(dataset, file_id)


    def get_parameters(self):
        return self.parameters.get_parameters()


    def get_paramspace(self):
        return self.parameters.get_paramspace()


    def update_inputs(self, dataset: str, input_files: [str, dict]):
        self.config['DATASETS'][dataset]['input'][self.module_name] = input_files
        self.__init__(
            module_name=self.module_name,
            config=self.config,
            parameters=self.parameters_df,
            default_output=self.default_output
        )
