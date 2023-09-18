from typing import Union
from pathlib import Path
import pandas as pd
import hashlib
from snakemake.io import expand
from snakemake.utils import Paramspace
from snakemake.rules import Rule

from .config import _get_or_default_from_config
from .misc import unique_dataframe


class ModuleConfig:

    # module_name = ''
    # config = {}
    # datasets = {}
    # wildcard_names = []
    # input_files = {}
    # input_file_wildcards = dict(dataset=[], file_id=[], file_name=[])
    # wildcards_df = pd.DataFrame()
    # parameters_df = pd.DataFrame()
    # paramspace = Paramspace()
    # default_output = None
    # out_dir = Path()
    # image_dir = Path()


    def __init__(self, module_name: str, config: dict, parameters: pd.DataFrame = None, default_output: [str, Rule] = None):
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
        
        self.input_files = {}
        self.input_file_wildcards = dict(dataset=[], file_id=[], file_name=[])
        for dataset in self.datasets:
            self.set_defaults_per_dataset(dataset)
            self.set_input_file_per_dataset(dataset)


        # determine wildcards
        self.set_wildcards()
        self.wildcard_names = ['dataset', 'file_id']
        
        # get parameters
        if parameters is not None:
            self.parameters_df = parameters
            if parameters.shape[1] > 0:
                self.wildcards_df = self.wildcards_df.merge(self.parameters_df, how='left')
        else:
            self.parameters_df = pd.DataFrame()
        
        # add input file wildcards
        self.wildcards_df = self.wildcards_df.merge(
            pd.DataFrame(self.input_file_wildcards, dtype='object'),
            on='dataset',
            how='left',
        )
        
        self.paramspace = Paramspace(self.wildcards_df[self.wildcard_names])


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


    @staticmethod
    def parse_input_files(input_files: str, digest_size: int = 5) -> bool:
        if isinstance(input_files, str):
            input_files = [input_files]
        
        if isinstance(input_files, list):
            input_files = {
                hashlib.blake2b(file.encode('utf-8'), digest_size=digest_size).hexdigest(): file
                for file in input_files
            }
        
        if not isinstance(input_files, dict):
            raise ValueError(f'input_files must be a list or dict, but is {type(input_files)}')
        
        return input_files


    def set_input_file_per_dataset(self, dataset: str, digest_size: int = 5):
        """Get input files for a given module and dataset
        This function maps an input file with its unique identifier, if no file ID is specified
        
        :param dataset: dataset key in config['DATASETS']
        """
        input_files = self.get_for_dataset(dataset, query=['input', self.module_name])
        input_files = self.parse_input_files(input_files, digest_size=digest_size)
        self.input_files[dataset] = input_files
        
        # Reshape the file ID to dataset mapping to a lists of wildcards
        for file_id, file in input_files.items():
            self.input_file_wildcards['dataset'].append(dataset)
            self.input_file_wildcards['file_id'].append(file_id)
            self.input_file_wildcards['file_name'].append(file)
        
        # write to file
        pd.DataFrame(self.input_files).to_csv(
            self.out_dir / 'input_files.tsv',
            sep='\t',
            index=False
        )


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
        # start at top level
        value = self.config['DATASETS'][dataset]

        # walk down query
        for q in query:
            if warn and q not in value:
                warnings.warn(f'key {q} not found in config for query {query}, returning default')
            value = value.get(q, default)
        return value


    def set_wildcards(
        self,
        config_params: list = None,
        wildcard_names: list = None,
        explode_by: list = None,
        config_keys: list = None,
        warn: bool = False,
    ):
        """
        Collect wildcards and parameters from an configuration instance (e.g. dataset) for a given module.
        This function assumes that the keys of the given config keys are the different instances that contain specific parameters for different modules.

        :param config_params: List of parameters for each config entry of a module
            e.g. ['integration', 'label', 'batch']
        :param wildcard_names: names of wildcards to be extracted.
            Must map to config keys, and prepended by a wildcard name for the config entries
            e.g. ['dataset', 'method', 'label', 'batch']
        :param explode_by: column to explode by, expecting list entry for that column
        :param config_keys: list of entries to subset the config by, otherwise use all keys
        """
        if len(self.datasets) == 0 or not config_params:
            config_params = []
        
        if not wildcard_names:
            wildcard_names = config_params
        
        if len(config_params) != len(wildcard_names):
            raise ValueError('config_params and wildcard_names must be of same length')
            
        if config_keys is None:
            config_keys = self.datasets.keys()

        # collect entries for dataframe
        records = [
            (
                key,
                *[
                    _get_or_default_from_config(
                        config=self.datasets[key],
                        defaults=self.config['defaults'].get(self.module_name),
                        key=self.module_name,
                        value=param,
                        update=True,
                        warn=warn,
                    )
                    for param in config_params
                ]
            )
            for key in config_keys
        ]
        df = pd.DataFrame.from_records(records, columns=[*['dataset']+wildcard_names])
        if explode_by is not None:
            explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
            for column in explode_by:
                df = df.explode(column)
        self.wildcards_df = df.reset_index(drop=True)


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


    def get_wildcards(self, exclude: list = None) -> dict:
        """
        Retrieve wildcard instances for this module
        
        :param exclude: list of wildcard names to exclude
        :return: dictionary of wildcards that can be applied directly for expanding target files
        """
        if exclude is None:
            exclude = []
        wildcard_names = [wildcard for wildcard in self.wildcard_names if wildcard not in exclude]
        return unique_dataframe(self.wildcards_df[wildcard_names]).to_dict('list')


    def get_output_files(self, pattern: [str, Rule] = None, exclude: list = None) -> list:
        if pattern is None:
            pattern = self.default_output
        if exclude is None:
            exclude = []
        return expand(pattern, zip, **self.get_wildcards(exclude=exclude))


    def get_input_file_wildcards(self):
        return self.input_file_wildcards


    def get_input_files(self):
        return self.input_files


    def get_input_files_per_dataset(self, dataset):
        return self.input_files[dataset]


    def get_input_file(self, dataset, file_id):
        return self.get_input_files_per_dataset(dataset)[file_id]


    def get_parameters(self):
        return self.parameters_df


    def get_paramspace(self):
        return self.paramspace


    def update_inputs(self, dataset: str, input_files: [str, dict]):
        self.config['DATASETS'][dataset]['input'][self.module_name] = input_files
        # print(self.input_files)
        # print(self.input_file_wildcards)
        self.input_files = {}
        self.input_file_wildcards = dict(dataset=[], file_id=[], file_name=[])
        self.__init__(self.module_name, self.config, self.parameters_df, self.default_output)
