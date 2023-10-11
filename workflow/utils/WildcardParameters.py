from pprint import pformat
import numpy as np
import pandas as pd
from snakemake.utils import Paramspace
from snakemake.io import Wildcards

from .config import _get_or_default_from_config
from .misc import unique_dataframe


class WildcardParameters:
    
    def __init__(
        self,
        module_name: str,
        parameters: pd.DataFrame,
        input_file_wildcards: dict,
        dataset_config: dict,
        default_config: dict,
        wildcard_names: list,
        config_params: list = None,
        explode_by: [str, list] = None,
        paramspace_kwargs: dict = None,
    ):
        """
        :param module_name: name of module
        :param parameters: dataframe with parameters for each dataset
        :param input_file_wildcards: dictionary with input file wildcards for each dataset
        :param dataset_config: dictionary with dataset configuration for module
        :param default_config: dictionary with default configuration for module
        :param wildcard_names: list of wildcard names for expanding rules
        :param config_params: list of parameters that a module should consider as wildcards, order and length must match wildcard_names, by default will take wildcard_names
        :param explode_by: column(s) to explode wildcard_names extracted from config by
        """
        self.module_name = module_name
        
        # determine wildcards
        self.dataset_config = dataset_config
        self.default_config = default_config
        self.set_wildcards(
            config_params=config_params,
            wildcard_names=wildcard_names,
            explode_by=explode_by,
        )
        mandatory_wildcards = ['dataset', 'file_id']
        for wildcard in mandatory_wildcards:
            while wildcard in wildcard_names:
                wildcard_names.remove(wildcard)
        self.wildcard_names = mandatory_wildcards + wildcard_names
        
        if parameters is not None:
            self.parameters_df = parameters
            if parameters.shape[1] > 0:
                self.wildcards_df = self.wildcards_df.merge(self.parameters_df, how='left')
        else:
            self.parameters_df = pd.DataFrame()
        
        # add input file wildcards
        self.wildcards_df = self.wildcards_df.merge(
            pd.DataFrame(input_file_wildcards, dtype='object'),
            on='dataset',
            how='left',
        )
        
        # set paramspace
        if paramspace_kwargs is None:
            paramspace_kwargs = {}
        self.paramspace_kwargs = paramspace_kwargs
        self.paramspace = Paramspace(
            self.wildcards_df[self.wildcard_names],
            **paramspace_kwargs
        )


    def subset_by_query(
        self,
        query_dict: [dict, Wildcards],
        columns: list = None
    ) -> pd.DataFrame:
        """
        Helper function to subset self.wildcards_df by query dictionary
        
        :param query_dict: Mapping of column and value to subset by
        :param columns: columns to subset wildcards_df by
        :return: subset of wildcards_df
        """
        if not query_dict and not columns:
            return self.wildcards_df
        
        if columns is None:
            columns = self.wildcards_df.columns
        
        df = self.wildcards_df.astype(str).copy()
        # subset by query
        for key, value in query_dict.items():
            if key not in df.columns:
                continue
            df = df[df[key] == value]
        
        # subset by columns
        df = unique_dataframe(df[columns]).reset_index(drop=True)
        assert df.shape[0] > 0, f'No wildcard combination found in wildcards df {query_dict}\n{self.wildcards_df}'
        
        return df


    def set_wildcards(
        self,
        config_params: list = None,
        wildcard_names: list = None,
        explode_by: list = None,
        config_entries: list = None,
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
        :param config_entries: list of entries to subset the config by, otherwise use all keys
        """
        if not wildcard_names:
            wildcard_names = []
        
        if not config_params:  # len(self.dataset_config) == 0 or 
            config_params = wildcard_names
        
        if len(config_params) != len(wildcard_names):
            raise ValueError('config_params and wildcard_names must be of same length')
            
        if config_entries is None:
            config_entries = self.dataset_config.keys()

        # collect entries for dataframe
        defaults = self.default_config.get(self.module_name)
        records = [
            (
                key,
                *[
                    _get_or_default_from_config(
                        config=self.dataset_config.get(key),
                        defaults=defaults,
                        key=self.module_name,
                        value=param,
                        update=True,
                        warn=warn,
                    )
                    for param in config_params
                ]
            )
            for key in config_entries
        ]
        df = pd.DataFrame.from_records(records, columns=[*['dataset']+wildcard_names])
        if explode_by is not None:
            explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
            for column in explode_by:
                df = df.explode(column)
        self.wildcards_df = df.reset_index(drop=True)


    def update(self, wildcards_df=None, parameters_df=None, wildcard_names=None, **kwargs):
        if wildcards_df is not None:
            self.wildcards_df = wildcards_df
        if parameters_df is not None:
            self.parameters_df = parameters_df
        if wildcard_names is not None:
            self.wildcard_names = wildcard_names
        self.paramspace_kwargs = kwargs
        self.paramspace = Paramspace(
            self.wildcards_df[self.wildcard_names],
            **kwargs
        )


    def get_wildcards(
        self,
        subset_dict: [dict, Wildcards] = None,
        exclude: [list, str] = None,
        wildcard_names: list = None,
        all_params: bool = False,
        as_df: bool = False,
    ) -> [dict, pd.DataFrame]:
        """
        Retrieve wildcard instances as dictionary
        
        :param exclude: list of wildcard names to exclude
        :return: dictionary of wildcards that can be applied directly for expanding target files
        """
        if exclude is None:
            exclude = []
        elif isinstance(exclude, str):
            exclude = [exclude]
        if subset_dict is None:
            subset_dict = {}
        if wildcard_names is None:
            wildcard_names = self.wildcards_df.columns if all_params else self.wildcard_names
        
        df = self.subset_by_query(
            query_dict=subset_dict,
            columns=[w for w in wildcard_names if w not in exclude]
        )
        return df if as_df else df.to_dict('list')


    def get_parameters(self) -> pd.DataFrame:
        return self.parameters_df


    def get_paramspace(
        self,
        wildcard_names: list = None,
        exclude: list = None,
        **kwargs
    ) -> Paramspace:
        """
        :param default: whether to return the default paramspace (all wildcards)
        :param wildcard_names: list of wildcard names to subset the paramspace by
        :param kwargs: additional arguments for snakemake.utils.Paramspace
        :return: snakemake.utils.Paramspace object
        """
        if not wildcard_names and not exclude and not kwargs:
            return self.paramspace
        if exclude is None:
            exclude = []
        if wildcard_names is None:
            wildcard_names = [w for w in self.wildcard_names if w not in exclude]
        if not kwargs:
            kwargs = self.paramspace_kwargs
        return Paramspace(
            self.wildcards_df[wildcard_names],
            **kwargs
        )


    def get_from_parameters(
        self,
        query_dict: [dict, Wildcards],
        parameter_key: str,
        wildcards_sub: [list, None] = None,
    ):
        """
        Get entries from parameters dataframe

        :param query_dict: dictionary with column (must be present in parameters_df) to value mapping
        :param parameter_key: key of parameter
        :param wildcards_sub: list of wildcards used for subsetting the parameters
        :return: single parameter value or list of parameters as specified by column
        """
        if wildcards_sub is None:
            wildcards_sub = self.wildcards_df.columns.tolist()
        
        try:
            assert parameter_key in self.wildcards_df.columns, f'"{parameter_key}" not in wildcards_df.columns'
            assert np.all([key in wildcards_sub for key in query_dict.keys()]), 'Not all query keys in wildcards_sub'
        except AssertionError as e:
            raise AssertionError(
                f'{e} for parameter_key={parameter_key}, wildcards_sub={wildcards_sub}'
                f'\nquery_dict:\n{pformat(query_dict)}'
                f'\n{self.wildcards_df}'
            ) from e

        wildcards_sub = list(set(wildcards_sub + [parameter_key]))
        params_sub = self.subset_by_query(
            query_dict={k: v for k, v in query_dict.items() if k in wildcards_sub},
            columns=[parameter_key]
        )
        assert params_sub.shape[0] <= 1, f'More than 1 row after subsetting\n{params_sub}'
        
        return params_sub[parameter_key].tolist()[0]
