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
        input_file_wildcards: pd.DataFrame,
        dataset_config: dict,
        default_config: dict,
        wildcard_names: list,
    ):
        self.module_name = module_name
        
        # determine wildcards
        self.dataset_config = dataset_config
        self.default_config = default_config
        self.set_wildcards()
        self.wildcard_names = wildcard_names
        
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
        
        self.paramspace = Paramspace(self.wildcards_df[self.wildcard_names])


    @staticmethod
    def subset_by_wildcards(
        df: pd.DataFrame,
        wildcards: [dict, Wildcards]
    ) -> pd.DataFrame:
        """
        Helper function to subset df by wildcards
        
        :param df: dataframe with wildcard names in column and wildcard values in rows
        :param wildcards: wildcards object from Snakemake
        :return: subset dataframe
        """
        if not wildcards:
            return df
        wildcards = {k: wildcards[k] for k in wildcards.keys() if k in df.columns}
        query = ' and '.join([f'{k} == "{v}"' for k, v in wildcards.items()])
        df = df.query(query)
        if df.shape[0] == 0:
            raise ValueError(f'no wildcard combination found in wildcards df {query}')
        return df


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
        if len(self.dataset_config) == 0 or not config_params:
            config_params = []
        
        if not wildcard_names:
            wildcard_names = config_params
        
        if len(config_params) != len(wildcard_names):
            raise ValueError('config_params and wildcard_names must be of same length')
            
        if config_keys is None:
            config_keys = self.dataset_config.keys()

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
            for key in config_keys
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
        self.paramspace = Paramspace(
            self.wildcards_df[self.wildcard_names],
            **kwargs
        )


    def get_wildcards(
        self,
        exclude: list = None,
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
        wildcard_names = self.wildcards_df.columns if all_params else self.wildcard_names
        wildcard_names = [wildcard for wildcard in wildcard_names if wildcard not in exclude]
        df = unique_dataframe(self.wildcards_df[wildcard_names])
        return df if as_df else df.to_dict('list')


    def get_parameters(self) -> pd.DataFrame:
        return self.parameters_df


    def get_paramspace(self) -> Paramspace:
        return self.paramspace


    def get_from_parameters(
        self,
        query_dict: dict,
        parameter_key: str,
        wildcards_sub: [list, None] = None
    ):
        """
        Get entries from parameters dataframe

        :param query_dict: dictionary with column (must be present in parameters_df) to value mapping
        :param parameter_key: key of parameter
        :param wildcards_sub: list of wildcards used for subsetting the parameters
        :return: single parameter value or list of parameters as specified by column
        """
        try:
            assert column in self.parameters_df.columns

            if wildcards_sub is None:
                wildcards_sub = self.parameters_df.columns
            assert np.all([key in query_dict.keys() for key in query_dict.keys()])

            # quickfix: ignore hyperparameters
            # if 'hyperparams' in wildcards_sub:
            #     wildcards_sub.remove('hyperparams')

            query_dict = {k: v for k, v in query_dict.items() if k in wildcards_sub}
            params_sub = self.subset_by_wildcards(self.parameters_df, query_dict)
            columns = list(set(wildcards_sub + [column]))
            params_sub = unique_dataframe(params_sub[columns])
        except Exception as e:
            raise ValueError(
                f'Error for wildcards={wildcards}, column={column}, wildcards_keys={wildcards_keys}'
                f'\n{parameters_df}\nError message: {e}'
            ) from e
        
        try:
            assert params_sub.shape[0] == 1
        except AssertionError as e:
            raise ValueError(f'More than 1 row after subsetting\n{params_sub}') from e
        
        param = params_sub[column].tolist()
        return param[0] if len(param) == 1 else param[wildcards_keys]
