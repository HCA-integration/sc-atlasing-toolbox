from pathlib import Path
import numpy as np
import pandas as pd
from snakemake.rules import Rule

from integration.IntegrationConfig import IntegrationConfig
from utils.ModuleConfig import ModuleConfig
from utils.WildcardParameters import WildcardParameters
from utils.config import _get_or_default_from_config
from utils.misc import expand_dict_and_serialize, unique_dataframe


class SampleRepresentationConfig(IntegrationConfig, ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        # parse output type info
        parameters = kwargs.get('parameters')
        if isinstance(parameters, str):
            parameters = pd.read_table(parameters)
        if isinstance(parameters, pd.DataFrame):
            parameters['input_type'] = parameters['input_type'].str.split(',')
            parameters = parameters.explode('input_type')
            kwargs['parameters'] = parameters
        
        ModuleConfig.__init__(self, **kwargs, write_output_files=False)
        
        # set hyperparameters
        self.set_hyperparams()

        # prune entries depending on input type
        wildcards_df = self.parameters.wildcards_df
        # set use_rep to embedding or counts depending on input type
        wildcards_df['use_rep'] = np.where(
            wildcards_df['input_type'] == 'embed',
            wildcards_df['use_rep'],
            wildcards_df['raw_counts'] # TODO: case for normalized counts
        )
        # set var_mask only for counts representations
        wildcards_df['var_mask'] = np.where(
            (
                (wildcards_df['use_rep'] == 'X') |
                wildcards_df['use_rep'].str.startswith(('layers/', 'raw/'))
            ),
            wildcards_df['var_mask'],
            'None'
        )
        # drop duplicates
        wildcards_df = wildcards_df.drop_duplicates(
            ['dataset', 'file_id', 'method', 'input_type', 'use_rep', 'var_mask', 'hyperparams']
        ).reset_index(drop=True)
        
        # set default paramspace arguments
        if kwargs.get('paramspace_kwargs') is None:
            paramspace_kwargs = dict(
                filename_params=['method', 'hyperparams', 'var_mask', 'use_rep'],
                filename_sep='--',
            )
        else:
            paramspace_kwargs = {}
        self.update_parameters(
            wildcards_df=wildcards_df,
            wildcard_names=self.parameters.wildcard_names + ['hyperparams'],
            **paramspace_kwargs,
        )

        # write output files
        self.write_output_files()
