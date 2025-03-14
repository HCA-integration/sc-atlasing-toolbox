import pandas as pd

from utils.ModuleConfig import ModuleConfig
from utils.misc import unique_dataframe, expand_dict_and_serialize


class PreprocessingConfig(ModuleConfig):
    
    HVG_PARAMS = [
        'n_top_genes',
        'min_disp',
        'max_disp',
        'min_mean',
        'max_mean',
        'span',
        'n_bins',
        'flavor',
        'batch_key'
    ]
    
    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        super().__init__(**kwargs)
        self.set_extra_hvgs()

    
    def set_extra_hvgs(self):
        wildcards_df = self.parameters.wildcards_df
        records = []
        for dataset, dataset_dict in self.parameters.dataset_config.items():
            ehvg_config = wildcards_df.query('dataset == @dataset')['extra_hvgs'].iloc[0]
            if ehvg_config is None or not isinstance(ehvg_config.get('overwrite_args'), dict):
                records.append((dataset, 'None', None))
            else:
                _dict = {
                    k: v for k, v in ehvg_config['overwrite_args'].items()
                    if k in self.HVG_PARAMS
                }
                records.extend(
                    (dataset, *rec) for rec in
                    expand_dict_and_serialize(_dict)
                )
            
        extra_hvgs_df = pd.DataFrame(
            records,
            columns=['dataset', 'overwrite_args', 'overwrite_args_dict']
        )
        
        wildcards_df = unique_dataframe(
            wildcards_df.merge(extra_hvgs_df, on='dataset', how='left')
        )
        self.update_parameters(wildcards_df=wildcards_df)
                
            