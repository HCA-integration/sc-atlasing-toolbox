from utils.ModuleConfig import ModuleConfig
from utils.WildcardParameters import WildcardParameters


class ClusteringConfig(ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        super().__init__(**kwargs)
        
        # parse hierarchical clustering parameters
        def parse_levels(x):
            if isinstance(x, dict):
                max_level = max(x.keys())
            elif isinstance(x, list):
                max_level = len(x)
            elif isinstance(x, int):
                max_level = x
            elif isinstance(x, str):
                max_level = int(x)
            else:
                max_level = 1
            return list(range(1, max_level+1))

        wildcards_df = self.parameters.wildcards_df
        wildcards_df['level'] = wildcards_df['hierarchy'].apply(parse_levels)
        wildcards_df = wildcards_df.explode('level').reset_index(drop=True)
        wildcards_df['level'] = wildcards_df['level'].astype(str)
        
        self.update_parameters(
            wildcards_df=wildcards_df,
            wildcard_names=self.parameters.wildcard_names+['level'],
        )
