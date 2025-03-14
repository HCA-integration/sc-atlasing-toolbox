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
        
        def get_max_level(x):
            if isinstance(x, dict):
                return max(x.keys())
            elif isinstance(x, list):
                return len(x)
            elif isinstance(x, int):
                return x
            elif isinstance(x, str):
                return int(x)
            return 1
        
        def parse_levels(x):
            return list(range(1, get_max_level(x)+1))

        def parse_kwargs(x):
            levels = parse_levels(x)
            if isinstance(x, dict):
                kwargs = [x.get(level, {}) for level in levels]
                return [
                    {'resolution': kw} if not isinstance(kw, dict) else kw
                    for kw in kwargs
                ]
            return [{}] * len(levels)


        wildcards_df = self.parameters.wildcards_df

        # check if all algorithms in config are configured
        available_algorihtms = ['leiden', 'louvain']
        for algorithm in wildcards_df['algorithm'].unique():
            assert algorithm in available_algorihtms, f'Algorithm not available: "{algorithm}"'
        
        # parse hierarchical clustering parameters
        wildcards_df['level'] = wildcards_df['hierarchy'].apply(parse_levels)
        wildcards_df['kwargs'] = wildcards_df['hierarchy'].apply(parse_kwargs)

        wildcards_df = wildcards_df.explode(['level', 'kwargs']).reset_index(drop=True)
        wildcards_df['level'] = wildcards_df['level'].astype(str)
        
        self.update_parameters(
            wildcards_df=wildcards_df,
            wildcard_names=self.parameters.wildcard_names#+['level'],
        )
