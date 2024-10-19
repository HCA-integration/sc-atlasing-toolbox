from matplotlib import pyplot as plt
import seaborn as sns


def plot_stacked_bar(df, category_key, covariate_key, count_type='count', palette='tab20', title_suffix=''):
    colors = palette
    if palette is not None:
        colors = plt.get_cmap(palette).colors
    df = df[[category_key, covariate_key]].copy()
    df[category_key] = df[category_key].astype(str)
    counts_df = df[[category_key, covariate_key]].value_counts().to_frame()
    if count_type == 'proportion':
        counts_df = counts_df / counts_df.groupby(category_key, observed=True).transform('sum')
        counts_df.columns = ['proportion']
    counts_df = counts_df.reset_index(drop=False)
    
    num_categories = counts_df[category_key].nunique()
    plt.figure(figsize=(8, max(5, num_categories * 0.3)))
    ax = counts_df.pivot(
        index=category_key,
        columns=covariate_key,
        values=count_type
    ).fillna(0) \
    .plot.barh(
        stacked=True,
        color=colors,
        ax=plt.gca()
    )
        
    if df[covariate_key].nunique() < 30:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        ax.legend_ = None
    title = f'{count_type.capitalize()} of {covariate_key} by {category_key}{title_suffix}'
    ax.set_title(title)


def plot_violin(df, category_key, covariate_key):
    clusters = sorted(df[category_key].unique())
    grouped_data = [
        df[df[category_key] == cluster][covariate_key]
        for cluster in clusters
    ]
    num_categories = len(grouped_data)
    plt.figure(figsize=(8, max(5, num_categories * 0.3)))
    plt.violinplot(
        grouped_data,
        showmeans=False,
        showmedians=True,
        vert=False,
    )
    plt.yticks(ticks=range(1, len(grouped_data) + 1), labels=clusters)
    plt.title(f'Distribution of {covariate_key} by {category_key}')
    plt.xlabel(covariate_key)
    plt.ylabel(category_key)
