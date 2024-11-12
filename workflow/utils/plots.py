import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def plot_stacked_bar(
    df,
    category_key,
    covariate_key,
    count_type='count',
    palette='tab20',
    title_suffix='',
    fig_width=8,
    fig_height_min=5,
):
    colors = palette
    if palette is not None:
        colors = plt.get_cmap(palette).colors
    
    # determine number of cells per group and category
    df = df[[category_key, covariate_key]].copy()
    df[category_key] = df[category_key].astype(str)
    counts_df = df[[category_key, covariate_key]].value_counts(dropna=False).to_frame()
    counts_df['total_count'] = counts_df.groupby(category_key, observed=True).transform('sum')
    
    if count_type == 'proportion':
        counts_df['proportion'] = counts_df['count'] / counts_df['total_count']
    counts_df = counts_df.reset_index(drop=False)
    
    # determine sort order of categories
    total_counts = counts_df.drop_duplicates(category_key) \
        .set_index(category_key)['total_count'].sort_values(ascending=True)
    counts_df[category_key] = pd.Categorical(
        counts_df[category_key].astype(str),
        categories=total_counts.index,
        ordered=True
    )

    # create plot
    num_categories = counts_df[category_key].nunique()
    fig_height = max(fig_height_min, num_categories * 0.3)
    plt.figure(figsize=(fig_width, fig_height))
    
    ax = counts_df.pivot(
        index=category_key,
        columns=covariate_key,
        values=count_type
    ).fillna(0) \
    .sort_index() \
    .plot.barh(
        stacked=True,
        color=colors,
        ax=plt.gca()
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    longest_label = 0
    for y_position, (category, total_value) in enumerate(total_counts.items()):
        x_position = 1 if count_type == 'proportion' else total_value
        label = f'{total_value:.2e}' if total_value >= 1000 else str(int(total_value))
        # label = '1.23e+06'  # for debugging
        ax.text(
            x_position,
            y_position,
            label,
            ha='left',
            va='center',
            color='black'
        )
        longest_label = max(longest_label, len(label))
        
    if df[covariate_key].nunique() < 30:
        # adjust figure size to avoid overlapping labels with legend
        label_length_factor = longest_label * 0.06
        plt.gcf().set_size_inches(fig_width + label_length_factor, fig_height)
        ax.legend(loc='center left', bbox_to_anchor=(1 + label_length_factor / fig_width, 0.5))
    else:
        ax.legend_ = None
    
    title = f'{count_type.capitalize()} of {covariate_key} by {category_key}{title_suffix}'
    ax.set_title(title)
    
    return plt.gcf()


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
    
    return plt.gcf()
