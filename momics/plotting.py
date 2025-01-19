import panel as pn
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

from .diversity import (
    diversity_input,
    alpha_diversity_parametrized,
    beta_diversity_parametrized,
)


# Plot the PCoA with optional coloring
# TODO: color_by does not work for categorical data
def plot_pcoa_black(pcoa_df, color_by=None):
    
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    if color_by is not None:
        scatter = plt.scatter(
            pcoa_df["PC1"], pcoa_df["PC2"],
            c=pcoa_df[color_by],
            cmap="RdYlGn", edgecolor="k"
        )
        # plt.colorbar(scatter, label=color_by.name)
        plt.colorbar(scatter, label=color_by)
    else:
        plt.scatter(pcoa_df["PC1"], pcoa_df["PC2"], color="black")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCoA Plot")
    plt.tight_layout()
    plt.close(plot)
    return plot


def mpl_alpha_diversity(alpha_df: pd.DataFrame, factor: str = None):
    alpha_df = alpha_df.sort_values(by=factor)
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)
    
    sns.barplot(
        data=alpha_df,
        x='ref_code',
        y='Shannon',
        hue=factor,
        palette="coolwarm",
        )

    ax.set_title(f"Shannon Index Grouped by {factor}")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Shannon Index")

    # TODO: check length for conditional rotation
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.close(plot)
    return plot


def mpl_average_per_factor(df: pd.DataFrame, factor: str = None):
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)
    
    sns.barplot(
        data=df,
        x=factor,
        y='Shannon',
        hue=factor,
        capsize=0.1,
        palette="coolwarm",
        )
    
    ax.set_title(f"Average Shannon Index Grouped by {factor}")
    ax.set_xlabel(factor)
    ax.set_ylabel("Average Shannon Index")
    plt.tight_layout()
    plt.close(plot)
    return plot

##################
# Plot for panel #
##################
# Alpha diversity
def alpha_plot(table_list, table_name, factor, metadata):
    alpha = alpha_diversity_parametrized(table_list, table_name, metadata)
    fig = pn.pane.Matplotlib(
        mpl_alpha_diversity(alpha, factor=factor),
                            sizing_mode="stretch_both",
                            name="Alpha div",
                            )
    return fig


def av_alpha_plot(table_list, table_name, factor, metadata):
    alpha = alpha_diversity_parametrized(table_list, table_name, metadata)
    fig = pn.pane.Matplotlib(
        mpl_average_per_factor(alpha, factor=factor),
        sizing_mode="stretch_both",
        name="AV Alpha div",
        )
    return fig


def beta_plot(table_list, table_name, taxon: str = "ncbi_tax_id"):
    beta = beta_diversity_parametrized(table_list[table_name],
                                       taxon=taxon,
                                       metric="braycurtis")

    fig = pn.pane.Matplotlib(
        mpl_plot_heatmap(beta.to_data_frame(), taxon=taxon),
                         sizing_mode="stretch_both",
                         name="Beta div",
                         )
    return fig


def beta_plot_pc(table_list, metadata, table_name, factor, taxon: str = "ncbi_tax_id"):
    beta = beta_diversity_parametrized(table_list[table_name],
                                       taxon=taxon,
                                       metric="braycurtis")
    pcoa_result = pcoa(beta, method='eigh', number_of_dimensions=3)
    pcoa_df = pd.merge(
        pcoa_result.samples,
        metadata,
        left_index=True,
        right_on="ref_code",
        how="inner",
    )

    # df_beta = pd.merge(beta.to_data_frame(),
    #                    table_list["metadata"],
    #                    left_index=True, right_on='ref_code')
    fig = pn.pane.Matplotlib(
        plot_pcoa_black(pcoa_df, color_by=factor),
        sizing_mode="stretch_both",
        name="Beta PCoA",
        )
    return fig


def mpl_plot_heatmap(df, taxon):
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    sns.heatmap(df, cmap="viridis")
    plt.title(f"Beta diversity for {taxon}")
    plt.tight_layout()
    plt.close(plot)
    return plot
