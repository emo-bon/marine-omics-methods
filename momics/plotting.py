import panel as pn
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

from skbio.stats.ordination import pcoa

from .diversity import (
    alpha_diversity_parametrized,
    beta_diversity_parametrized,
)


def plot_pcoa_black(pcoa_df: pd.DataFrame, color_by: str = None) -> plt.Figure:
    """
    Plots a PCoA plot with optional coloring.

    Args:
        pcoa_df (pd.DataFrame): A DataFrame containing PCoA results.
        color_by (str, optional): The column name to color the points by. Defaults to None.

    Returns:
        plt.Figure: The PCoA plot.
    """
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    if color_by is not None:
        scatter = plt.scatter(
            pcoa_df["PC1"],
            pcoa_df["PC2"],
            c=pcoa_df[color_by],
            cmap="RdYlGn",
            edgecolor="k",
        )
        plt.colorbar(scatter, label=color_by)
    else:
        plt.scatter(pcoa_df["PC1"], pcoa_df["PC2"], color="black")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCoA Plot")
    plt.tight_layout()
    plt.close(plot)
    return plot


def mpl_alpha_diversity(alpha_df: pd.DataFrame, factor: str = None) -> plt.Figure:
    """
    Plots the Shannon index grouped by a factor.

    Args:
        alpha_df (pd.DataFrame): A DataFrame containing alpha diversity results.
        factor (str, optional): The column name to group by. Defaults to None.

    Returns:
        plt.Figure: The Shannon index plot.
    """
    alpha_df = alpha_df.sort_values(by=factor)
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    sns.barplot(
        data=alpha_df,
        x="ref_code",
        y="Shannon",
        hue=factor,
        palette="coolwarm",
    )

    ax.set_title(f"Shannon Index Grouped by {factor}")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Shannon Index")

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.close(plot)
    return plot


def mpl_average_per_factor(df: pd.DataFrame, factor: str = None) -> plt.Figure:
    """
    Plots the average Shannon index grouped by a factor.

    Args:
        df (pd.DataFrame): A DataFrame containing alpha diversity results.
        factor (str, optional): The column name to group by. Defaults to None.

    Returns:
        plt.Figure: The average Shannon index plot.
    """
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    sns.barplot(
        data=df,
        x=factor,
        y="Shannon",
        hue=factor,
        capsize=0.1,
        palette="coolwarm",
    )

    ax.set_title(f"Average Shannon Index Grouped by {factor}")
    ax.set_xlabel(factor)
    ax.set_ylabel("Shannon Index")

    plt.tight_layout()
    plt.close(plot)
    return plot


##################
# Plot for panel #
##################
# Alpha diversity
def alpha_plot(
    tables_dict: Dict[str, pd.DataFrame],
    table_name: str,
    factor: str,
    metadata: pd.DataFrame,
) -> pn.pane.Matplotlib:
    """
    Creates an alpha diversity plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        factor (str): The column name to group by.
        metadata (pd.DataFrame): A DataFrame containing metadata.

    Returns:
        pn.pane.Matplotlib: A Matplotlib pane containing the alpha diversity plot.
    """
    alpha = alpha_diversity_parametrized(tables_dict, table_name, metadata)
    fig = pn.pane.Matplotlib(
        mpl_alpha_diversity(alpha, factor=factor),
        sizing_mode="stretch_both",
        name="Alpha div",
    )
    return fig


def av_alpha_plot(
    tables_dict: Dict[str, pd.DataFrame],
    table_name: str,
    factor: str,
    metadata: pd.DataFrame,
) -> pn.pane.Matplotlib:
    """
    Creates an average alpha diversity plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        factor (str): The column name to group by.
        metadata (pd.DataFrame): A DataFrame containing metadata.

    Returns:
        pn.pane.Matplotlib: A Matplotlib pane containing the average alpha diversity plot.
    """
    alpha = alpha_diversity_parametrized(tables_dict, table_name, metadata)
    fig = pn.pane.Matplotlib(
        mpl_average_per_factor(alpha, factor=factor),
        sizing_mode="stretch_both",
        name="AV Alpha div",
    )
    return fig


def beta_plot(
    tables_dict: Dict[str, pd.DataFrame], table_name: str, taxon: str = "ncbi_tax_id"
) -> pn.pane.Matplotlib:
    """
    Creates a beta diversity heatmap plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        taxon (str, optional): The taxon level for beta diversity calculation. Defaults to "ncbi_tax_id".

    Returns:
        pn.pane.Matplotlib: A Matplotlib pane containing the beta diversity heatmap plot.
    """
    beta = beta_diversity_parametrized(
        tables_dict[table_name], taxon=taxon, metric="braycurtis"
    )

    fig = pn.pane.Matplotlib(
        mpl_plot_heatmap(beta.to_data_frame(), taxon=taxon),
        sizing_mode="stretch_both",
        name="Beta div",
    )
    return fig


def beta_plot_pc(
    tables_dict: Dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    table_name: str,
    factor: str,
    taxon: str = "ncbi_tax_id",
) -> pn.pane.Matplotlib:
    """
    Creates a beta diversity PCoA plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        metadata (pd.DataFrame): A DataFrame containing metadata.
        table_name (str): The name of the table to process.
        factor (str): The column name to color the points by.
        taxon (str, optional): The taxon level for beta diversity calculation. Defaults to "ncbi_tax_id".

    Returns:
        pn.pane.Matplotlib: A Matplotlib pane containing the beta diversity PCoA plot.
    """
    beta = beta_diversity_parametrized(
        tables_dict[table_name], taxon=taxon, metric="braycurtis"
    )
    pcoa_result = pcoa(beta, method="eigh", number_of_dimensions=3)
    pcoa_df = pd.merge(
        pcoa_result.samples,
        metadata,
        left_index=True,
        right_on="ref_code",
        how="inner",
    )

    fig = pn.pane.Matplotlib(
        plot_pcoa_black(pcoa_df, color_by=factor),
        sizing_mode="stretch_both",
        name="Beta PCoA",
    )
    return fig


def mpl_plot_heatmap(df: pd.DataFrame, taxon: str) -> plt.Figure:
    """
    Creates a heatmap plot for beta diversity.

    Args:
        df (pd.DataFrame): A DataFrame containing beta diversity distances.
        taxon (str): The taxon level for beta diversity calculation.

    Returns:
        plt.Figure: The heatmap plot.
    """
    plot = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    ax = plot.add_subplot(111)

    sns.heatmap(df, cmap="viridis")
    plt.title(f"Beta diversity for {taxon}")
    plt.tight_layout()
    plt.close(plot)
    return plot
