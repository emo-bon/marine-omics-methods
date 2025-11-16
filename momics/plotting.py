"""

Constants
---------
- PLOT_FACE_COLOR : str
    The face color for the plot.

TODO:
- Returns should be plt.figure and not pn.pane.Matplotlib, as already implemented for beta_plot_pc() function.
"""

import logging
from typing import List, Tuple, Dict, Union
from textwrap import fill
import numpy as np
import panel as pn
import pandas as pd
import networkx as nx

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import holoviews as hv
import hvplot.pandas  # noqa

from skbio.diversity import beta_diversity
from bokeh.models import CategoricalColorMapper, ContinuousColorMapper, LogColorMapper
from bokeh.palettes import Category20, viridis

from skbio.stats.ordination import pcoa
from .diversity import (
    alpha_diversity_parametrized,
    beta_diversity_parametrized,
)
from .utils import (
    check_index_names,
)

PLOT_FACE_COLOR = "#e6e6e6"
MARKER_SIZE = 16

# logger setup
FORMAT = "%(levelname)s | %(name)s | %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


##########
# HVplot #
##########
def hvplot_heatmap(
    df: pd.DataFrame, taxon: str, norm: bool = False
) -> hv.element.HeatMap:
    """
    Creates a heatmap plot for beta diversity using hvplot.

    Args:
        df (pd.DataFrame): DataFrame containing beta diversity distances.
        taxon (str): The taxon level for beta diversity calculation.
        norm (bool): Whether to normalize the data.

    Returns:
        hv.element.HeatMap: A heatmap plot of beta diversity.
    """
    # Create the heatmap using hvplot
    heatmap = df.hvplot.heatmap(
        cmap="viridis",
        colorbar=True,
        xlabel="Sample",
        ylabel="Sample",
        title=f"Beta Diversity ({taxon})",
    ).opts(
        xticks=0,  # remove xticks labels
        yticks=0,  # remove yticks labels
        show_legend=False,  # hide legend
    )
    if norm:
        heatmap.opts(clim=(0, 1.0))  # Set color limits for normalization
    else:
        heatmap.opts(
            clim=(df.min().min(), df.max().max())
        )  # Set color limits based on data
    return heatmap


def hvplot_alpha_diversity(alpha: pd.DataFrame, factor: str) -> hv.element.Bars:
    """
    Creates a horizontal bar plot for alpha diversity using hvplot.

    Args:
        alpha (pd.DataFrame): DataFrame containing alpha diversity data.
        factor (str): The column name to group by.

    Returns:
        hv.element.Bars: A horizontal bar plot of alpha diversity.
    """
    # Define the color mapper using Bokeh's CategoricalColorMapper
    if 2 < len(alpha[factor].unique()) <= 20:
        pal = Category20[
            len(alpha[factor].unique())
        ]  # Use the correct number of colors
    else:
        pal = viridis(len(alpha[factor].unique()))

    color_mapper = CategoricalColorMapper(
        factors=alpha[factor]
        .unique()
        .tolist(),  # Unique categories in the factor column
        palette=pal,
    )

    # Create the horizontal bar plot using hvplot
    fig = alpha.hvplot.barh(
        y="Shannon",
        xlabel="Sample",
        ylabel="Shannon Index",
        title=f"Alpha Diversity ({factor})",
        color=factor,  # Use the factor column for coloring
    ).opts(
        yticks=0,  # remove yticks labels
        xaxis="top",
        cmap=color_mapper.palette,  # Apply the color mapper's palette
        legend_position="top_right",  # Adjust legend position
        tools=["hover"],  # Add hover tool for interactivity
        backend_opts={"plot.toolbar.autohide": True},
    )
    return fig


def hvplot_average_per_factor(alpha: pd.DataFrame, factor: str) -> hv.element.Bars:
    """
    Creates a horizontal bar plot for alpha diversity using hvplot.

    Args:
        alpha (pd.DataFrame): DataFrame containing alpha diversity data.
        factor (str): The column name to group by.

    Returns:
        hv.element.Bars: A horizontal bar plot of alpha diversity.
    """
    # Define the color mapper using Bokeh's CategoricalColorMapper
    if len(alpha[factor].unique()) <= 20:
        pal = Category20[
            len(alpha[factor].unique())
        ]  # Use the correct number of colors
    else:
        pal = viridis(len(alpha[factor].unique()))

    color_mapper = CategoricalColorMapper(
        factors=alpha[factor]
        .unique()
        .tolist(),  # Unique categories in the factor column
        palette=pal,
    )

    # Create the horizontal bar plot using hvplot
    fig = alpha.hvplot.bar(
        x=factor,
        y="Shannon",
        xlabel=factor,
        ylabel="Shannon Index",
        title=f"Average Shannon Index Grouped by {factor}",
        color=factor,  # Use the factor column for coloring
        hover_cols=[alpha.index.name],
    ).opts(
        cmap=color_mapper.palette,  # Apply the color mapper's palette
        show_legend=False,
        # legend_position="top_right",  # Adjust legend position
        # tools=["hover"],  # Add hover tool for interactivity
        backend_opts={"plot.toolbar.autohide": True},
    )
    return fig


def hvplot_plot_pcoa_black(
    pcoa_df: pd.DataFrame,
    color_by: str = None,
    explained_variance: Tuple[float, float] = None,
    **kwargs,
) -> hv.element.Scatter:
    """
    Plots a PCoA plot with optional coloring using hvplot.

    Args:
        pcoa_df (pd.DataFrame): A DataFrame containing PCoA results.
        color_by (str, optional): The column name to color the points by. Defaults to None.

    Returns:
        hv.element.Scatter: The PCoA plot.
    """
    log_scale = kwargs.get('log_scale', False)
    palette = kwargs.get('palette', "Turbo256")
    index_name = pcoa_df.index.name if pcoa_df.index.name else "sample"
    pcoa_df = pcoa_df.reset_index(names=index_name)  # Ensure index is a column for hvplot

    if color_by is None:
        # No coloring specified, use black
        fig = pcoa_df.hvplot.scatter(
            x="PC1",
            y="PC2",
            color="black",
            hover_cols=[index_name, "PC1", "PC2"],
        )
        valid_perc = 100.0
        title = "PCoA (no coloring applied)"
        color_palette = None
    else:
        valid_perc = pcoa_df[color_by].count() / len(pcoa_df[color_by]) * 100
        # Handle logarithmic scaling for continuous data
        if log_scale and pcoa_df[color_by].dtype != "object":
            # Handle zeros and negative values for log scaling
            color_data = pcoa_df[color_by].copy()
            color_data = color_data.where(color_data > 0, 1e-2)

            
            # Update the DataFrame with processed data
            pcoa_df[f'{color_by}_log'] = color_data
            color_column = f'{color_by}_log'
            
            # Create logarithmic color mapper
            if color_data.min() > 0:
                color_mapper = LogColorMapper(
                    palette=palette,
                    low=color_data.min(),
                    high=color_data.max(),
                )
            else:
                # Fallback to linear if log scaling fails
                log_scale = False
                color_column = color_by
                logger.info(f"Warning: Cannot use log scale due to non-positive values. Falling back to linear scale.")
        else:
            color_column = color_by
        if not log_scale:  # Original logic for non-log scaling
            if 2 < len(pcoa_df[color_by].unique()) <= 20:
                if pcoa_df[color_by].dtype == "object":
                    pal = Category20[
                        len(pcoa_df[color_by].unique())
                    ]  # Use the correct number of colors
                    color_mapper = CategoricalColorMapper(
                        factors=pcoa_df[color_by]
                        .unique()
                        .tolist(),  # Unique categories in the factor column
                        palette=pal,
                    )
                else:
                    color_mapper = ContinuousColorMapper(
                        palette="Turbo256",
                        low=pcoa_df[color_by].min(),
                        high=pcoa_df[color_by].max(),
                    )
            else:
                if pcoa_df[color_by].dtype == "object":
                    pal = viridis(len(pcoa_df[color_by].unique()))
                    color_mapper = CategoricalColorMapper(
                        factors=pcoa_df[color_by]
                        .unique()
                        .tolist(),  # Unique categories in the factor column
                        palette=pal,
                    )
                else:
                    color_mapper = ContinuousColorMapper(
                        palette="Turbo256",
                        low=pcoa_df[color_by].min(),
                        high=pcoa_df[color_by].max(),
                    )

        if pcoa_df[color_by].count() >= 0:
            # Create the scatter plot using hvplot
            hvplot_kwargs = {
                "x": "PC1",
                "y": "PC2",
                "color": color_column if log_scale else color_by,
                "hover_cols": [index_name, "PC1", "PC2"],
            }
            
            # Add log-specific hvplot options
            if log_scale and 'color_mapper' in locals():
                hvplot_kwargs["logz"] = True
                
            fig = pcoa_df.hvplot.scatter(**hvplot_kwargs)
        else:
            fig = pcoa_df.hvplot.scatter(
                x="PC1",
                y="PC2",
                color="black",  # Use black for coloring
                hover_cols=[index_name, "PC1", "PC2"],
            )
        
        # Update title to indicate log scale
        scale_info = " (log scale)" if log_scale else ""
        title = f"PCoA colored by {color_by}{scale_info}, valid values: ({valid_perc:.2f}%)"

        if 'color_mapper' in locals():
            color_palette = color_mapper.palette
        else:
            color_palette = None

    if explained_variance:
        var_perc = explained_variance[0] * 100, explained_variance[1] * 100
        fig = fig.opts(
            xlabel=f"PC1 ({var_perc[0]:.2f}%)",
            ylabel=f"PC2 ({var_perc[1]:.2f}%)",
        )
    else:
        fig = fig.opts(
            xlabel="PC1",
            ylabel="PC2",
        )
    
    assert "PC1" in pcoa_df.columns, f"Missing 'PC1' column in PCoA DataFrame"
    assert "PC2" in pcoa_df.columns, f"Missing 'PC2' column in PCoA DataFrame"

    opts = {
        "title": title,
        "size": MARKER_SIZE,
        "fill_alpha": 0.5,
        "show_legend": False,
        "backend_opts": {"plot.toolbar.autohide": True},
    }
    
    # Add log-specific options
    if log_scale:
        opts["logz"] = True
    
    if color_palette is not None:
        opts["cmap"] = color_palette
    
    fig = fig.opts(**opts)
    return fig


##############
# Matplotlib #
##############
def plot_pcoa_black(pcoa_df: pd.DataFrame, color_by: str = None) -> plt.Figure:
    """
    Plots a PCoA plot with optional coloring.

    Args:
        pcoa_df (pd.DataFrame): A DataFrame containing PCoA results.
        color_by (str, optional): The column name to color the points by. Defaults to None.

    Returns:
        plt.Figure: The PCoA plot.
    """
    flag_massage = False
    fig = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    fig.patch.set_facecolor(PLOT_FACE_COLOR)
    ax = fig.add_subplot(111)

    if color_by is not None:
        labels = fold_legend_labels_from_series(pcoa_df[color_by], 35)
        # BETA created now only for numerical
        if pcoa_df[color_by].dtype == "object":
            sns.scatterplot(
                data=pcoa_df,
                x="PC1",
                y="PC2",
                hue=color_by,
                palette="coolwarm",
                edgecolor="black",
            )
            ax = change_legend_labels(ax, labels)
        else:
            flag_massage = True
            perc = pcoa_df[color_by].count() / len(pcoa_df[color_by]) * 100
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
    if flag_massage:
        ax.set_title(f"PCoA Plot with valid {color_by} values: ({perc:.2f}%)")
    else:
        ax.set_title("PCoA Plot")
    plt.tight_layout()
    plt.close(fig)
    return fig


def mpl_alpha_diversity(alpha_df: pd.DataFrame, factor: str = None) -> plt.Figure:
    """Plots the Shannon index grouped by a factor.

    Args:
        alpha_df (pd.DataFrame): A DataFrame containing alpha diversity results.
        factor (str, optional): The column name to group by. Defaults to None.

    Returns:
        plt.Figure: The Shannon index plot.
    """
    fig = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    fig.patch.set_facecolor(PLOT_FACE_COLOR)
    labels = fold_legend_labels_from_series(alpha_df[factor], 35)

    ax = fig.add_subplot(111)
    sns.barplot(
        data=alpha_df,
        x=alpha_df.index,
        y="Shannon",
        hue=factor,
        palette="coolwarm",
    )

    # check axes and find which is have legend
    ax = change_legend_labels(ax, labels)
    ax.tick_params(axis="x", which="major", labelsize=np.log(4e5 / len(alpha_df)))

    ax.set_title(f"Shannon Index Grouped by {factor}")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Shannon Index")

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.close(fig)
    return fig


def mpl_average_per_factor(df: pd.DataFrame, factor: str = None) -> plt.Figure:
    """Plots the average Shannon index grouped by a factor.

    Args:
        df (pd.DataFrame): A DataFrame containing alpha diversity results.
        factor (str, optional): The column name to group by. Defaults to None.

    Returns:
        plt.Figure: The average Shannon index plot.
    """
    fig = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    fig.patch.set_facecolor(PLOT_FACE_COLOR)
    ax = fig.add_subplot(111)

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
    ax = cut_xaxis_labels(ax, 15)

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.close(fig)
    return fig


## gecco analysis ##
####################
def mpl_bgcs_violin(df: pd.DataFrame, normalize: bool = False) -> plt.Figure:
    fig = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    fig.patch.set_facecolor(PLOT_FACE_COLOR)
    ax = fig.add_subplot(111)

    df["type"] = df["type"].fillna("Unknown")

    sns.swarmplot(
        data=df, x="type", y="average_p", hue="max_p", size=7, palette="coolwarm", ax=ax
    )
    sns.violinplot(
        data=df,
        x="type",
        y="average_p",
        inner=None,
        fill=False,
        color="black",
        ax=ax,
    )
    if normalize:
        ax.set_ylim(-0.1, 1.1)

    ax.set_title(f"Probabilities of identified BGCs by type")
    ax.set_xlabel("BGC type")
    ax.set_ylabel("Average probability")
    ax = cut_xaxis_labels(ax, 15)

    plt.tight_layout()
    plt.close(fig)
    return fig


def hvplot_bgcs_violin(df: pd.DataFrame, normalize: bool = False) -> hv.Overlay:
    """
    Creates a violin plot for BGC probabilities by type using hvplot.

    Args:
        df (pd.DataFrame): A DataFrame containing BGC data with columns 'type', 'average_p', and 'max_p'.
        normalize (bool): Whether to normalize the y-axis to the range [0, 1].

    Returns:
        hv.Overlay: An overlay of swarm and violin plots.
    """
    df["type"] = df["type"].fillna("Unknown")

    # Create the swarm plot
    swarm = df.hvplot.scatter(
        x="type",
        y="average_p",
        c="max_p",
        size=MARKER_SIZE,
        cmap="coolwarm",
        colorbar=True,
        title="Probabilities of identified BGCs by type",
        xlabel="BGC type",
        ylabel="Average probability",
        tools=["hover"],
    )

    # Create the violin plot
    violin = df.hvplot.violin(
        # x="type",
        y="average_p",
        by="type",
    ).opts(
        violin_fill_color=PLOT_FACE_COLOR,
        violin_alpha=0.5,
        violin_line_width=0.2,
        violin_visible=True,
    )

    # Combine the swarm and violin plots
    combined = violin * swarm

    # Apply normalization if needed
    if normalize:
        combined = combined.opts(yaxis="log", ylim=(-0.1, 1.1))

    return combined


def plot_domain_abundance(
    filtered_domains: pd.Series, abundance_min: int
) -> hv.element.Bars:
    """
    Plot the histogram of the number of pfam domains from the feature table using hvplot.

    Args:
        filtered_domains (pd.Series): A Series containing domain names as the index and their abundances as values.
        abundance_min (int): The minimum abundance threshold for domains to be included in the plot.

    Returns:
        hv.element.Bars: A horizontal bar plot of domain abundances.
    """
    plot = filtered_domains.hvplot.bar(
        xlabel="Domain",
        ylabel="Abundance",
        title=f"Domain Abundance (at least {abundance_min})",
        height=600,
        width=1400,
        cmap="viridis",
        tools=["hover"],
    ).opts(
        xticks=0,  # Remove xtick labels if too many
        show_legend=False,
    )
    return plot


def plot_tsne(X_embedded: np.ndarray, kmeans) -> hv.element.Scatter:
    """
    Plot the t-SNE embedding of the clusters using hvplot.

    Args:
        X_embedded (np.ndarray): The t-SNE embedded coordinates.
        kmeans: The fitted KMeans object containing cluster labels.

    Returns:
        hv.element.Scatter: The t-SNE plot of domain clusters.
    """
    # Create a DataFrame for hvplot
    tsne_df = pd.DataFrame(
        {
            "t-SNE 1": X_embedded[:, 0],
            "t-SNE 2": X_embedded[:, 1],
            "Cluster": kmeans.labels_,
        }
    )

    # Create the scatter plot using hvplot
    plot = tsne_df.hvplot.scatter(
        x="t-SNE 1",
        y="t-SNE 2",
        c="Cluster",
        cmap="viridis",
        colorbar=True,
        size=MARKER_SIZE,
        title="t-SNE of Domain Clusters",
        xlabel="t-SNE 1",
        ylabel="t-SNE 2",
        tools=["hover"],
    )
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
    order: str = "factor",  # or values
    backend: str = "hvplot",  # Options: "matplotlib" or "hvplot"
) -> Union[pn.pane.Matplotlib, pn.pane.HoloViews]:
    """
    Creates an alpha diversity plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        factor (str): The column name to group by.
        metadata (pd.DataFrame): A DataFrame containing metadata.
        order (str): The order of sorting the data. Can be "factor" or "value".
        backend (str): The plotting backend to use. Can be "matplotlib" or "hvplot".

    Returns:
        Union[pn.pane.Matplotlib, pn.pane.HoloViews]: A pane containing the alpha diversity plot.
    """
    alpha = alpha_diversity_parametrized(tables_dict, table_name, metadata)
    hash_sort = {"factor": factor, "values": "Shannon"}
    alpha = alpha.sort_values(by=hash_sort[order])

    if backend == "matplotlib":
        fig = pn.pane.Matplotlib(
            mpl_alpha_diversity(alpha, factor=factor),
            sizing_mode="stretch_both",
            name="Alpha div",
        )
    elif backend == "hvplot":
        fig = pn.pane.HoloViews(
            hvplot_alpha_diversity(alpha, factor=factor),
            name="Alpha div",
            width=900,
            height=1500,
        )
    else:
        raise ValueError(f"Unknown backend: {backend}")
    return fig


def av_alpha_plot(
    tables_dict: Dict[str, pd.DataFrame],
    table_name: str,
    factor: str,
    metadata: pd.DataFrame,
    order: str = "factor",  # or values
    backend: str = "hvplot",  # Options: "matplotlib" or "hvplot"
) -> Union[pn.pane.Matplotlib, pn.pane.HoloViews]:
    """
    Creates an average alpha diversity plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        factor (str): The column name to group by.
        metadata (pd.DataFrame): A DataFrame containing metadata.

    Returns:
        Union[pn.pane.Matplotlib, pn.pane.HoloViews]: A pane containing the average alpha diversity plot.
    """
    alpha = alpha_diversity_parametrized(tables_dict, table_name, metadata)
    # TODO: this will not work, because it gets grouped in mpl_average_per_factor I think
    hash_sort = {"factor": factor, "values": "Shannon"}
    alpha = alpha.sort_values(by=hash_sort[order])

    if backend == "matplotlib":
        fig = pn.pane.Matplotlib(
            mpl_average_per_factor(alpha, factor=factor),
            sizing_mode="stretch_both",
            name="AV Alpha div",
        )
    elif backend == "hvplot":
        fig = pn.pane.HoloViews(
            hvplot_average_per_factor(alpha, factor=factor),
            sizing_mode="stretch_both",
            name="AV Alpha div",
        )
    else:
        raise ValueError(f"Unknown backend: {backend}")
    return fig


def beta_plot(
    tables_dict: Dict[str, pd.DataFrame],
    table_name: str,
    norm: bool,
    taxon: str = "ncbi_tax_id",
    backend: str = "hvplot",  # Options: "matplotlib" or "hvplot"
) -> Union[pn.pane.Matplotlib, pn.pane.HoloViews]:
    """
    Creates a beta diversity heatmap plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.
        taxon (str, optional): The taxon level for beta diversity calculation. Defaults to "ncbi_tax_id".
        norm (bool): Whether to normalize the data.
        backend (str): The plotting backend to use. Can be "matplotlib" or "hvplot".

    Returns:
        Union[pn.pane.Matplotlib, pn.pane.HoloViews]: A pane containing the beta diversity heatmap plot.
    """
    beta = beta_diversity_parametrized(
        tables_dict[table_name], taxon=taxon, metric="braycurtis"
    )

    if backend == "matplotlib":
        fig = pn.pane.Matplotlib(
            mpl_plot_heatmap(beta.to_data_frame(), taxon=taxon, norm=norm),
            sizing_mode="stretch_both",
            name="Beta div",
        )
    elif backend == "hvplot":
        fig = pn.pane.HoloViews(
            hvplot_heatmap(beta.to_data_frame(), taxon=taxon, norm=norm),
            sizing_mode="stretch_both",
            name="Beta div",
        )
    else:
        raise ValueError(f"Unknown backend: {backend}")
    return fig


def beta_plot_pc(
    tables_dict: Dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    table_name: str,
    factor: str,
    taxon: str = "ncbi_tax_id",
) -> Tuple[hv.element.Scatter, Tuple[float, float]]:
    """
    Creates a beta diversity PCoA plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        metadata (pd.DataFrame): A DataFrame containing metadata.
        table_name (str): The name of the table to process.
        factor (str): The column name to color the points by.
        taxon (str, optional): The taxon level for beta diversity calculation. Defaults to "ncbi_tax_id".

    Returns:
        Tuple[hv.element.Scatter, Tuple[float, float]]: A tuple containing the beta diversity PCoA plot and the explained variance for PC1 and PC2.
    """
    beta = beta_diversity_parametrized(
        tables_dict[table_name], taxon=taxon, metric="braycurtis"
    )
    pcoa_result = pcoa(beta, method="eigh")
    explained_variance = (
        pcoa_result.proportion_explained[0],
        pcoa_result.proportion_explained[1],
    )
    if not set(pcoa_result.samples.index) == set(metadata.index):
        raise ValueError("Metadata index name does not match PCoA result.")

    pcoa_df = pd.merge(
        pcoa_result.samples,
        metadata,
        left_index=True,
        right_index=True,
        how="inner",
    )
    return (
        hvplot_plot_pcoa_black(
            pcoa_df, color_by=factor, explained_variance=explained_variance
        ),
        explained_variance,
    )


def beta_plot_pc_granular(
    filtered_data: pd.DataFrame,
    metadata: pd.DataFrame,
    factor: str,
) -> Tuple[hv.element.Scatter, Tuple[float, float]]:
    """
    Creates a beta diversity PCoA plot.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        metadata (pd.DataFrame): A DataFrame containing metadata.
        table_name (str): The name of the table to process.
        factor (str): The column name to color the points by.
        taxon (str, optional): The taxon level for beta diversity calculation. Defaults to "ncbi_tax_id".

    Returns:
        Tuple[plt.figure, float]: A tuple containing the beta diversity PCoA plot and the explained variance.
    """
    # beta = beta_diversity("braycurtis", filtered_data.iloc[:, 1:].T)
    beta = beta_diversity("braycurtis", filtered_data.T)
    pcoa_result = pcoa(beta, method="eigh")  # , number_of_dimensions=3)
    explained_variance = (
        pcoa_result.proportion_explained[0],
        pcoa_result.proportion_explained[1],
    )

    # Check if metadata index matches PCoA result
    if not set(pcoa_result.samples.index) == set(metadata.index):
        logger.debug(
            f"Metadata index name does not match PCoA result. {set(pcoa_result.samples.index) ^ set(metadata.index)}"
        )
    pcoa_df = pd.merge(
        pcoa_result.samples,
        metadata,
        left_index=True,
        right_index=True,
        how="inner",
    )

    return (
        hvplot_plot_pcoa_black(
            pcoa_df, color_by=factor, explained_variance=explained_variance
        ),
        explained_variance,
    )


def mpl_plot_heatmap(df: pd.DataFrame, taxon: str, norm=False) -> plt.Figure:
    """
    Creates a heatmap plot for beta diversity.

    Args:
        df (pd.DataFrame): A DataFrame containing beta diversity distances.
        taxon (str): The taxon level for beta diversity calculation.
        norm (bool): Whether to normalize the data.

    Returns:
        plt.Figure: The heatmap plot.
    """
    fig = plt.figure(figsize=(10, 6), facecolor=(0, 0, 0, 0))
    fig.patch.set_facecolor(PLOT_FACE_COLOR)
    _ = fig.add_subplot(111)
    if norm:
        sns.heatmap(df, vmin=0, vmax=1.0, cmap="viridis")
    else:
        sns.heatmap(df, cmap="viridis")
    plt.title(f"Beta diversity for {taxon}")
    plt.tight_layout()
    plt.close(fig)
    return fig


####################
# Helper functions #
####################
def fold_legend_labels_from_series(df: pd.Series, max_len: int = 30) -> List[str]:
    """Folds a list of labels to a maximum length from a Series.

    Args:
        df (pd.Series): The series to extract unique labels.
        max_len (int, optional): The maximum length of a label. Defaults to 30.

    Returns:
        List[str]: The folded list of labels.
    """
    return [
        fill(x, max_len) if isinstance(x, str) and len(x) > max_len else str(x)
        for x in df.unique()
    ]


def change_legend_labels(ax: plt.axis, labels: List[str]) -> plt.axis:
    """Changes the labels of a legend on a given matplotlib axis.

        ax (plt.axis): The matplotlib axis object whose legend labels need to be changed.
        labels (List[str]): A list of new labels to be set for the legend.

    Returns:
        plt.axis: The matplotlib axis object with updated legend labels.
    """
    leg = ax.get_legend()
    for t, l in zip(leg.texts, labels):
        t.set_text(l)
    return ax


def cut_xaxis_labels(ax: plt.axis, n: int = 15) -> plt.axis:
    """Changes the x-tick labels by cutting them short.

    Args:
        ax: The axes to change the x-axis of.
        n: cutoff for max number of characters for the xtick label.

    Returns:
        plt.axis: The axes with the new x-tick labels.
    """
    ticks = ax.get_xticklabels()
    new_ticks = []
    for tick in ticks:
        tick.set_text(tick.get_text()[: min(len(tick.get_text()), n)])
        new_ticks.append(tick)
    ax.set_xticklabels(new_ticks)
    return ax


############
## Plotly ##
############
def get_sankey(df, cat_cols=[], value_cols="", title="Sankey Diagram"):
    # Colors
    colorPalette = [
        "rgba(31, 119, 180, 0.8)",
        "rgba(255, 127, 14, 0.8)",
        "rgba(44, 160, 44, 0.8)",
        "rgba(214, 39, 40, 0.8)",
        "rgba(148, 103, 189, 0.8)",
        "rgba(140, 86, 75, 0.8)",
        "rgba(227, 119, 194, 0.8)",
        "rgba(127, 127, 127, 0.8)",
    ]
    labelList = []
    colorNumList = []

    for catCol in cat_cols:
        labelListTemp = list(set(df[catCol].values))
        colorNumList.append(len(labelListTemp))
        labelList = labelList + labelListTemp

    # remove duplicates from labelList
    labelList = list(dict.fromkeys(labelList))

    # define colors based on number of levels
    colorList = []
    for idx, colorNum in enumerate(colorNumList):
        colorList = colorList + [colorPalette[idx]] * colorNum

    # transform df into a source-target pair
    for i in range(len(cat_cols) - 1):
        if i == 0:
            sourceTargetDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
            sourceTargetDf.columns = ["source", "target", "count"]
        else:
            tempDf = df[[cat_cols[i], cat_cols[i + 1], value_cols]]
            tempDf.columns = ["source", "target", "count"]
            sourceTargetDf = pd.concat([sourceTargetDf, tempDf])
        sourceTargetDf = (
            sourceTargetDf.groupby(["source", "target"])
            .agg({"count": "sum"})
            .reset_index()
        )

    # add index for source-target pair
    sourceTargetDf["sourceID"] = sourceTargetDf["source"].apply(
        lambda x: labelList.index(x)
    )
    sourceTargetDf["targetID"] = sourceTargetDf["target"].apply(
        lambda x: labelList.index(x)
    )

    # creating data for the sankey diagram
    data = dict(
        type="sankey",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labelList,
            color=colorList,
        ),
        link=dict(
            source=sourceTargetDf["sourceID"],
            target=sourceTargetDf["targetID"],
            value=sourceTargetDf["count"],
        ),
    )

    # override gray link colors with 'source' colors
    opacity = 0.4
    # change 'magenta' to its 'rgba' value to add opacity
    data["node"]["color"] = [
        "rgba(255,0,255, 0.8)" if color == "magenta" else color
        for color in data["node"]["color"]
    ]
    data["link"]["color"] = [
        data["node"]["color"][src].replace("0.8", str(opacity))
        for src in data["link"]["source"]
    ]

    fig = go.Figure(
        data=[
            go.Sankey(
                # Define nodes
                node=dict(
                    pad=15,
                    thickness=15,
                    line=dict(color="black", width=0.5),
                    label=data["node"]["label"],
                    color=data["node"]["color"],
                ),
                # Add links
                link=dict(
                    source=data["link"]["source"],
                    target=data["link"]["target"],
                    value=data["link"]["value"],
                    color=data["link"]["color"],
                ),
            )
        ]
    )

    fig.update_layout(title_text=title, font_size=10)

    return fig


def plot_network(network_results, association_data, alpha=0.5):
    _, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, factor in zip(axes, list(association_data.keys())):
        G = network_results[factor]["graph"]
        colors = nx.get_edge_attributes(G, "color")
        pos = nx.spring_layout(G, k=0.2, iterations=50, seed=42)
        nx.draw_networkx_nodes(
            G, pos, ax=ax, alpha=alpha, node_color="grey", node_size=17
        )
        nx.draw_networkx_edges(
            G, pos, ax=ax, alpha=alpha - 0.3, edge_color=list(colors.values())
        )
        ax.set_title(factor)
        ax.axis("off")

    plt.tight_layout()
    plt.show()
