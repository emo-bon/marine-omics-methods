import pytest
import pandas as pd
import matplotlib.pyplot as plt
from momics.plotting import *
from momics.diversity import alpha_diversity_parametrized, beta_diversity_parametrized
from skbio.stats.ordination import pcoa


@pytest.fixture
def sample_metadata():
    """
    Fixture that provides sample metadata for testing.
    """
    metadata = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "factor1": ["A", "B", "C"],
    }
    return pd.DataFrame(metadata)

@pytest.fixture
def sample_pcoa_df():
    """
    Fixture that provides sample PCoA data for testing.
    """
    data = {
        "PC1": [0.1, 0.2, 0.3],
        "PC2": [0.4, 0.5, 0.6],
        "color_by": [1, 2, 3],
    }
    return pd.DataFrame(data)


def test_plot_pcoa_black_no_color(sample_pcoa_df):
    """
    Tests the plot_pcoa_black function without coloring.
    """
    fig = plot_pcoa_black(sample_pcoa_df)

    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"

    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == "PC1", "The x-axis label should be 'PC1'"
    assert ax.get_ylabel() == "PC2", "The y-axis label should be 'PC2'"
    assert ax.get_title() == "PCoA Plot", "The title should be 'PCoA Plot'"

    # Check if the scatter plot is created with black color
    scatter = ax.collections[0]
    assert scatter.get_facecolor()[0, :3].tolist() == [
        0,
        0,
        0,
    ], "The points should be black"


def test_plot_pcoa_black_with_color(sample_pcoa_df):
    """
    Tests the plot_pcoa_black function with coloring.
    """
    fig = plot_pcoa_black(sample_pcoa_df, color_by="color_by")

    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"

    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == "PC1", "The x-axis label should be 'PC1'"
    assert ax.get_ylabel() == "PC2", "The y-axis label should be 'PC2'"
    assert ax.get_title() == "PCoA Plot", "The title should be 'PCoA Plot'"

    # Check if the scatter plot is created with the correct colormap
    scatter = ax.collections[0]
    assert scatter.get_cmap().name == "RdYlGn", "The colormap should be 'RdYlGn'"

    # Check if the colorbar is added
    colorbar = fig.axes[1]
    assert colorbar.get_label() == "color_by", "The colorbar label should be 'color_by'"


@pytest.fixture
def sample_alpha_df():
    """Fixture that provides sample alpha diversity data for testing."""
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "Shannon": [1.0, 2.0, 3.0],
        "factor": ["A", "B", "A"],
    }
    return pd.DataFrame(data)


def test_mpl_alpha_diversity(sample_alpha_df):
    """
    Tests the mpl_alpha_diversity function.
    """
    factor = "factor"
    fig = mpl_alpha_diversity(sample_alpha_df, factor=factor)

    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"

    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Sample", "The x-axis label should be 'Sample'"
    assert (
        ax.get_ylabel() == "Shannon Index"
    ), "The y-axis label should be 'Shannon Index'"
    assert (
        ax.get_title() == f"Shannon Index Grouped by {factor}"
    ), f"The title should be 'Shannon Index Grouped by {factor}'"

    # Check if the bar plot is created with the correct data
    bars = ax.patches
    expected_heights = sample_alpha_df["Shannon"].tolist()
    actual_heights = [bar.get_height() for bar in bars]
    assert (
        actual_heights == expected_heights
    ), f"Expected bar heights {expected_heights}, but got {actual_heights}"

    # Check if the hue is applied correctly
    expected_hues = sample_alpha_df[factor].tolist()
    actual_hues = [bar.get_facecolor() for bar in bars]
    unique_hues = sns.color_palette("coolwarm", len(set(expected_hues)))
    expected_colors = [
        unique_hues[sorted(set(expected_hues)).index(hue)] for hue in expected_hues
    ]
    assert all(
        [
            actual_hue[:3] == expected_color
            for actual_hue, expected_color in zip(actual_hues, expected_colors)
        ]
    ), "The hues are not applied correctly"


def test_mpl_average_per_factor(sample_alpha_df):
    """
    Tests the mpl_average_per_factor function.
    """
    factor = "factor"
    fig = mpl_average_per_factor(sample_alpha_df, factor=factor)
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == factor, f"The x-axis label should be '{factor}'"
    assert ax.get_ylabel() == "Shannon Index", "The y-axis label should be 'Shannon Index'"
    assert ax.get_title() == f"Average Shannon Index Grouped by {factor}", f"The title should be 'Average Shannon Index Grouped by {factor}'"
    
    # Check if the bar plot is created with the correct data
    bars = ax.patches
    expected_heights = sample_alpha_df.groupby(factor)["Shannon"].mean().tolist()
    actual_heights = [bar.get_height() for bar in bars]
    assert actual_heights == expected_heights, f"Expected bar heights {expected_heights}, but got {actual_heights}"
    
    # Check if the hue is applied correctly
    expected_hues = sample_alpha_df[factor].unique().tolist()
    actual_hues = [bar.get_facecolor() for bar in bars]
    unique_hues = sns.color_palette("coolwarm", len(set(expected_hues)))
    expected_colors = [unique_hues[sorted(set(expected_hues)).index(hue)] for hue in expected_hues]
    assert all([actual_hue[:3] == expected_color for actual_hue, expected_color in zip(actual_hues, expected_colors)]), "The hues are not applied correctly"


@pytest.fixture
def sample_tables_dict():
    """
    Fixture that provides a dictionary of sample tables for testing.
    """
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "GO:0001": [10, 0, 5],
        "GO:0002": [20, 0, 5],
        "IPR0001": [30, 0, 5],
        "K0001": [40, 0, 5],
        "PF0001": [50, 0, 5],
        "abundance": [1, 1, 1],  # Add the 'abundance' column
    }
    return {"sample_table": pd.DataFrame(data)}

@pytest.fixture
def sample_metadata():
    """
    Fixture that provides sample metadata for testing.
    """
    metadata = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "factor1": ["A", "B", "C"],
    }
    return pd.DataFrame(metadata)

def test_alpha_plot(sample_tables_dict, sample_metadata):
    """
    Tests the alpha_plot function.
    """
    table_name = "sample_table"
    factor = "factor1"
    fig_pane = alpha_plot(sample_tables_dict, table_name, factor, sample_metadata)
    
    # Check if the result is a panel Matplotlib pane
    assert isinstance(fig_pane, pn.pane.Matplotlib), "The result should be a panel Matplotlib pane"
    
    # Extract the figure from the pane
    fig = fig_pane.object
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Sample", "The x-axis label should be 'Sample'"
    assert ax.get_ylabel() == "Shannon Index", "The y-axis label should be 'Shannon Index'"
    assert ax.get_title() == f"Shannon Index Grouped by {factor}", f"The title should be 'Shannon Index Grouped by {factor}'"
    
    # Check if the bar plot is created with the correct data
    alpha_df = alpha_diversity_parametrized(sample_tables_dict, table_name, sample_metadata)
    bars = ax.patches
    expected_heights = alpha_df["Shannon"].tolist()
    actual_heights = [bar.get_height() for bar in bars]
    assert actual_heights == expected_heights, f"Expected bar heights {expected_heights}, but got {actual_heights}"

    # Check if the hue is applied correctly
    expected_hues = alpha_df[factor].tolist()
    actual_hues = [bar.get_facecolor() for bar in bars]
    unique_hues = sns.color_palette("coolwarm", len(set(expected_hues)))
    expected_colors = [unique_hues[sorted(set(expected_hues)).index(hue)] for hue in expected_hues]
    assert all([actual_hue[:3] == expected_color for actual_hue, expected_color in zip(actual_hues, expected_colors)]), "The hues are not applied correctly"


def test_av_alpha_plot(sample_tables_dict, sample_metadata):
    """
    Tests the av_alpha_plot function.
    """
    table_name = "sample_table"
    factor = "factor1"
    fig_pane = av_alpha_plot(sample_tables_dict, table_name, factor, sample_metadata)
    
    # Check if the result is a panel Matplotlib pane
    assert isinstance(fig_pane, pn.pane.Matplotlib), "The result should be a panel Matplotlib pane"
    
    # Extract the figure from the pane
    fig = fig_pane.object
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == factor, f"The x-axis label should be '{factor}'"
    assert ax.get_ylabel() == "Shannon Index", "The y-axis label should be 'Shannon Index'"
    assert ax.get_title() == f"Average Shannon Index Grouped by {factor}", f"The title should be 'Average Shannon Index Grouped by {factor}'"
    
    # Check if the bar plot is created with the correct data
    alpha_df = alpha_diversity_parametrized(sample_tables_dict, table_name, sample_metadata)
    expected_heights = alpha_df.groupby(factor)["Shannon"].mean().tolist()
    bars = ax.patches
    actual_heights = [bar.get_height() for bar in bars]
    assert actual_heights == expected_heights, f"Expected bar heights {expected_heights}, but got {actual_heights}"
    
    # Check if the hue is applied correctly
    expected_hues = alpha_df[factor].unique().tolist()
    actual_hues = [bar.get_facecolor() for bar in bars]
    unique_hues = sns.color_palette("coolwarm", len(set(expected_hues)))
    expected_colors = [unique_hues[sorted(set(expected_hues)).index(hue)] for hue in expected_hues]
    assert all([actual_hue[:3] == expected_color for actual_hue, expected_color in zip(actual_hues, expected_colors)]), "The hues are not applied correctly"


def test_beta_plot(sample_tables_dict):
    """
    Tests the beta_plot function.
    """
    table_name = "sample_table"
    taxon = "GO:0001"
    fig_pane = beta_plot(sample_tables_dict, table_name, taxon)
    
    # Check if the result is a panel Matplotlib pane
    assert isinstance(fig_pane, pn.pane.Matplotlib), "The result should be a panel Matplotlib pane"
    
    # Extract the figure from the pane
    fig = fig_pane.object
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_title() == f"Beta diversity for {taxon}", f"The title should be 'Beta diversity for {taxon}'"
    
    # Check if the heatmap is created with the correct data
    beta_df = beta_diversity_parametrized(sample_tables_dict[table_name], taxon=taxon, metric="braycurtis")
    heatmap_data = beta_df.to_data_frame()
    heatmap_values = ax.collections[0].get_array().reshape(heatmap_data.shape)
    expected_values = heatmap_data.values
    assert (heatmap_values == expected_values).all(), "The heatmap values are not correct"


def test_beta_plot_pc(sample_tables_dict, sample_metadata):
    """
    Tests the beta_plot_pc function.
    """
    table_name = "sample_table"
    factor = "factor1"
    taxon = "GO:0001"
    fig_pane = beta_plot_pc(sample_tables_dict, sample_metadata, table_name, factor, taxon)
    
    # Check if the result is a panel Matplotlib pane
    assert isinstance(fig_pane, pn.pane.Matplotlib), "The result should be a panel Matplotlib pane"
    
    # Extract the figure from the pane
    fig = fig_pane.object
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct labels and title
    ax = fig.axes[0]
    assert ax.get_xlabel() == "PC1", "The x-axis label should be 'PC1'"
    assert ax.get_ylabel() == "PC2", "The y-axis label should be 'PC2'"
    assert ax.get_title() == "PCoA Plot", "The title should be 'PCoA Plot'"
    
    # Check if the scatter plot is created with the correct data
    beta_df = beta_diversity_parametrized(sample_tables_dict[table_name], taxon=taxon, metric="braycurtis")
    pcoa_result = pcoa(beta_df, method="eigh", number_of_dimensions=3)
    pcoa_df = pd.merge(pcoa_result.samples, sample_metadata, left_index=True, right_on="ref_code", how="inner")
    scatter = ax.collections[0]
    assert all(scatter.get_offsets()[:, 0] == pcoa_df["PC1"]), "The PC1 values are not correct"
    assert all(scatter.get_offsets()[:, 1] == pcoa_df["PC2"]), "The PC2 values are not correct"
    
    # Check if the color is applied correctly
    if factor:
        assert scatter.get_array().tolist() == pcoa_df[factor].tolist(), "The color values are not correct"


@pytest.fixture
def sample_beta_df():
    """
    Fixture that provides sample beta diversity data for testing.
    """
    data = {
        "sample1": [0.0, 0.5, 0.8],
        "sample2": [0.5, 0.0, 0.3],
        "sample3": [0.8, 0.3, 0.0],
    }
    return pd.DataFrame(data, index=["sample1", "sample2", "sample3"])

def test_mpl_plot_heatmap(sample_beta_df):
    """
    Tests the mpl_plot_heatmap function.
    """
    taxon = "GO:0001"
    fig = mpl_plot_heatmap(sample_beta_df, taxon)
    
    # Check if the result is a matplotlib Figure
    assert isinstance(fig, plt.Figure), "The result should be a matplotlib Figure"
    
    # Check if the plot contains the correct title
    ax = fig.axes[0]
    assert ax.get_title() == f"Beta diversity for {taxon}", f"The title should be 'Beta diversity for {taxon}'"
    
    # Check if the heatmap is created with the correct data
    heatmap_values = ax.collections[0].get_array().reshape(sample_beta_df.shape)
    expected_values = sample_beta_df.values
    assert (heatmap_values == expected_values).all(), "The heatmap values are not correct"