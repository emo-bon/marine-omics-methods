import pandas as pd
import pytest
from momics.plotting import *


def sample_tables_dict(name="sample_table", add_abundance=False):
    """Fixture that provides a dictionary of sample tables for testing.
    These need to start with 'GO', 'IPR', 'K', 'PF', which are hardcoded in
    the calculate_alpha_diversity function for instance.
    TODO: that is probably not that great, should raise warning for other columns
    Put in a separate test?
    """
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "GO:0001": [10, 0, 5],
        "GO:0002": [20, 0, 5],
        "IPR0001": [30, 0, 5],
        "K0001": [40, 0, 5],
        "PF0001": [50, 0, 5],
    }
    if add_abundance:
        data["abundance"] = [50, 0, 5]

    out = pd.DataFrame(data)
    out = out.set_index("ref_code")
    return {name: out}


def sample_data(name="sample_table", add_abundance=False):
    return sample_tables_dict(name, add_abundance)[name]


def add_column(df, name):
    """Add a column to a DataFrame."""
    df[name] = [str(k) for k in range(len(df))]
    return df


@pytest.fixture
def sample_factors():
    """Fixture that provides sample factors for testing."""
    factors = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "factor1": ["A", "B", "C"],
    }
    out = pd.DataFrame(factors)
    out = out.set_index("ref_code")

    return out


@pytest.fixture
def sample_metadata():
    """
    Fixture that provides sample metadata for testing.
    """
    metadata = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "factor1": ["A", "B", "C"],
    }
    metadata = pd.DataFrame(metadata)
    metadata = metadata.set_index("ref_code")
    return metadata


# @pytest.fixture
# def sample_tables_dict():
#     """
#     Fixture that provides a dictionary of sample tables for testing.
#     """
#     data = {
#         "ref_code": ["sample1", "sample2", "sample3"],
#         "GO:0001": [10, 0, 5],
#         "GO:0002": [20, 0, 5],
#         "IPR0001": [30, 0, 5],
#         "K0001": [40, 0, 5],
#         "PF0001": [50, 0, 5],
#         "abundance": [1, 1, 1],  # Add the 'abundance' column
#     }
#     return {"sample_table": pd.DataFrame(data)}


@pytest.fixture
def sample_alpha_df():
    """Fixture that provides sample alpha diversity data for testing."""
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "Shannon": [1.0, 2.0, 3.0],
        "factor": ["A", "B", "A"],
    }
    return pd.DataFrame(data)


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
    assert (
        ax.get_title() == f"Beta diversity for {taxon}"
    ), f"The title should be 'Beta diversity for {taxon}'"

    # Check if the heatmap is created with the correct data
    heatmap_values = ax.collections[0].get_array().reshape(sample_beta_df.shape)
    expected_values = sample_beta_df.values
    assert (
        heatmap_values == expected_values
    ).all(), "The heatmap values are not correct"
