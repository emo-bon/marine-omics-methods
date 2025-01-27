import pytest
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity

from momics.diversity import *


@pytest.fixture
def sample_data():
    """Fixture that provides sample data for testing."""
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "GO:0001": [10, 0, 5],
        "GO:0002": [20, 0, 5],
        "IPR0001": [30, 0, 5],
        "K0001": [40, 0, 5],
        "abundance": [50, 0, 5],
        # "ncbi_tax_id": ["taxon1", "taxon2", "taxon3"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_tables_dict():
    """Fixture that provides a dictionary of sample tables for testing."""
    data = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "GO:0001": [10, 0, 5],
        "GO:0002": [20, 0, 5],
        "IPR0001": [30, 0, 5],
        "K0001": [40, 0, 5],
        "PF0001": [50, 0, 5],
    }
    return {"sample_table": pd.DataFrame(data)}


@pytest.fixture
def sample_factors():
    """Fixture that provides sample factors for testing."""
    factors = {
        "ref_code": ["sample1", "sample2", "sample3"],
        "factor1": ["A", "B", "C"],
    }

    return pd.DataFrame(factors)


def test_shannon_index(sample_data):
    """Tests the shannon_index function."""
    # Test with a row containing species abundances
    row = sample_data.iloc[0]
    result = shannon_index(row[1:])
    expected = -sum(
        (x / row[1:].sum()) * np.log(x / row[1:].sum()) for x in row[1:] if x > 0
    )
    assert np.isclose(result, expected), f"Expected {expected}, but got {result}"

    # Test with a row[1:] containing zero abundances
    row = sample_data.iloc[1]
    result = shannon_index(row[1:])
    assert np.isnan(result), f"Expected NaN, but got {result}"

    # Test with a row[1:] containing equal abundances
    row = sample_data.iloc[2]
    result = shannon_index(row[1:])
    expected = -sum(
        (x / row[1:].sum()) * np.log(x / row[1:].sum()) for x in row[1:] if x > 0
    )
    assert np.isclose(result, expected), f"Expected {expected}, but got {result}"


def test_calculate_alpha_diversity(sample_data, sample_factors):
    """Tests the calculate_alpha_diversity function."""
    result = calculate_alpha_diversity(sample_data, sample_factors)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected columns
    expected_columns = ["ref_code", "Shannon", "factor1"]
    assert all(
        col in result.columns for col in expected_columns
    ), f"Expected columns {expected_columns}, but got {result.columns.tolist()}"

    # Check if the Shannon index values are calculated correctly
    expected_shannon = sample_data.apply(lambda row: shannon_index(row[1:]), axis=1)
    assert all(
        result["Shannon"].round(6) == expected_shannon.round(6)
    ), "The Shannon index values are not calculated correctly"

    # Check if the factors are merged correctly
    assert all(
        result["factor1"] == sample_factors["factor1"]
    ), "The factors are not merged correctly"


def test_alpha_diversity_parametrized(sample_tables_dict, sample_factors):
    """Tests the alpha_diversity_parametrized function."""
    result = alpha_diversity_parametrized(
        sample_tables_dict, "sample_table", sample_factors
    )

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected columns
    expected_columns = ["ref_code", "Shannon", "factor1"]
    assert all(
        col in result.columns for col in expected_columns
    ), f"Expected columns {expected_columns}, but got {result.columns.tolist()}"

    # Check if the Shannon index values are calculated correctly
    expected_shannon = sample_tables_dict["sample_table"].apply(
        lambda row: shannon_index(row[1:]), axis=1
    )
    assert all(
        result["Shannon"].round(6) == expected_shannon.round(6)
    ), "The Shannon index values are not calculated correctly"

    # Check if the factors are merged correctly
    assert all(
        result["factor1"] == sample_factors["factor1"]
    ), "The factors are not merged correctly"


def test_beta_diversity_parametrized(sample_data):
    """
    Tests the beta_diversity_parametrized function.
    """
    taxon = "GO:0001"
    metric = "braycurtis"

    # Mock the diversity_input function to return the input DataFrame
    def mock_diversity_input(df, kind, taxon):
        return df.set_index("ref_code")

    # Replace the diversity_input function with the mock
    from momics.diversity import diversity_input

    diversity_input = mock_diversity_input

    result = beta_diversity_parametrized(sample_data, taxon, metric)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected distances
    expected_beta = beta_diversity(metric, sample_data.set_index("ref_code"))
    pd.testing.assert_frame_equal(
        result, expected_beta.to_data_frame(), check_like=True
    )


def test_diversity_input_alpha(sample_data):
    """
    Tests the diversity_input function for alpha diversity.
    """
    result = diversity_input(sample_data, kind="alpha", taxon="ncbi_tax_id")

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ref_codes and 3 taxons
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # Check if the values are correct
    expected_values = {
        "taxon1": [10, 0, 0],
        "taxon2": [0, 20, 0],
        "taxon3": [0, 0, 30],
    }
    for taxon, values in expected_values.items():
        assert all(
            result[taxon] == values
        ), f"Expected values {values} for {taxon}, but got {result[taxon].tolist()}"


def test_diversity_input_beta(sample_data):
    """
    Tests the diversity_input function for beta diversity.
    """
    result = diversity_input(sample_data, kind="beta", taxon="ncbi_tax_id")

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ref_codes and 3 taxons
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # Check if the values are normalized correctly
    expected_values = {
        "taxon1": [1.0, 0.0, 0.0],
        "taxon2": [0.0, 1.0, 0.0],
        "taxon3": [0.0, 0.0, 1.0],
    }
    for taxon, values in expected_values.items():
        assert all(
            result[taxon].round(6) == values
        ), f"Expected values {values} for {taxon}, but got {result[taxon].tolist()}"


def test_get_key_column():
    """
    Tests the get_key_column function.
    """
    # Test known table names
    assert get_key_column("go") == "id", "Expected 'id' for table 'go'"
    assert get_key_column("go_slim") == "id", "Expected 'id' for table 'go_slim'"
    assert get_key_column("ips") == "accession", "Expected 'accession' for table 'ips'"
    assert get_key_column("ko") == "entry", "Expected 'entry' for table 'ko'"
    assert get_key_column("pfam") == "entry", "Expected 'entry' for table 'pfam'"

    # Test unknown table name
    with pytest.raises(ValueError, match="Unknown table: unknown_table"):
        get_key_column("unknown_table")


def test_alpha_input(sample_tables_dict):
    """
    Tests the alpha_input function.
    """
    table_name = "go"
    result = alpha_input(sample_tables_dict, table_name)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ids and 3 ref_codes
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # Check if the values are correct
    expected_values = {
        "sample1": [10, 0, 0],
        "sample2": [0, 20, 0],
        "sample3": [0, 0, 30],
    }
    for ref_code, values in expected_values.items():
        assert all(
            result[ref_code] == values
        ), f"Expected values {values} for {ref_code}, but got {result[ref_code].tolist()}"
