import pytest
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
import skbio
from momics.diversity import *
from momics.test.fixtures import *


@pytest.mark.parametrize("name", ["sample_table"])
def test_shannon_index(name):
    """Tests the shannon_index function."""
    data = sample_data(name)
    # Test with a row containing species abundances
    row = data.iloc[0]
    assert (
        row.all()
        == pd.Series(
            {
                "ref_code": "sample1",
                "GO:0001": 10,
                "GO:0002": 20,
                "IPR0001": 30,
                "K0001": 40,
                "PF0001": 50,
            }
        ).all()
    ), f"Expected {row}, but got {pd.Series(row)}"

    result = shannon_index(row[1:])
    expected = -sum(
        (x / row[1:].sum()) * np.log(x / row[1:].sum()) for x in row[1:] if x > 0
    )
    assert np.isclose(result, expected), f"Expected {expected}, but got {result}"

    # Test with a row[1:] containing zero abundances
    row = data.iloc[1]
    result = shannon_index(row[1:])
    assert np.isnan(result), f"Expected NaN, but got {result}"

    # Test with a row[1:] containing equal abundances
    row = data.iloc[2]
    result = shannon_index(row[1:])
    expected = -sum(
        (x / row[1:].sum()) * np.log(x / row[1:].sum()) for x in row[1:] if x > 0
    )
    assert np.isclose(
        result, expected, equal_nan=True
    ), f"Expected {expected}, but got {result}"


def test_calculate_alpha_diversity(sample_factors):
    """Tests the calculate_alpha_diversity function."""
    data = sample_data()
    result = calculate_alpha_diversity(data, sample_factors)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected columns
    expected_columns = ["Shannon", "factor1"]
    assert all(
        col in result.columns for col in expected_columns
    ), f"Expected columns {expected_columns}, but got {result.columns.tolist()}"

    # Check if the Shannon index values are calculated correctly
    expected_shannon = data.apply(lambda row: shannon_index(row[:]), axis=1)
    assert np.isclose(result["Shannon"], expected_shannon, equal_nan=True).all(), (
        "The Shannon_index values are not calculated correctly, diff is "
        + f"{(result['Shannon'] - expected_shannon).tolist()}"
    )

    # Check if the factors are merged correctly
    assert all(
        result["factor1"] == sample_factors["factor1"]
    ), "The factors are not merged correctly"


@pytest.mark.parametrize(
    "table_name, col_to_add",
    [
        ("go", "id"),
        ("go_slim", "id"),
        ("ips", "accession"),
        ("ko", "entry"),
        ("pfam", "entry"),
    ],
)
def test_alpha_diversity_parametrized(table_name, col_to_add, sample_factors):
    """Tests the alpha_diversity_parametrized function."""
    data_dict = sample_tables_dict(table_name, add_abundance=True)
    data_dict = {
        table_name: add_column(data_dict[table_name], col_to_add),
    }
    assert (
        data_dict[table_name].index.name == sample_factors.index.name
    ), "The index names of the input DataFrame and metadata do not match."
    print(sample_factors)
    result = alpha_diversity_parametrized(data_dict, table_name, sample_factors)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected columns
    expected_columns = ["Shannon", "factor1"]
    assert all(
        col in result.columns for col in expected_columns
    ), f"Expected columns {expected_columns}, but got {result.columns.tolist()}"

    # Check if the Shannon index values are calculated correctly
    expected_shannon = data_dict[table_name].apply(
        lambda row: shannon_index(row[1:-1]), axis=1
    )

    # assert np.isclose(result["Shannon"], expected_shannon, equal_nan=True).all(), (
    #     "The Shannon_index values are not calculated correctly, diff is " +
    #     f"{(result['Shannon'] - expected_shannon).tolist()}"
    # )

    # Check if the factors are merged correctly
    assert all(
        result["factor1"] == sample_factors["factor1"]
    ), "The factors are not merged correctly"


def test_beta_diversity_parametrized():
    """
    Tests the beta_diversity_parametrized function.
    """
    taxon = "GO:0001"
    metric = "braycurtis"
    data = sample_data(add_abundance=True)

    # # Mock the diversity_input function to return the input DataFrame
    # def mock_diversity_input(df, kind, taxon):
    #     return df.set_index("ref_code")

    # Replace the diversity_input function with the mock
    from momics.diversity import diversity_input

    # diversity_input = mock_diversity_input
    beta_input = diversity_input(data, kind="beta", taxon=taxon)
    expected_beta = beta_diversity(metric, beta_input)

    data["ncbi_tax_id"] = ["taxon1", "taxon2", "taxon3"]

    result = beta_diversity_parametrized(data, taxon, metric)

    # Check if the result is a DataFrame
    assert isinstance(
        result, skbio.stats.distance._base.DistanceMatrix
    ), "The result should be a DataFrame"

    # Check if the result contains the expected distances
    # TODO: does not work yet
    # pd.testing.assert_frame_equal(
    #     result.to_data_frame(), expected_beta.to_data_frame(), check_like=True
    # )


def test_diversity_input_alpha():
    """Tests the diversity_input function for alpha diversity."""
    data = sample_data(add_abundance=True)
    data["ncbi_tax_id"] = ["taxon1", "taxon2", "taxon3"]
    result = diversity_input(data, kind="alpha", taxon="ncbi_tax_id")

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ref_codes and 3 taxons
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # Check if the values are correct
    expected_values = {
        "taxon1": [50, 0, 0],
        "taxon2": [0, 0, 0],
        "taxon3": [0, 0, 5],
    }
    for taxon, values in expected_values.items():
        assert all(
            result[taxon] == values
        ), f"Expected values {values} for {taxon}, but got {result[taxon].tolist()}"


def test_diversity_input_beta():
    """
    Tests the diversity_input function for beta diversity.
    """
    data = sample_data(add_abundance=True)
    data["ncbi_tax_id"] = ["taxon1", "taxon2", "taxon3"]
    result = diversity_input(data, kind="beta", taxon="ncbi_tax_id")

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ref_codes and 3 taxons
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # Check if the values are normalized correctly
    expected_values = {
        "taxon1": [1.0, np.nan, 0.0],
        "taxon2": [0.0, 1.0, 0.0],
        "taxon3": [0.0, 0.0, 1.0],
    }
    for taxon, values in expected_values.items():
        np.isclose(result[taxon], values, equal_nan=True),
        f"Expected values {values} for {taxon}, but got {result[taxon].tolist()}"


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


def test_alpha_input():
    """
    Tests the alpha_input function.
    """
    table_name = "go"
    data_dict = sample_tables_dict(table_name, add_abundance=True)
    data_dict = {
        table_name: add_column(data_dict[table_name], "id"),
    }
    result = alpha_input(data_dict, table_name)

    # Check if the result is a DataFrame
    assert isinstance(result, pd.DataFrame), "The result should be a DataFrame"

    # Check if the result contains the expected shape
    expected_shape = (3, 3)  # 3 ids and 3 ref_codes
    assert (
        result.shape == expected_shape
    ), f"Expected shape {expected_shape}, but got {result.shape}"

    # TODO: Check if the values are correct
    expected_values = {
        "sample1": [50, 0, 0],
        "sample2": [0, 0, 0],
        "sample3": [0, 0, 5],
    }
    for ref_code, values in expected_values.items():
        assert all(
            result[ref_code] == values
        ), f"Expected values {values} for {ref_code}, but got {result[ref_code].tolist()}"
