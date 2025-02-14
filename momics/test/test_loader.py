import pytest
import pandas as pd

# import fastparquet

from momics.loader.parquets import load_parquets


@pytest.fixture
def folder(tmp_path):
    # Create a temporary folder
    folder = tmp_path / "data"
    folder.mkdir()
    # Create a sample DataFrame and save it as a parquet file
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
    df.to_parquet(folder / "test.parquet")
    return folder


def test_load_parquets(folder):
    # Load the parquet file
    data = load_parquets(folder)
    # Check that the data is loaded correctly
    assert "test" in data
    assert isinstance(data["test"], pd.DataFrame)
    assert isinstance(data, dict)
    assert data["test"].equals(pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}))
