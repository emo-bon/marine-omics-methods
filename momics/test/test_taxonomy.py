import pytest
import pandas as pd
import os

from momics.metadata import (
    get_metadata_udal,
    enhance_metadata,
    clean_metadata,
)
from momics.taxonomy import (
    pivot_taxonomic_data,
    normalize_abundance,
    split_metadata,
    taxon_in_table,
    map_taxa_up,
    split_taxonomic_data_pivoted,
)
from momics.constants import COL_NAMES_HASH_EMO_BON_VRE as COL_NAMES_HASH
from momics.constants import TAXONOMY_RANKS


def test_pivot_taxonomic_data():
    """
    Tests the pivot_taxonomic_data function for the case of simple ref_code index
    """
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=0,
    )
    # this pivot table has multiindex implemented
    pivot = pivot_taxonomic_data(ssu)
    assert isinstance(pivot, pd.DataFrame), "The result should be a DataFrame"
    assert not pivot.empty, "The pivoted DataFrame should not be empty"
    assert pivot.columns.name == "ref_code", "Columns should be named 'ref_code'"
    assert pivot.index.names == [
        "ncbi_tax_id",
        "taxonomic_concat",
    ], "Index names should be 'taxonomic_concat' and 'ncbi_tax_id'"

    expected_columns = ["EMOBON00084"]

    assert all(
        col in pivot.columns for col in expected_columns
    ), f"The pivoted DataFrame should contain the expected columns: {expected_columns}"

    assert pivot.shape[1] == len(
        expected_columns
    ), f"The pivoted DataFrame should have {len(expected_columns)} columns, but got {pivot.shape[1]}"

    assert (
        pivot.shape[0] == 10
    ), "The pivoted DataFrame should have one row for the test data"

    # assert one specific row
    tax = "338190;sk_Archaea;k_;p_Thaumarchaeota;c_;o_Nitrosopumilales;f_Nitrosopumilaceae;g_;s_"
    index_vals = pivot.index.get_level_values("taxonomic_concat")

    assert (
        tax in index_vals
    ), f"The taxonomic_concat column should contain the taxonomic string '{tax}'"
    ncbi_tax_id = pivot.index[pivot.index.get_level_values("taxonomic_concat") == tax][
        0
    ][0]

    assert (
        ncbi_tax_id == 338190
    ), f"The ncbi_tax_id for taxonomic_concat '{tax}' should be 338190"

    # Check if the EMOBON00084 column contains integer values
    assert (
        pivot["EMOBON00084"].dtype == int
    ), "The EMOBON00084 column should contain integer values"


def test_pivot_taxonomic_data_multiindex():
    """
    Tests the pivot_taxonomic_data function for the case of multiindex
    """
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=[0, 1],
    )
    # this pivot table has multiindex implemented
    pivot = pivot_taxonomic_data(ssu)
    assert isinstance(pivot, pd.DataFrame), "The result should be a DataFrame"
    assert not pivot.empty, "The pivoted DataFrame should not be empty"
    assert pivot.columns.name == "ref_code", "Columns should be named 'ref_code'"
    assert pivot.index.names == [
        "ncbi_tax_id",
        "taxonomic_concat",
    ], "Index names should be 'taxonomic_concat' and 'ncbi_tax_id'"

    expected_columns = ["EMOBON00084"]

    assert all(
        col in pivot.columns for col in expected_columns
    ), f"The pivoted DataFrame should contain the expected columns: {expected_columns}"

    assert pivot.shape[1] == len(
        expected_columns
    ), f"The pivoted DataFrame should have {len(expected_columns)} columns, but got {pivot.shape[1]}"

    assert (
        pivot.shape[0] == 10
    ), "The pivoted DataFrame should have one row for the test data"

    # assert one specific row
    tax = "338190;sk_Archaea;k_;p_Thaumarchaeota;c_;o_Nitrosopumilales;f_Nitrosopumilaceae;g_;s_"
    index_vals = pivot.index.get_level_values("taxonomic_concat")
    assert (
        tax in index_vals
    ), f"The taxonomic_concat column should contain the taxonomic string '{tax}'"
    ncbi_tax_id = pivot.index[pivot.index.get_level_values("taxonomic_concat") == tax][
        0
    ][0]

    assert (
        ncbi_tax_id == 338190
    ), f"The ncbi_tax_id for taxonomic_concat '{tax}' should be 338190"

    # Check if the EMOBON00084 column contains integer values
    assert (
        pivot["EMOBON00084"].dtype == int
    ), "The EMOBON00084 column should contain integer values"


@pytest.mark.parametrize("reset_index", [False, True])
def test_normalize_abundance(reset_index):
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=0,
    )
    pivot = pivot_taxonomic_data(ssu)
    normalized = normalize_abundance(pivot)

    assert isinstance(normalized, pd.DataFrame), "The result should be a DataFrame"
    assert not normalized.empty, "The normalized DataFrame should not be empty"
    assert (
        normalized.sum(axis=1) - 1
    ).abs().max() < 1e-6, "Each row in the normalized DataFrame should sum to 1"

    # Check if the original data is preserved
    if not reset_index:
        assert all(
            pivot.index == normalized.index
        ), "The index of the normalized DataFrame should match the original pivoted DataFrame"

    if reset_index:
        pivot = pivot.reset_index()
        with pytest.raises(
            IndexError,
            match="DataFrame must have a multiindex with 'taxonomic_concat' and 'ncbi_tax_id'.",
        ):
            normalize_abundance(pivot)


# TODO: parametrize with different normalization methods
@pytest.mark.parametrize("method", ["tss_sqrt", "rarefy"])
def test_normalize_abundance_methods(method):
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=0,
    )
    pivot = pivot_taxonomic_data(ssu)
    normalized = normalize_abundance(pivot, method=method)

    assert isinstance(normalized, pd.DataFrame), "The result should be a DataFrame"
    assert not normalized.empty, "The normalized DataFrame should not be empty"
    if method != "rarefy":
        assert (
            normalized.sum(axis=1) - 1
        ).abs().max() < 1e-6, "Each row in the normalized DataFrame should sum to 1"


def test_taxon_in_table():
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
    )

    ssu.set_index(["ref_code", "ncbi_tax_id"], inplace=True)
    assert taxon_in_table(ssu, TAXONOMY_RANKS, "Euryarchaeota", "phylum") == 28890
    assert (
        taxon_in_table(ssu, TAXONOMY_RANKS, "Candidatus_Woesearchaeota", "phylum")
        == 1801616
    )
    assert taxon_in_table(ssu, TAXONOMY_RANKS, "Cenarchaeum", "genus") == 46769
    assert taxon_in_table(ssu, TAXONOMY_RANKS, "Archaea", "superkingdom") == 2157
    assert (
        taxon_in_table(ssu, TAXONOMY_RANKS, "Candidatus_Woesearchaeota", "superkingdom")
        == -1
    )

    assert taxon_in_table(ssu, TAXONOMY_RANKS, "Euryarchaeota", "species") == None


def test_map_taxa_up():
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=[0, 1],
    )

    # Map all rows with phylum 'Euryarchaeota' up to the tax_id 28890
    result = map_taxa_up(
        ssu.copy(), taxon="Euryarchaeota", tax_level="phylum", tax_id=28890
    )

    # After mapping, only one row for Euryarchaeota (tax_id=28890) should remain, with summed abundance
    filtered = result[result["phylum"] == "Euryarchaeota"]
    assert (
        filtered.shape[0] == 1
    ), "Should only have one row for Euryarchaeota after mapping"
    assert filtered["abundance"].iloc[0] == 2, "Abundance should be summed (1 + 1)"

    result = map_taxa_up(
        ssu.copy(), taxon="Nanoarchaeota", tax_level="phylum", tax_id=192989
    )
    # The Nanoarchaeota row should remain unchanged
    nano = result[result["phylum"] == "Nanoarchaeota"]
    assert nano.shape[0] == 1
    assert nano["abundance"].iloc[0] == 3, "Abundance for Nanoarchaeota should remain 3"

    result = map_taxa_up(
        ssu.copy(), taxon="Thaumarchaeota", tax_level="phylum", tax_id=651137
    )
    thaum = result[result["phylum"] == "Thaumarchaeota"]
    assert thaum.shape[0] == 1
    assert (
        thaum["abundance"].iloc[0] == 139
    ), "Abundance for Thaumarchaeota should remain 3"


def test_split_metadata():
    """
    Tests the split_metadata function.
    """
    metadata = get_metadata_udal()
    assert isinstance(metadata, pd.DataFrame), "The result should be a DataFrame"

    metadata, added_columns = enhance_metadata(metadata)
    assert isinstance(metadata, pd.DataFrame), "The enhanced metadata should be a DataFrame"
    assert (added_columns == ['year', 'month', 'month_name', 'day', 'season', 'replicate_info'], 
            "Unexpected added columns")
    
    # convert added_columns to a dictionary
    added_columns = {col: col.replace("_", " ") for col in added_columns}
    # extend COL_NAMES_HASH with added columns
    COL_NAMES_HASH.update(added_columns)
    metadata = clean_metadata(metadata, COL_NAMES_HASH)

    # Identify object columns
    object_cols = metadata.select_dtypes(include='object').columns

    # Convert them all at once to category
    metadata = metadata.astype({col: 'category' for col in object_cols})

    filtered_metadata = metadata.drop_duplicates(subset='replicate info', keep='first')
    
    groups = split_metadata(
        filtered_metadata,
        'season',
    )

    assert isinstance(groups, dict), "The result should be a dictionary"
    assert len(groups) == 4, "There should be 4 groups for each season"
    assert set(groups.keys()) == {'Spring', 'Summer', 'Autumn', 'Winter'}, "The groups should be the seasons"
    assert len(groups['Spring']) == 3, "There should be 3 samples in Spring group"
    assert len(groups['Summer']) == 39, "There should be 39 samples in Summer group"
    assert len(groups['Autumn']) == 49, "There should be 49 samples in Autumn group"
    assert len(groups['Winter']) == 9, "There should be 9 samples in Winter group"


def test_split_taxonomic_data_pivoted():
    import numpy as np

    # Create a DataFrame with two "sample" columns and a multiindex
    taxonomy = pd.DataFrame({
        "A": [1, 0, 0],
        "B": [0, 2, 0],
        "C": [0, 0, 0],
    }, index=pd.MultiIndex.from_tuples([
        (1, "tax1"),
        (2, "tax2"),
        (3, "tax3"),
    ], names=["ncbi_tax_id", "taxonomic_concat"]))

    groups = {
        "group1": ["A", "B"],
        "group2": ["C"],
    }

    result = split_taxonomic_data_pivoted(taxonomy, groups)

    # group1 should have two rows (tax1 and tax2), group2 should be empty (all zeros)
    assert "group1" in result
    assert result["group1"].shape == (2, 2)
    assert set(result["group1"].columns) == {"A", "B"}
    assert 1 in result["group1"].index.get_level_values(0)
    assert 2 in result["group1"].index.get_level_values(0)

    # group2 should not be in result (all zeros, filtered out)
    assert "group2" not in result

    # Should raise ValueError if groups is not a dict
    with pytest.raises(ValueError):
        split_taxonomic_data_pivoted(taxonomy, "not_a_dict")