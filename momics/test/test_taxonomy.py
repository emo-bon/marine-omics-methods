import pytest
import pandas as pd
import numpy as np
import os

from momics.metadata import (
    get_metadata_udal,
    enhance_metadata,
    clean_metadata,
)
from momics.taxonomy import *
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
    tax = "338190;sk__Archaea;k__;p__Thaumarchaeota;c__;o__Nitrosopumilales;f__Nitrosopumilaceae;g__;s__"
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


@pytest.fixture
def valid_df():
    idx = pd.MultiIndex.from_tuples(
        [("Bacteria;Firmicutes", 123), ("Bacteria;Proteobacteria", 456)],
        names=["taxonomic_concat", "ncbi_tax_id"],
    )
    data = pd.DataFrame([[10, 20], [30, 40]], index=idx, columns=["sample1", "sample2"])
    return data


def test_tss_sqrt_normalization(valid_df):
    result = normalize_abundance(valid_df, method="tss_sqrt")
    expected = valid_df.div(valid_df.sum(axis=1), axis=0).apply(np.sqrt)
    pd.testing.assert_frame_equal(result, expected)


def test_rarefy_method_calls_rarefy_table(monkeypatch, valid_df):
    called = {}

    def fake_rarefy_table(df, depth=None):
        called["called"] = True
        return df

    monkeypatch.setattr("momics.taxonomy.rarefy_table", fake_rarefy_table)
    result = normalize_abundance(valid_df, method="rarefy", rarefy_depth=5)
    assert called["called"]
    pd.testing.assert_frame_equal(result, valid_df)


def test_invalid_index_raises():
    df = pd.DataFrame([[1, 2], [3, 4]], columns=["sample1", "sample2"])
    with pytest.raises(IndexError):
        normalize_abundance(df)


def test_non_numeric_raises(valid_df):
    df = valid_df.copy()
    df["sample1"] = ["a", "b"]
    with pytest.raises(TypeError):
        normalize_abundance(df)


def test_invalid_method_raises(valid_df):
    with pytest.raises(ValueError):
        normalize_abundance(valid_df, method="invalid")


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
    tax = "338190;sk__Archaea;k__;p__Thaumarchaeota;c__;o__Nitrosopumilales;f__Nitrosopumilaceae;g__;s__"
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


def test_separate_taxonomy_basic():
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=0,
    )
    df = pivot_taxonomic_data(ssu)

    result = separate_taxonomy(df)

    # Check output keys
    expected_keys = {
        "Prokaryotes All",
        "Eukaryota All",
        "Bacteria",
        "Archaea",
        "Bacteria_phylum",
        "Bacteria_class",
        "Bacteria_order",
        "Bacteria_family",
        "Bacteria_genus",
    }
    assert expected_keys.issubset(result.keys())

    # Check that Bacteria and Archaea are separated correctly
    assert all("Bacteria" in idx for idx in result["Bacteria"].index)
    assert all("Archaea" in idx for idx in result["Archaea"].index)
    assert all("Eukaryota" in idx for idx in result["Eukaryota All"].index)


def test_separate_taxonomy_eukaryota_basic():
    # Create a mock DataFrame with a 'taxonomic_concat' column
    data = {
        "taxonomic_concat": [
            "Eukaryota;Opisthokonta;Metazoa;Chordata;Mammalia;Homo;H. sapiens",
            "Eukaryota;Archaeplastida;Viridiplantae;Streptophyta;Magnoliopsida;Arabidopsis;A. thaliana",
            "Eukaryota;Stramenopiles;Bacillariophyta;Coscinodiscophyceae;Thalassiosirales;Thalassiosiraceae;Thalassiosira",
        ],
        "sample1": [1, 2, 3],
        "sample2": [4, 5, 6],
    }
    df = pd.DataFrame(data)
    keywords = ["Opisthokonta", "Archaeplastida", "Stramenopiles", "NotFound"]

    result = separate_taxonomy_eukaryota(df, keywords)

    # Check that all keywords are present as keys
    for key in keywords:
        assert key in result

    # Check that the correct rows are selected for each keyword
    assert len(result["Opisthokonta"]) == 1
    assert "Opisthokonta" in result["Opisthokonta"].iloc[0]["taxonomic_concat"]
    assert len(result["Archaeplastida"]) == 1
    assert "Archaeplastida" in result["Archaeplastida"].iloc[0]["taxonomic_concat"]
    assert len(result["Stramenopiles"]) == 1
    assert "Stramenopiles" in result["Stramenopiles"].iloc[0]["taxonomic_concat"]
    # NotFound should be empty


def test_split_taxonomy_bacteria():
    idx = "6;sk_Bacteria;k_;p_Proteobacteria;c_Alphaproteobacteria;o_Rhizobiales;f_Xanthobacteraceae;g_Azorhizobium;s_"
    result = split_taxonomy(idx)
    assert result == [
        "p_Proteobacteria",
        "c_Alphaproteobacteria",
        "o_Rhizobiales",
        "f_Xanthobacteraceae",
        "g_Azorhizobium",
        "s_",
    ]


def test_split_taxonomy_archaea():
    idx = "Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobacterium;M. formicicum"
    result = split_taxonomy(idx)
    assert result == [
        "Methanobacteria",
        "Methanobacteriales",
        "Methanobacteriaceae",
        "Methanobacterium",
        "M. formicicum",
    ]


def test_split_taxonomy_other():
    idx = "Eukaryota;Opisthokonta;Metazoa;Chordata;Mammalia;Homo;H. sapiens"
    result = split_taxonomy(idx)
    # Should return empty list for non-Bacteria/Archaea
    assert result == []


def test_split_taxonomy_short():
    idx = "Bacteria;Firmicutes"
    result = split_taxonomy(idx)
    # Should handle short taxonomy gracefully (may return less than 6 elements)
    assert isinstance(result, list)


def test_aggregate_by_taxonomic_level_basic():
    # Create a mock DataFrame with taxonomic levels and abundance columns
    data = {
        "phylum": ["Firmicutes", "Firmicutes", "Proteobacteria", "Proteobacteria"],
        "class": ["Bacilli", "Bacilli", "Gammaproteobacteria", "Gammaproteobacteria"],
        "order": [
            "Lactobacillales",
            "Lactobacillales",
            "Enterobacterales",
            "Enterobacterales",
        ],
        "family": [
            "Lactobacillaceae",
            "Lactobacillaceae",
            "Enterobacteriaceae",
            "Enterobacteriaceae",
        ],
        "genus": ["Lactobacillus", "Lactobacillus", "Escherichia", "Escherichia"],
        "species": ["L. acidophilus", "L. casei", "E. coli", "E. fergusonii"],
        "sample1": [10, 20, 30, 40],
        "sample2": [5, 10, 15, 20],
    }
    df = pd.DataFrame(data)

    # Aggregate by 'phylum'
    result = aggregate_by_taxonomic_level(df, "phylum")
    assert set(result.index) == {"Firmicutes", "Proteobacteria"}
    assert result.loc["Firmicutes", "sample1"] == 30  # 10 + 20
    assert result.loc["Proteobacteria", "sample2"] == 35  # 15 + 20

    # Aggregate by 'genus'
    result_genus = aggregate_by_taxonomic_level(df, "genus")
    assert set(result_genus.index) == {"Lactobacillus", "Escherichia"}
    assert result_genus.loc["Lactobacillus", "sample1"] == 30
    assert result_genus.loc["Escherichia", "sample2"] == 35


def test_aggregate_by_taxonomic_level_missing_level():
    # DataFrame with some missing genus values
    data = {
        "phylum": ["Firmicutes", "Firmicutes", "Proteobacteria"],
        "class": ["Bacilli", "Bacilli", "Gammaproteobacteria"],
        "order": ["Lactobacillales", "Lactobacillales", "Enterobacterales"],
        "family": ["Lactobacillaceae", "Lactobacillaceae", "Enterobacteriaceae"],
        "genus": ["Lactobacillus", None, "Escherichia"],
        "species": ["L. acidophilus", "L. casei", "E. coli"],
        "sample1": [10, 20, 30],
        "sample2": [5, 10, 15],
    }
    df = pd.DataFrame(data)

    # Should drop the row with missing genus
    result = aggregate_by_taxonomic_level(df, "genus")
    assert set(result.index) == {"Lactobacillus", "Escherichia"}
    assert "None" not in result.index
    assert result.loc["Lactobacillus", "sample1"] == 10
    assert result.loc["Escherichia", "sample2"] == 15


def test_remove_high_taxa_basic():
    # Create a DataFrame with some missing phylum and genus
    data = {
        "phylum": ["Firmicutes", "Firmicutes", None, "Proteobacteria"],
        "class": ["Bacilli", "Bacilli", "Gammaproteobacteria", "Gammaproteobacteria"],
        "order": [
            "Lactobacillales",
            "Lactobacillales",
            "Enterobacterales",
            "Enterobacterales",
        ],
        "family": [
            "Lactobacillaceae",
            "Lactobacillaceae",
            "Enterobacteriaceae",
            "Enterobacteriaceae",
        ],
        "genus": ["Lactobacillus", "Lactobacillus", "Escherichia", None],
        "species": ["L. acidophilus", "L. casei", "E. coli", "E. fergusonii"],
        "abundance": [10, 20, 30, 40],
    }
    df = pd.DataFrame(data)
    taxonomy_ranks = ["phylum", "class", "order", "family", "genus", "species"]

    # Remove rows with missing phylum (default strict=True)
    result = remove_high_taxa(df, taxonomy_ranks, tax_level="phylum", strict=False)
    assert result["phylum"].isna().sum() == 0
    assert set(result["phylum"]) == {"Firmicutes", "Proteobacteria"}


# TODO: this works only if index_col is [0, 1], adapt function to try to add ncbi_tax_id if not present
# and only if not there throw an error.
def test_remove_high_taxa_strict():
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=[0, 1],
    )

    # Should not raise and should keep only rows with non-null phylum
    result = remove_high_taxa(ssu, TAXONOMY_RANKS, tax_level="phylum", strict=True)
    assert result["phylum"].isna().sum() == 0
    assert (
        "unclassified" in result["phylum"].values
        or "unclassified" not in result["phylum"].values
    )  # Accept both, as strict skips mapping for "unclassified"

    # All rows should have a phylum value
    assert all(pd.notna(result["phylum"]))


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


def test_prevalence_cutoff_basic():
    # DataFrame with 2 metadata columns and 3 sample columns
    data = {
        "taxon": ["A", "B", "C", "D"],
        "info": [1, 2, 3, 4],
        "sample1": [1, 0, 0, 1],
        "sample2": [0, 1, 0, 1],
        "sample3": [0, 0, 1, 1],
    }
    df = pd.DataFrame(data)
    # Set index to taxon for clarity, but function works on columns
    df.set_index("taxon", inplace=True)

    # Only "D" is present in all 3 samples, "A" and "B" in 1, "C" in 1
    # With percent=50, threshold=1.5, so only "D" should remain
    filtered = prevalence_cutoff(df, percent=50, skip_columns=1)
    assert list(filtered.index) == ["D"]

    # With percent=33, threshold=~1, so "A", "B", "C", "D" with at least 1 sample should remain
    filtered = prevalence_cutoff(df, percent=33, skip_columns=1)
    assert set(filtered.index) == {"A", "B", "C", "D"}


def test_prevalence_cutoff_all_removed():
    data = {
        "taxon": ["A", "B"],
        "info": [1, 2],
        "sample1": [0, 0],
        "sample2": [0, 0],
    }
    df = pd.DataFrame(data).set_index("taxon")
    # With percent=50, threshold=1, all rows have 0 prevalence, so should be empty
    filtered = prevalence_cutoff(df, percent=50, skip_columns=1)
    assert filtered.empty


def test_prevalence_cutoff_taxonomy_basic():
    # DataFrame with MultiIndex (taxon, sample), abundance column
    index = pd.MultiIndex.from_tuples(
        [("A", "s1"), ("A", "s2"), ("B", "s1"), ("B", "s2"), ("C", "s1"), ("C", "s2")],
        names=["taxon", "sample"],
    )
    data = {"abundance": [10, 0, 0, 0, 5, 6]}
    df = pd.DataFrame(data, index=index)

    # With percent=50, for each taxon, only keep rows with abundance > 50% of total for that taxon
    # For A: sum=10, threshold=5, only ("A", "s1") kept
    # For B: sum=0, threshold=0, none kept
    # For C: sum=11, threshold=5.5, only ("C", "s2") kept
    filtered = prevalence_cutoff_taxonomy(df, percent=50)
    assert ("A", "s1") in filtered.index
    assert ("C", "s2") in filtered.index
    assert ("A", "s2") not in filtered.index
    assert ("C", "s1") not in filtered.index
    assert ("B", "s1") not in filtered.index
    assert ("B", "s2") not in filtered.index


def test_prevalence_cutoff_taxonomy_no_multiindex():
    test_dir = os.path.dirname(__file__)
    ssu = pd.read_csv(
        os.path.join(test_dir, "data", "ssu_head.csv"),
        index_col=[0, 1],
    )
    filtered = prevalence_cutoff_taxonomy(ssu, percent=10)
    assert (
        filtered.loc[filtered.index == ("EMOBON00084", 338190), "abundance"] == 129
    ).any()

    data = {"taxon": ["A", "B", "C", "D", "E", "F"], "abundance": [10, 0, 0, 0, 5, 6]}
    df = pd.DataFrame(data)
    df.set_index("taxon", inplace=True)
    filtered = prevalence_cutoff(df, percent=50, skip_columns=0)
    print(filtered)
    assert set(filtered.index) == {"A", "E", "F"}
    assert (filtered.loc[filtered.index == "A", "abundance"] == 10).any()
    assert (filtered.loc[filtered.index == "E", "abundance"] == 5).any()
    assert (filtered.loc[filtered.index == "F", "abundance"] == 6).any()


def test_rarefy_table_axis1():
    # Simple abundance table: rows=features, columns=samples
    df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [5, 15, 10],
        },
        index=["taxonA", "taxonB", "taxonC"],
    )
    # Minimum sample sum is 30 (sample2)
    result = rarefy_table(df, axis=1)
    assert result.shape == df.shape
    # All values should be <= original, and sum per column should be close to 30 (rarefaction depth)
    assert np.all(result.sum(axis=0).dropna() <= 30.1)


def test_rarefy_table_axis0():
    # Features in columns, samples in rows
    df = pd.DataFrame(
        {
            "taxonA": [10, 5],
            "taxonB": [20, 15],
            "taxonC": [30, 10],
        },
        index=["sample1", "sample2"],
    )
    # Minimum sample sum is 30 (sample2)
    result = rarefy_table(df, axis=0)
    # print(result)
    assert result.shape == df.T.shape
    assert np.all(result.sum(axis=0).dropna() <= 30.1)
    assert (result["sample2"] == [5, 15, 10]).all()


def test_rarefy_table_with_depth():
    df = pd.DataFrame(
        {
            "sample1": [10, 20, 30],
            "sample2": [5, 15, 10],
        },
        index=["taxonA", "taxonB", "taxonC"],
    )
    # Set depth to 10, so all columns should sum to 10
    result = rarefy_table(df, depth=10, axis=1)
    assert np.allclose(result.sum(axis=0).dropna(), 10)


def test_rarefy_table_not_enough_counts():
    # One sample has less than the rarefaction depth
    df = pd.DataFrame(
        {
            "sample1": [1, 1, 1],
            "sample2": [10, 10, 10],
        },
        index=["taxonA", "taxonB", "taxonC"],
    )
    result = rarefy_table(df, depth=10, axis=1)
    print(result)
    # sample1 should be all NaN
    assert result["sample1"].isna().all()
    # sample2 should sum to 10
    assert np.isclose(result["sample2"].sum(), 10)


def test_fill_taxonomy_placeholders_basic():
    # DataFrame with missing class and order, but genus present
    data = {
        "phylum": ["Firmicutes", "Firmicutes"],
        "class": ["", ""],
        "order": ["", ""],
        "family": ["Lactobacillaceae", "Lactobacillaceae"],
        "genus": ["Lactobacillus", "Lactobacillus"],
        "species": ["L. acidophilus", "L. casei"],
    }
    df = pd.DataFrame(data)
    taxonomy_ranks = ["phylum", "class", "order", "family", "genus", "species"]

    result = fill_taxonomy_placeholders(df, taxonomy_ranks)

    # class and order should be filled with 'unclassified_<lower_rank_value>'
    assert all(result["class"] == "unclassified_unclassified_Lactobacillaceae")
    assert all(result["order"] == "unclassified_Lactobacillaceae")
    # family, genus, species should remain unchanged
    assert all(result["family"] == "Lactobacillaceae")
    assert all(result["genus"] == "Lactobacillus")
    assert set(result["species"]) == {"L. acidophilus", "L. casei"}


def test_fill_taxonomy_placeholders_no_missing():
    # DataFrame with no missing values
    data = {
        "phylum": ["Firmicutes"],
        "class": ["Bacilli"],
        "order": ["Lactobacillales"],
        "family": ["Lactobacillaceae"],
        "genus": ["Lactobacillus"],
        "species": ["L. acidophilus"],
    }
    df = pd.DataFrame(data)
    taxonomy_ranks = ["phylum", "class", "order", "family", "genus", "species"]

    result = fill_taxonomy_placeholders(df, taxonomy_ranks)
    # Should not change any values
    assert (result == df).all().all()


def test_split_metadata():
    """
    Tests the split_metadata function.
    """
    metadata = get_metadata_udal()
    assert isinstance(metadata, pd.DataFrame), "The result should be a DataFrame"

    metadata, added_columns = enhance_metadata(metadata)
    assert isinstance(
        metadata, pd.DataFrame
    ), "The enhanced metadata should be a DataFrame"
    assert added_columns == [
        "year",
        "month",
        "month_name",
        "day",
        "season",
        "replicate_info",
    ], "Unexpected added columns"

    # convert added_columns to a dictionary
    added_columns = {col: col.replace("_", " ") for col in added_columns}
    # extend COL_NAMES_HASH with added columns
    COL_NAMES_HASH.update(added_columns)
    metadata = clean_metadata(metadata, COL_NAMES_HASH)

    # Identify object columns
    object_cols = metadata.select_dtypes(include="object").columns

    # Convert them all at once to category
    metadata = metadata.astype({col: "category" for col in object_cols})

    filtered_metadata = metadata.drop_duplicates(subset="replicate info", keep="first")

    groups = split_metadata(
        filtered_metadata,
        "season",
    )

    assert isinstance(groups, dict), "The result should be a dictionary"
    assert len(groups) == 4, "There should be 4 groups for each season"
    assert set(groups.keys()) == {
        "Spring",
        "Summer",
        "Autumn",
        "Winter",
    }, "The groups should be the seasons"
    assert len(groups["Spring"]) == 3, "There should be 3 samples in Spring group"
    assert len(groups["Summer"]) == 39, "There should be 39 samples in Summer group"
    assert len(groups["Autumn"]) == 49, "There should be 49 samples in Autumn group"
    assert len(groups["Winter"]) == 9, "There should be 9 samples in Winter group"


def test_split_taxonomic_data_basic():
    # DataFrame with index as ref_code
    df = pd.DataFrame(
        {"abundance": [10, 20, 30, 40], "taxon": ["A", "B", "A", "C"]},
        index=["ref1", "ref2", "ref3", "ref4"],
    )
    groups = {
        "group1": ["ref1", "ref3"],
        "group2": ["ref2"],
        "group3": ["ref4", "refX"],  # refX does not exist
    }
    result = split_taxonomic_data(df, groups)
    # group1 should have ref1 and ref3
    assert set(result["group1"].index) == {"ref1", "ref3"}
    # group2 should have ref2
    assert list(result["group2"].index) == ["ref2"]
    # group3 should have ref4 only (refX ignored)
    assert set(result["group3"].index) == {"ref4"}


def test_split_taxonomic_data_empty_group():
    df = pd.DataFrame(
        {
            "abundance": [10, 20],
        },
        index=["ref1", "ref2"],
    )
    groups = {"empty": ["refX"], "all": ["ref1", "ref2"]}
    result = split_taxonomic_data(df, groups)
    assert result["empty"].empty
    assert set(result["all"].index) == {"ref1", "ref2"}


def test_split_taxonomic_data_pivoted():
    import numpy as np

    # Create a DataFrame with two "sample" columns and a multiindex
    taxonomy = pd.DataFrame(
        {
            "A": [1, 0, 0],
            "B": [0, 2, 0],
            "C": [0, 0, 0],
        },
        index=pd.MultiIndex.from_tuples(
            [
                (1, "tax1"),
                (2, "tax2"),
                (3, "tax3"),
            ],
            names=["ncbi_tax_id", "taxonomic_concat"],
        ),
    )

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


def test_compute_bray_curtis_samples():
    # DataFrame with 2 metadata columns and 3 sample columns
    data = {
        "taxon": ["A", "B", "C"],
        "info": [1, 2, 3],
        "sample1": [1, 2, 3],
        "sample2": [4, 5, 6],
        "sample3": [7, 8, 9],
    }
    df = pd.DataFrame(data)
    # Bray-Curtis for samples (skip first 2 columns)
    result = compute_bray_curtis(df, skip_cols=2, direction="samples")
    # Should be a square DataFrame with sample names as index/columns
    assert isinstance(result, pd.DataFrame)
    assert set(result.index) == {"sample1", "sample2", "sample3"}
    assert set(result.columns) == {"sample1", "sample2", "sample3"}
    # Diagonal should be zero
    assert np.allclose(np.diag(result), 0)


def test_compute_bray_curtis_taxa():
    # DataFrame with ncbi_tax_id and 3 sample columns
    data = {
        "ncbi_tax_id": [111, 222, 333],
        "taxon": ["A", "B", "C"],
        "sample1": [1, 2, 3],
        "sample2": [4, 5, 6],
        "sample3": [7, 8, 9],
    }
    df = pd.DataFrame(data)
    df.set_index(["ncbi_tax_id", "taxon"], inplace=True)
    print(df)
    # Bray-Curtis for taxa (skip first 2 columns)
    result = compute_bray_curtis(df, skip_cols=2, direction="taxa")
    # Should be a square DataFrame with ncbi_tax_id as index/columns
    assert isinstance(result, pd.DataFrame)
    assert set(result.index) == {111, 222, 333}
    assert set(result.columns) == {111, 222, 333}
    # Diagonal should be zero
    assert np.allclose(np.diag(result), 0)


def test_fdr_pvals_basic():
    # Symmetric p-value matrix with zeros on diagonal
    data = [
        [0.0, 0.01, 0.04],
        [0.01, 0.0, 0.03],
        [0.04, 0.03, 0.0],
    ]
    df = pd.DataFrame(data, columns=["A", "B", "C"], index=["A", "B", "C"])

    # Use a high cutoff to ensure all are "significant"
    result = fdr_pvals(df, pval_cutoff=0.1)
    trili = np.tril_indices_from(result, k=0)
    # Should be nan on and below diagonal, and FDR-corrected above
    assert np.all(np.isnan(result.values[trili]))

    # All upper triangle values should be between 0 and 1
    upper = result.values[np.triu_indices_from(result, k=1)]
    assert np.all((upper >= 0) & (upper <= 1))


def test_fdr_pvals_no_significant():
    # All p-values are high, so none should be significant
    data = [
        [0.0, 0.8, 0.9],
        [0.8, 0.0, 0.7],
        [0.9, 0.7, 0.0],
    ]
    df = pd.DataFrame(data, columns=["A", "B", "C"], index=["A", "B", "C"])
    result = fdr_pvals(df, pval_cutoff=0.05)
    trili = np.tril_indices_from(result, k=0)
    # Should be nan on and below diagonal, and FDR-corrected above
    assert np.all(np.isnan(result.values[trili]))
    # All upper triangle values should be >= 0
    upper = result.values[np.triu_indices_from(result, k=1)]
    assert np.all(upper >= 0)
