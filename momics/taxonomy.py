import pandas as pd
import numpy as np
from typing import List, Dict


"""
inspired by Andrqej Tkacz, must be cleaned up and documented
"""


def pivot_taxonomic_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepares the taxonomic data (LSU and SSU tables) for analysis. Apart from
    pivoting, it also normalizes and calculates square root of the abundances,
    ie Total Sum Scaling (TSS) followed by Square Root Transformation.

    TODO: refactor scaling to a new method and offer different options.

    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.

    Returns:
        pd.DataFrame: A pivot table with taxonomic data.
    """
    # Select relevant columns
    df["taxonomic_concat"] = (
        df["ncbi_tax_id"].astype(str)
        + ";sk_"
        + df["superkingdom"].fillna("")
        + ";k_"
        + df["kingdom"].fillna("")
        + ";p_"
        + df["phylum"].fillna("")
        + ";c_"
        + df["class"].fillna("")
        + ";o_"
        + df["order"].fillna("")
        + ";f_"
        + df["family"].fillna("")
        + ";g_"
        + df["genus"].fillna("")
        + ";s_"
        + df["species"].fillna("")
    )
    pivot_table = df.pivot_table(
        index=["ncbi_tax_id", "taxonomic_concat"],
        columns="ref_code",
        values="abundance",
    ).fillna(0)
    pivot_table = pivot_table.reset_index()
    # change inex name
    pivot_table.columns.name = None

    # normalize values
    pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: x / x.sum())
    pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: np.sqrt(x))
    return pivot_table


def separate_taxonomy(
    df: pd.DataFrame, eukaryota_keywords: List[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    Separate the taxonomic data into different categories based on the index names.
    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information (LSU/SSU tables).
        eukaryota_keywords (List[str]): List of keywords to filter Eukaryota data.
    Returns:
        Dict[str, pd.DataFrame]: A dictionary containing separate DataFrames for Prokaryotes and Eukaryota.
    """
    # Separate rows based on "Bacteria", "Archaea", and "Eukaryota" entries
    prokaryotes_all = df[df.index.str.contains("Bacteria|Archaea", regex=True)]
    eukaryota_all = df[df.index.str.contains("Eukaryota", regex=True)]

    # Further divide "Prokaryotes all" into "Bacteria" and "Archaea"
    bacteria = prokaryotes_all[prokaryotes_all.index.str.contains("Bacteria")]
    archaea = prokaryotes_all[prokaryotes_all.index.str.contains("Archaea")]

    # Apply taxonomy splitting to the index
    taxonomy_levels = bacteria.index.to_series().apply(split_taxonomy)
    taxonomy_df = pd.DataFrame(
        taxonomy_levels.tolist(),
        columns=["phylum", "class", "order", "family", "genus", "species"],
        index=bacteria.index,
    )

    # Combine taxonomy with the abundance data
    bacteria_data = pd.concat([taxonomy_df, bacteria], axis=1)

    # Aggregate at each taxonomic level and save to CSV
    taxonomic_levels = ["phylum", "class", "order", "family", "genus"]
    bacteria_levels_dict = {}
    for level in taxonomic_levels:
        aggregated_df = aggregate_by_taxonomic_level(bacteria_data, level)
        # Standardize the values so each column sums to 100
        aggregated_df_normalized = (
            aggregated_df.div(aggregated_df.sum(axis=0), axis=1) * 100
        )
        bacteria_levels_dict[f"Bacteria_{level}"] = aggregated_df_normalized

    all_data = {
        "Prokaryotes All": prokaryotes_all,
        "Eukaryota All": eukaryota_all,
        "Bacteria": bacteria,
        "Archaea": archaea,
    }
    # all_data.update(eukaryota_dict)
    all_data.update(bacteria_levels_dict)

    # If eukaryota keywords are provided, separate Eukaryota data
    if eukaryota_keywords:
        eukaryota_dict = separate_taxonomy_eukaryota(eukaryota_all, eukaryota_keywords)
        all_data.update(eukaryota_dict)

    return all_data


def separate_taxonomy_eukaryota(df: pd.DataFrame, eukaryota_keywords: List[str]):
    """
    Separate Eukaryota data into different files based on specific keywords.
    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information (LSU/SSU tables).
        eukaryota_keywords (List[str]): List of keywords to filter Eukaryota data.

    Example keywords:
        eukaryota_keywords = ['Discoba', 'Stramenopiles', 'Rhizaria', 'Alveolata',
                              'Amorphea', 'Archaeoplastida', 'Excavata']
    """
    # Further divide "Eukaryota all" by specific keywords
    eukaryota_dict = {}
    for keyword in eukaryota_keywords:
        subset = df[df["taxonomic_concat"].str.contains(keyword)]
        eukaryota_dict[keyword] = subset

    return eukaryota_dict


def split_taxonomy(index_name: str) -> List[str]:
    """
    Splits the taxonomic string into its components and removes prefixes.
    Args:
        index_name (str): The taxonomic string to split.
    Returns:
        List[str]: A list of taxonomic levels.
    """
    # Remove anything before "Bacteria" or "Archaea"
    if "Bacteria" in index_name:
        taxonomy = index_name.split("Bacteria;", 1)[1].split(";")
    elif "Archaea" in index_name:
        taxonomy = index_name.split("Archaea;", 1)[1].split(";")
    else:
        taxonomy = []
    # Return a list with taxonomic levels up to species
    return taxonomy[1:7]  # ['phylum', 'class', 'order', 'family', 'genus', 'species']


def aggregate_by_taxonomic_level(df: pd.DataFrame, level: str) -> pd.DataFrame:
    """
    Aggregates the DataFrame by a specific taxonomic level and sums abundances across samples.
    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.
        level (str): The taxonomic level to aggregate by (e.g., 'phylum', 'class', etc.).
    Returns:
        pd.DataFrame: A DataFrame aggregated by the specified taxonomic level.
    """
    # Drop rows where the level is missing
    df_level = df.dropna(subset=[level])
    # Group by the specified level and sum abundances across samples (columns)
    df_grouped = df_level.groupby(level).sum(numeric_only=True)
    return df_grouped
