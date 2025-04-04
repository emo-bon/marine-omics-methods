import pandas as pd
import numpy as np


"""
inspired by Andrqej Tkacz, must be cleaned up and documented
"""


def pivot_taxonomic_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepares the taxonomic data (LSU and SSU tables) for analysis. Apart from
    pivoting, it also normalizes and calculates square root of the abundances.

    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.

    Returns:
        pd.DataFrame: A pivot table with taxonomic data.
    """
    # Select relevant columns
    df['taxonomic_concat'] = (
        df['ncbi_tax_id'].astype(str) + 
        ';sk_' + df['superkingdom'].fillna('') + 
        ';k_' + df['kingdom'].fillna('') + 
        ';p_' + df['phylum'].fillna('') + 
        ';c_' + df['class'].fillna('') + 
        ';o_' + df['order'].fillna('') + 
        ';f_' + df['family'].fillna('') + 
        ';g_' + df['genus'].fillna('') + 
        ';s_' + df['species'].fillna('')
    )
    pivot_table = df.pivot_table(
        index=['ncbi_tax_id','taxonomic_concat'], 
        columns='ref_code', 
        values='abundance',
    ).fillna(0)
    pivot_table = pivot_table.reset_index()
    # change inex name
    pivot_table.columns.name = None

    # normalize values
    pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: x / x.sum())
    pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: np.sqrt(x))
    return pivot_table


def separate_taxonomy(df):
    prokaryotes_keywords = ['Bacteria', 'Archaea']
    eukaryota_keywords = ['Discoba', 'Stramenopiles', 'Rhizaria', 'Alveolata', 'Amorphea', 'Archaeoplastida', 'Excavata']

    # Separate rows based on "Bacteria", "Archaea", and "Eukaryota" entries
    prokaryotes_all = df[df['taxonomic_concat'].str.contains("Bacteria|Archaea", regex=True)]
    eukaryota_all = df[df['taxonomic_concat'].str.contains("Eukaryota", regex=True)]

    # Save Prokaryotes and Eukaryota files
    prokaryotes_all.to_csv("Prokaryotes_all.csv")
    eukaryota_all.to_csv("Eukaryota_all.csv")

    # Further divide "Prokaryotes all" into "Bacteria" and "Archaea"
    bacteria = prokaryotes_all[prokaryotes_all['taxonomic_concat'].str.contains("Bacteria")]
    archaea = prokaryotes_all[prokaryotes_all['taxonomic_concat'].str.contains("Archaea")]

    # Save Bacteria and Archaea files
    bacteria.to_csv("Bacteria.csv")
    archaea.to_csv("Archaea.csv")

    # Further divide "Eukaryota all" by specific keywords
    eukaryota_dict = {}
    for keyword in eukaryota_keywords:
        subset = eukaryota_all[eukaryota_all['taxonomic_concat'].str.contains(keyword)]
        eukaryota_dict[keyword] = subset
        # Standardize each column to sum to 100 before saving the CSV
        subset_normalized = subset.div(subset.sum(axis=0), axis=1) * 100
        subset_normalized.to_csv(f"{keyword}.csv")


    # Apply taxonomy splitting to the index
    taxonomy_levels = bacteria['taxonomic_concat'].apply(split_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_levels.tolist(), columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
                               index=bacteria.index)

    # Combine taxonomy with the abundance data
    bacteria_data = pd.concat([taxonomy_df, bacteria], axis=1)

    return bacteria_data


# Function to extract taxonomic levels and remove prefixes before "Bacteria" or "Archaea"
def split_taxonomy(index_name):
    # Remove anything before "Bacteria" or "Archaea"
    if "Bacteria" in index_name:
        taxonomy = index_name.split("Bacteria;", 1)[1].split(";")
    elif "Archaea" in index_name:
        taxonomy = index_name.split("Archaea;", 1)[1].split(";")
    else:
        taxonomy = []
    # Return a list with taxonomic levels up to species
    return taxonomy[:7]  # ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


# Function to aggregate data by a specific taxonomic level
def aggregate_by_taxonomic_level(df, level):
    # Drop rows where the level is missing
    df_level = df.dropna(subset=[level])
    # Group by the specified level and sum abundances across samples (columns)
    df_grouped = df_level.groupby(level).sum(numeric_only=True)
    return df_grouped