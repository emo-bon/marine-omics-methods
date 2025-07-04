import logging
import pandas as pd
import numpy as np

from tqdm import tqdm
from typing import List, Dict
from statsmodels.stats.multitest import multipletests
from skbio.stats import subsample_counts
from skbio.diversity import beta_diversity


# logger setup
FORMAT = "%(levelname)s | %(name)s | %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


"""
Some functions were originally developed by Andrzej Tkacz at CCMAR-Algarve.
"""


def pivot_taxonomic_data(
    df: pd.DataFrame, normalize: str = None, rarefy_depth: int = None
) -> pd.DataFrame:
    """
    Prepares the taxonomic data (LSU and SSU tables) for analysis. Apart from
    pivoting.

    Normalization of the pivot is optional. Methods include:
    - None: no normalization.
    - 'tss_sqrt': Total Sum Scaling and Square Root Transformation.
    - 'rarefy': rarefaction to a specified depth, if None, min of sample sums is used.

    TODO: refactor scaling to a new method and offer different options.

    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.
        normalize (str, optional): Normalization method.
            Options: None, 'tss_sqrt', 'rarefy'. Defaults to None.
        rarefy_depth (int, optional): Depth for rarefaction. If None, uses min sample sum.
            Defaults to None.

    Returns:
        pd.DataFrame: A pivot table with taxonomic data.
    """
    # Select relevant columns
    # fmt: off
    df1 = df.copy()
    df1["taxonomic_concat"] = (
        df1["ncbi_tax_id"].astype(str) +
        ";sk_" + df1["superkingdom"].fillna("") +
        ";k_" + df1["kingdom"].fillna("") +
        ";p_" + df1["phylum"].fillna("") +
        ";c_" + df1["class"].fillna("") +
        ";o_" + df1["order"].fillna("") +
        ";f_" + df1["family"].fillna("") +
        ";g_" + df1["genus"].fillna("") +
        ";s_" + df1["species"].fillna("")
    )
    # fmt: on
    pivot_table = (
        df1.pivot_table(
            index=["ncbi_tax_id", "taxonomic_concat"],
            columns="ref_code",
            values="abundance",
        )
        .fillna(0)
        .astype(int)
    )
    pivot_table = pivot_table.reset_index()
    # change inex name
    pivot_table.columns.name = None

    # normalize values
    if normalize == "tss_sqrt":
        # Total Sum Scaling and Square Root Transformation
        pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: x / x.sum())
        pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: np.sqrt(x))
    elif normalize == "rarefy":
        # rarefy
        pivot_table.iloc[:, 2:] = rarefy_table(
            pivot_table.iloc[:, 2:], depth=rarefy_depth
        )
    else:
        # no normalization
        pass
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


def remove_high_taxa(df: pd.DataFrame, taxonomy_ranks: list, tax_level: str = "phylum", strict: bool = True) -> pd.DataFrame:
    """
    Remove high level taxa from the dataframe.

    Args:
        df (pd.DataFrame): DataFrame containing taxonomic data.
        tax_level (str): The taxonomic level to filter by (e.g., 'phylum', 'class', 'order', etc.).
        strict (bool): If True, the lower taxa are all mapped to the tax_level.
            For instance, tax_level='phylum' will map all the more granular assignments (class, order, etc)
            to the phylum level.

    Returns:
        pd.DataFrame: DataFrame with rows where the specified taxonomic level is not None.
    """
    if tax_level not in df.columns:
        raise ValueError(f"Taxonomic level '{tax_level}' not found in DataFrame.")

    # Filter out rows where the taxonomic level is None or NaN
    df = df[~df[tax_level].isna()].copy()
    if strict:
        # unique values at the tax_level
        unique_taxa = df[tax_level].unique()
        bad_count = 0
        unmapped_taxa = []
        for taxon in tqdm(unique_taxa):
            tax_id = taxon_in_table(df, taxonomy_ranks, taxon, tax_level)
            # now map all the lower taxa to the tax_level
            if tax_id is None:
                # No lower taxon for no mapping needed.
                continue

            if tax_id == -1 and 'unclassified' in taxon.lower():
                # If the taxon is 'unclassified', we do not map it up.
                continue

            if tax_id == -1:
                # bad no mapping possible because I do not have the ncbi tax_id in the current table
                unmapped_taxa.append(taxon)
                bad_count += 1
            else:
                # mapping here
                df = map_taxa_up(df, taxon, tax_level, tax_id)

        logger.info(f'Number of bad taxa at {tax_level}: {bad_count}')
        logger.info(f'Unmapped taxa at {tax_level}: {unmapped_taxa}')
    return df[df[tax_level].notna()].copy()


def taxon_in_table(df: pd.DataFrame, taxonomy_ranks: list, taxon: str, tax_level: str) -> int:
    """
    Check if a taxon exists in the DataFrame at the specified taxonomic level.

    Args:
        df (pd.DataFrame): DataFrame containing taxonomic data.
        taxon (str): The taxon to check for.
        tax_level (str): The taxonomic level to check against.

    Returns:
        int: The index of the taxon in the DataFrame, or -1 if not found.
    """
    if tax_level not in df.columns:
        raise ValueError(f"Taxonomic level '{tax_level}' not found in DataFrame.")

    lower_taxon = taxonomy_ranks[taxonomy_ranks.index(tax_level) + 1] if taxonomy_ranks.index(tax_level) + 1 < len(taxonomy_ranks) else None
    if lower_taxon is None:
        return None

    # Find the indices of all the taxons
    index = df[df[tax_level] == taxon]['ncbi_tax_id']
    # check which one has lower taxon None
    for i in index.index:
        if pd.isna(df.loc[i, lower_taxon]):
            # Taxon {taxon} has lower taxon None at index i
            return df.loc[i, 'ncbi_tax_id']
    
    return -1  # Return -1 if the taxon is not found, therefore unknown ncbi_tax_id to map to


def map_taxa_up(df: pd.DataFrame, taxon: str, tax_level: str, tax_id: int) -> pd.DataFrame:
    """
    Map all lower taxa to the specified taxonomic level in the DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame containing taxonomic data.
        taxon (str): The taxon to map up.
        tax_level (str): The taxonomic level to map to.
        tax_id (int): The NCBI taxonomic ID to map to.

    Returns:
        pd.DataFrame: DataFrame with lower taxa mapped to the specified taxonomic level.
    """
    # Find the indices of all the taxons
    index = df[df[tax_level] == taxon]['ncbi_tax_id']
    index = index[index != tax_id]

    # Aggregate abundance for all rows matching the taxon at the specified taxonomic level, grouped by ref_code
    abundance_by_ref_code = df[df[tax_level] == taxon].groupby('ref_code')['abundance'].sum()

    for i in abundance_by_ref_code.index:  # ref_codes
        if i not in df['ref_code'].values:
            # ref_code i not in the DataFrame, skipping
            continue

        # Update the abundance for the taxon at the specified tax_id
        df.loc[(df['ncbi_tax_id'] == tax_id) & (df['ref_code'] == i), 'abundance'] = abundance_by_ref_code[i]

    # remove rows which are equal to index but not tax_id
    df = df[~df['ncbi_tax_id'].isin(index)]
    return df


def prevalence_cutoff(
    df: pd.DataFrame, percent: float = 10, skip_columns: int = 2
) -> pd.DataFrame:
    """
    Apply a prevalence cutoff to the DataFrame, removing features that do not
    appear in at least a certain percentage of samples.
    This is useful for filtering out low-prevalence features that may not be
    biologically relevant.

    Args:
        df (pd.DataFrame): The input DataFrame containing feature abundances.
        percent (float): The prevalence threshold as a percentage.
        skip_columns (int): The number of columns to skip (e.g., taxonomic information).

    Returns:
        pd.DataFrame: A filtered DataFrame with low-prevalence features removed.
    """
    # Calculate the number of samples
    num_samples = df.shape[1] - skip_columns
    # Calculate the prevalence threshold
    threshold = (percent / 100) * num_samples
    # Filter features based on prevalence
    filtered = df.loc[df.iloc[:, skip_columns:].gt(0).sum(axis=1) >= threshold]
    # Reset the index
    filtered = filtered.reset_index(drop=True)
    return filtered


def prevalence_cutoff_taxonomy(
    df: pd.DataFrame, percent: float = 10
    ) -> pd.DataFrame:
    """
    Apply a prevalence cutoff to the taxonomy DataFrame, which is not pivoted, removing
    features taxa with low abundance in each of the samples separately.

    Args:
        df (pd.DataFrame): The input DataFrame containing feature abundances.
        percent (float): The prevalence threshold as a percentage.

    Returns:
        pd.DataFrame: A filtered DataFrame with low-prevalence features removed.
    """

    try:
        del result_df
    except NameError:
        pass
    for ref_code in df['ref_code'].unique():
        abundance_sum = df[df['ref_code'] == ref_code]['abundance'].sum()
        threshold = abundance_sum * (percent / 100)

        # new filtered DataFrame
        df_filtered = df[(df['ref_code'] == ref_code) & (df['abundance'] > threshold)]

        # Concatenate each filtered DataFrame to a result DataFrame
        if 'result_df' not in locals():
            result_df = df_filtered.copy()
        else:
            result_df = pd.concat([result_df, df_filtered], ignore_index=True)
    # reset index
    result_df.reset_index(drop=True, inplace=True)
    return result_df


def rarefy_table(df: pd.DataFrame, depth: int = None, axis: int = 1) -> pd.DataFrame:
    """
    Rarefy an abundance table to a given depth. If depth is None, uses the
    minimum sample sum across all samples.
    This function is a wrapper around the skbio.stats.subsample_counts function.

    Args:
        df: pd.DataFrame (rows: features, columns: samples)
        depth: int or None, rarefaction depth. If None, uses min sample sum.
        axis: int, 1 for samples in columns, 0 for samples in rows.

    Returns:
        pd.DataFrame: A rarefied DataFrame.
    """
    if axis == 1:
        sample_sums = df.sum(axis=0)
    else:
        sample_sums = df.sum(axis=1)
    if depth is None:
        depth = sample_sums.min()
    rarefied = {}
    print("Minimum rarefaction depth:", depth)
    for sample in df.columns if axis == 1 else df.index:
        counts = df[sample].values if axis == 1 else df.loc[sample].values
        if sample_sums[sample] < depth:
            # Not enough counts, fill with NaN or zeros
            rarefied[sample] = np.full_like(counts, np.nan)
        else:
            rarefied_counts = subsample_counts(counts.astype(int), int(depth))
            rarefied[sample] = rarefied_counts
    if axis == 1:
        return pd.DataFrame(rarefied, index=df.index)
    else:
        return pd.DataFrame(rarefied, columns=df.columns)



def fill_taxonomy_placeholders(df: pd.DataFrame, taxonomy_ranks: list) -> pd.DataFrame:
    """
    Fill higher missing taxonomy levels in a DataFrame with placeholders
    like 'unclassified_<lower_rank_value>'.

    No downwards propagation is done, only upwards filling.

    Parameters:
    - df: pandas DataFrame containing taxonomy columns.
    - taxonomy_ranks: ordered list of taxonomy column names from higher to lower rank.

    Returns:
    - df with placeholders filled.
    """
    df = df.copy()

    for i in range(2, len(taxonomy_ranks)):
        lower = taxonomy_ranks[- i + 1]  # lower rank column
        current = taxonomy_ranks[-i]
        
        # Fill missing current-level values using higher-level information
        df[current] = df.apply(
            lambda row: f'unclassified_{row[lower]}' if row[current] == '' else row[current],
            axis=1
        )

    return df


def split_metadata(metadata: pd.DataFrame, factor: str) -> Dict[str, list]:
    """
    Splits the metadata ref codes to dictionary of key being the factor value and
    value is a list of the ref codes.

    Args:
        metadata (pd.DataFrame): The input DataFrame containing metadata.
        factor (str): The column name to split the metadata by.
    Returns:
        Dict[str, list]: A dictionary with keys as unique values of the factor and
                         values as lists of ref codes.
    """
    if factor not in metadata.columns:
        raise ValueError(f"Factor '{factor}' not found in DataFrame columns.")
    # check if column is categorical
    if not isinstance(metadata[factor].dtype, pd.CategoricalDtype):
        raise ValueError(f"Column '{factor}' is not categorical (object dtype).")

    # for each unique value in the factor column, create a new table and append to a dictionary
    grouped_data = {}
    for value in metadata[factor].unique():
        # filter the dataframe
        filtered_df = metadata[metadata[factor] == value]
        # get the ref codes
        ref_codes = filtered_df["ref_code"].tolist()
        # add to the dictionary
        grouped_data[value] = ref_codes

    return grouped_data


def split_taxonomic_data(
    taxonomy: pd.DataFrame, groups: Dict[str, list]
) -> Dict[str, pd.DataFrame]:
    """
    Splits the taxonomic data into dictionary of DataFrames for each group.
    The split is based on the ref_code column, which needs to be present in the
    dataframes.

    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.
        groups (Dict[str, list]): A dictionary where keys are unique values of a factor
            which correspond to the groups to split by.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary with keys as unique values of the factor and
            values as DataFrames containing taxonomic data for each group.
    """
    if not isinstance(groups, dict):
        raise ValueError("Groups must be a dictionary.")

    # for each unique value in the factor column, create a new table and append to a dictionary
    grouped_data = {}
    for value, ref_codes in groups.items():
        # filter the dataframe
        filtered_df = taxonomy[taxonomy["ref_code"].isin(ref_codes)]
        # add to the dictionary
        grouped_data[value] = filtered_df

    return grouped_data


def split_taxonomic_data_pivoted(
    taxonomy: pd.DataFrame, groups: Dict[str, list]
) -> Dict[str, pd.DataFrame]:
    """
    Splits the taxonomic data into dictionary of DataFrames for each group.
    The split is based on the column names which need to match between the taxonomy DataFrame
    and the groups lists. The DataFrame should have a 'ncbi_tax_id' and 'taxonomic_concat' which
    will serve as index of the resulting DataFrames.

    Args:
        taxonomy (pd.DataFrame): The input DataFrame containing taxonomic information.
        groups (Dict[str, list]): A dictionary where keys are unique values of a factor
            which correspond to the groups to split by.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary with keys as unique values of the factor and
            values as DataFrames with separate columns for each taxonomic rank.
    """
    if not isinstance(groups, dict):
        raise ValueError("Groups must be a dictionary.")

    # for each unique value in the factor column, create a new table and append to a dictionary
    grouped_data = {}
    for value, ref_codes in groups.items():
        # filter the dataframe
        filtered_df = taxonomy[["ncbi_tax_id", "taxonomic_concat"] + ref_codes]
        # remove rows with all zeros and print how many rows were removed
        len_before = len(filtered_df)
        filtered_df = filtered_df[filtered_df[ref_codes].sum(axis=1) != 0]
        len_after = len(filtered_df)
        print(f"Removed {len_before - len_after} rows with all zeros for {value}.")
        # check if the dataframe is empty
        if filtered_df.empty:
            print(f"Warning: No data for {value} in the taxonomic data.")
            continue

        # add to the dictionary
        grouped_data[value] = filtered_df

    return grouped_data


def compute_bray_curtis(
    df: pd.DataFrame, skip_cols: int = 2, direction: str = "samples"
) -> pd.DataFrame:
    """
    Compute Bray-Curtis dissimilarity and return as a pandas DataFrame.
    This function computes the Bray-Curtis dissimilarity for samples in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing sample counts.
        skip_cols (int): Number of columns to skip (e.g., taxonomic information).
        direction (str): Direction of the dissimilarity calculation, 'samples' or 'taxa'.

    Returns:
        pd.DataFrame: A DataFrame containing the Bray-Curtis dissimilarity matrix.
    """
    if direction not in ["samples", "taxa"]:
        raise ValueError("Direction must be either 'samples' or 'taxa'.")

    if direction == "samples":
        # Use the sample IDs as the index
        ids = df.columns[skip_cols:].astype(str).tolist()
        result = beta_diversity(
            metric="braycurtis", counts=df.iloc[:, skip_cols:].T, ids=ids
        )
    elif direction == "taxa":
        ids = df["ncbi_tax_id"].astype(str).tolist()
        result = beta_diversity(
            metric="braycurtis", counts=df.iloc[:, skip_cols:], ids=ids
        )

    bray_curtis_df = pd.DataFrame(result.data, index=ids, columns=ids)
    return bray_curtis_df


def fdr_pvals(p_spearman_df: pd.DataFrame, pval_cutoff: float) -> pd.DataFrame:
    """
    Apply FDR correction to the p-values DataFrame using Benjamini/Hochberg (non-negative)
    method. This function extracts the upper triangle of the p-values DataFrame.

    Args:
        p_spearman_df (pd.DataFrame): DataFrame containing p-values.
        pval_cutoff (float): P-value cutoff for FDR correction.

    Returns:
        pd.DataFrame: DataFrame with FDR corrected p-values.
    """
    # Extract upper triangle p-values
    pval_array = (
        p_spearman_df.where(np.triu(np.ones(p_spearman_df.shape), k=1).astype(bool))
        .stack()
        .values
    )

    # Apply FDR correction
    _rejected, pvals_corrected, _, _ = multipletests(
        pval_array, alpha=pval_cutoff, method="fdr_bh"
    )

    # Map corrected p-values back to a DataFrame
    pvals_fdr = p_spearman_df.copy()
    pvals_fdr.values[np.triu_indices_from(p_spearman_df, k=1)] = pvals_corrected
    pvals_fdr.values[np.tril_indices_from(p_spearman_df, k=0)] = (
        np.nan
    )  # Optional: keep only upper triangle
    return pvals_fdr
