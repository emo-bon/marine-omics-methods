import logging
import pandas as pd
import numpy as np
from typing import Union, List, Dict

import skbio
from skbio.diversity import beta_diversity

from skbio.stats.distance import permanova
from sklearn.metrics import pairwise_distances
from .utils import (
    check_index_names,
)
from momics.constants import TAXONOMY_RANKS

# logger setup
FORMAT = "%(levelname)s | %(name)s | %(message)s"
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


#########################
# Statistical functions #
#########################
def run_permanova(
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    permanova_factor: str,
    permanova_group: List[str],
    permanova_additional_factors: List[str],
    permutations: int = 999,
    verbose: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Run PERMANOVA on the given data and metadata.
    Args:
        data (pd.DataFrame): DataFrame containing the abundance data.
        metadata (pd.DataFrame): DataFrame containing the metadata.
        permanova_factor (str): The factor to use for PERMANOVA.
        permanova_group (List[str]): List of groups to include in the analysis.
        permanova_additional_factors (List[str]): Additional factors to test.
        permutations (int): Number of permutations for PERMANOVA. Default is 999.
        verbose (bool): If True, print detailed output.
    Returns:
        Dict[str, pd.DataFrame]: Dictionary containing PERMANOVA results for each factor.
    """
    # Filter metadata based on selected groups
    if permanova_factor == "All":
        filtered_metadata = metadata.copy()
    else:
        filtered_metadata = metadata[metadata[permanova_factor].isin(permanova_group)]

    # Match data and metadata samples
    abundance_matrix = data[filtered_metadata.index].T

    permanova_results = {}
    # factors_to_test = permanova_additional_factors
    for remaining_factor in permanova_additional_factors:
        factor_metadata = filtered_metadata.dropna(subset=[remaining_factor])
        combined_abundance = abundance_matrix.loc[factor_metadata.index]

        # Calculate Bray-Curtis distance matrix
        dissimilarity_matrix = pairwise_distances(
            combined_abundance, metric="braycurtis"
        )
        distance_matrix_obj = skbio.DistanceMatrix(
            dissimilarity_matrix, ids=combined_abundance.index
        )

        factor_metadata = factor_metadata.loc[
            factor_metadata.index.intersection(distance_matrix_obj.ids)
        ]

        if remaining_factor not in factor_metadata.columns:
            continue

        group_vector = factor_metadata[remaining_factor]
        if group_vector.nunique() < len(group_vector):
            if set(distance_matrix_obj.ids) == set(group_vector.index):
                permanova_result = permanova(
                    distance_matrix_obj,
                    grouping=group_vector,
                    permutations=permutations,
                )
                permanova_results[remaining_factor] = permanova_result
                if verbose:
                    logger.info(f"Factor: {remaining_factor}")
                    logger.info(
                        f"  F-statistic: {permanova_result['test statistic']:.4f}"
                    )
                    logger.info(f"  p-value: {permanova_result['p-value']:.4f}\n")
        else:
            logger.info(
                f"Skipping factor '{remaining_factor}' due to unique values in grouping vector."
            )

    return permanova_results


def shannon_index(row: pd.Series) -> float:
    """
    Calculates the Shannon index for a given row of data.

    Args:
        row (pd.Series): A row of data containing species abundances.

    Returns:
        float: The Shannon index value.
    """
    row = pd.to_numeric(row, errors="coerce")
    total_abundance = row.sum()
    if total_abundance == 0:
        return np.nan
    relative_abundance = row / total_abundance
    ln_relative_abundance = np.log(relative_abundance)
    ln_relative_abundance[relative_abundance == 0] = 0
    multi = relative_abundance * ln_relative_abundance * -1
    return multi.sum()  # Shannon entropy


def calculate_shannon_index(df: pd.DataFrame) -> pd.Series:
    """
    Applies the Shannon index calculation to each row of a DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing species abundances.

    Returns:
        pd.Series: A Series containing the Shannon index for each row.
    """
    return df.apply(shannon_index, axis=1)


####################
# Search functions #
####################
def find_taxa_in_table(
        table: pd.DataFrame,
        tax_level: str,
        search_term: Union[str, int],
        ncbi_tax_id: bool=False,
        exact_match:bool=False,
    ) -> pd.DataFrame:
    """
    Find taxa in the given table at the specified taxonomic level matching the search term.

    args:
        table (pd.DataFrame): DataFrame containing taxonomic data.
        tax_level (str): Taxonomic level to search ('all' for all levels).
        search_term (str|int): Term to search for.
        ncbi_tax_id (bool): If True, search by NCBI taxonomic ID.
        exact_match (bool): If True, perform exact match; otherwise, use substring match.

    returns:
        pd.DataFrame: DataFrame containing matching taxa.
    """
    # ncbi_tax_id search
    index_names = getattr(table.index, "names", [])
    if ncbi_tax_id and ('ncbi_tax_id' not in table.columns and 'ncbi_tax_id' not in index_names):
        raise ValueError("The table does not contain 'ncbi_tax_id' column or index level.")

    # if ncbi_tax_id is an index level, bring it into a column for uniform handling
    if ncbi_tax_id and ('ncbi_tax_id' in index_names):
        table = table.reset_index()

    if ncbi_tax_id:
        # Search by NCBI taxonomic ID
        matching_taxa = table[table['ncbi_tax_id'].astype(str) == str(search_term)]
        return matching_taxa.set_index(index_names) if index_names else matching_taxa

    # search by taxonomic level, all ranks
    if tax_level == 'all':
        found = []
        for tax_level in TAXONOMY_RANKS:
            if exact_match:
                found.append(table[table[tax_level].str.lower().fillna('') == search_term.lower()])
            else:
                found.append(table[table[tax_level].str.contains(search_term, case=False, na=False)])
        matching_taxa = pd.concat(found)
    # specific taxonomic level
    else:
        if exact_match:
            matching_taxa = table[table[tax_level].str.lower().fillna('') == search_term.lower()]
        else:
            matching_taxa = table[table[tax_level].str.contains(search_term, case=False, na=False)]

    return matching_taxa

#######################
# diversity functions #
#######################

def calculate_alpha_diversity(df: pd.DataFrame, factors: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the alpha diversity (Shannon index) for a DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing species abundances.
        factors (pd.DataFrame): A DataFrame containing additional factors to merge.

    Returns:
        pd.DataFrame: A DataFrame containing the Shannon index and additional factors.
    """
    # Select columns that start with the appropriate prefix
    numeric_columns = [
        col
        for col in df.columns
        if col.startswith("GO:")
        or col.startswith("IPR")
        or col.startswith("K")
        or col.startswith("PF")
    ]

    # Calculate Shannon index only from the selected columns
    shannon_values = calculate_shannon_index(df[numeric_columns])

    # Create DataFrame with Shannon values and index of the input DataFrame
    alpha_diversity_df = pd.DataFrame(
        {df.index.name: df.index, "Shannon": shannon_values},
    )
    alpha_diversity_df.set_index(df.index.name, inplace=True)

    # Merge with factors
    alpha_diversity_df = alpha_diversity_df.merge(
        factors, left_index=True, right_index=True
    )

    return alpha_diversity_df


# alpha diversity
def alpha_diversity_parametrized(
    tables_dict: Dict[str, pd.DataFrame], table_name: str, metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculates the alpha diversity for a list of tables and merges with metadata.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table.
        metadata (pd.DataFrame): A DataFrame containing metadata.

    Returns:
        pd.DataFrame: A DataFrame containing the alpha diversity and metadata.

    Raises:
        ValueError: If the index names of the input DataFrame and metadata do not match.
    """
    df_alpha_input = alpha_input(tables_dict, table_name).T.sort_index()

    # Ensure the index name is set correctly
    if not check_index_names(df_alpha_input, metadata):
        raise ValueError(
            "The index names of the input DataFrame and metadata do not match."
        )

    # Merge the alpha input DataFrame with metadata
    df_alpha_input = pd.merge(
        df_alpha_input,
        metadata,
        left_index=True,
        right_index=True,
    )
    alpha = calculate_alpha_diversity(df_alpha_input, metadata)
    return alpha


def beta_diversity_parametrized(
    df: pd.DataFrame, taxon: str, metric: str = "braycurtis"
) -> pd.DataFrame:
    """
    Calculates the beta diversity for a DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing species abundances.
        taxon (str): The taxon to use for the beta diversity calculation.
        metric (str, optional): The distance metric to use. Defaults to "braycurtis".

    Returns:
        pd.DataFrame: A DataFrame containing the beta diversity distances.
    """
    df_beta_input = diversity_input(df, kind="beta", taxon=taxon)
    beta = beta_diversity(metric, df_beta_input)
    return beta


####################
# helper functions #
####################
def update_subset_indicator(indicator, df):
    """
    Update the subset indicator with the number of unique `index` ids.
    """
    indicator.value = df.index.nunique()


def update_taxa_count_indicator(indicator, df):
    """
    Update the taxa count indicator with the number of unique taxa.
    """
    indicator.value = df.index.nunique()


# I think this is only useful for beta, not alpha diversity
def diversity_input(
    df: pd.DataFrame, kind: str = "alpha", taxon: str = "ncbi_tax_id"
) -> pd.DataFrame:
    """
    Prepare input for diversity analysis.

    Args:
        df (pd.DataFrame): The input dataframe.
        kind (str): The type of diversity analysis. Either 'alpha' or 'beta'.
        taxon (str): The column name containing the taxon IDs.

    Returns:
        pd.DataFrame: The input for diversity analysis.
    """
    if isinstance(df.index, pd.MultiIndex):
        index_name = df.index.names[0]
    else:
        index_name = df.index.name
    df1 = df.reset_index()

    # Convert DF
    out = pd.pivot_table(
        df1,
        index=index_name,
        columns=taxon,
        values="abundance",
        fill_value=0,
    )

    # Normalize rows
    if kind == "beta":
        out = out.div(out.sum(axis=1), axis=0)

    assert df1[taxon].nunique(), out.shape[1]
    return out


# Function to get the appropriate column based on the selected table
# Valid table names: ['go', 'go_slim', 'ips', 'ko', 'pfam']
def get_key_column(table_name: str) -> str:
    """Returns the key column name based on the table name.

    Args:
        table_name (str): The name of the table.

    Returns:
        str: The key column name.

    Raises:
        ValueError: If the table name is unknown.
    """
    if table_name in ["go", "go_slim"]:
        return "id"
    elif table_name == "ips":
        return "accession"
    elif table_name in ["ko", "pfam"]:
        return "entry"
    else:
        raise ValueError(f"Unknown table: {table_name}")


def alpha_input(tables_dict: Dict[str, pd.DataFrame], table_name: str) -> pd.DataFrame:
    """
    Prepares the input data for alpha diversity calculation.

    Args:
        tables_dict (Dict[str, pd.DataFrame]): A dictionary of DataFrames containing species abundances.
        table_name (str): The name of the table to process.

    Returns:
        pd.DataFrame: A pivot table with species abundances indexed by the key column of the functional table
            and index column converted to columns.
    """
    key_column = get_key_column(table_name)

    df = tables_dict[table_name]
    index_name = df.index.name
    df = df.reset_index()

    # select distinct index_vals from the dataframe
    out = pd.pivot_table(
        df,
        values="abundance",
        index=[key_column],
        columns=[index_name],
        aggfunc="sum",
        fill_value=0,
    )
    return out
