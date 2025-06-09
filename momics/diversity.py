import pandas as pd
import numpy as np
import logging
from typing import List, Dict

import skbio
from skbio.diversity import beta_diversity

# from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from sklearn.metrics import pairwise_distances

logging.getLogger(__name__)


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

    assert "source_mat_id" in filtered_metadata.columns, "The 'source_mat_id' column is missing from filtered metadata."

    # TODO: fix presence of source_mat_id in data first
    # assert "source_mat_id" in data.columns, "The 'source_mat_id' column is missing from data."

    # Match data and metadata samples
    abundance_matrix = data[filtered_metadata["ref_code"]].T

    permanova_results = {}
    # factors_to_test = permanova_additional_factors
    for remaining_factor in permanova_additional_factors:
        factor_metadata = filtered_metadata.dropna(subset=[remaining_factor])
        combined_abundance = abundance_matrix.loc[factor_metadata["ref_code"]]

        # Calculate Bray-Curtis distance matrix
        dissimilarity_matrix = pairwise_distances(
            combined_abundance, metric="braycurtis"
        )
        distance_matrix_obj = skbio.DistanceMatrix(
            dissimilarity_matrix, ids=combined_abundance.index
        )

        factor_metadata = factor_metadata.set_index("ref_code")
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
                    print(f"Factor: {remaining_factor}")
                    logging.info(f"Factor: {remaining_factor}")
                    logging.info(
                        f"  F-statistic: {permanova_result['test statistic']:.4f}"
                    )
                    logging.info(f"  p-value: {permanova_result['p-value']:.4f}\n")
        else:
            logging.info(
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

    # Create DataFrame with Shannon index and source_mat_id
    alpha_diversity_df = pd.DataFrame(
        {"source_mat_id": df["source_mat_id"], "Shannon": shannon_values}
    )

    # Merge with factors
    alpha_diversity_df = alpha_diversity_df.merge(factors, on="source_mat_id")

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
    """
    df_alpha_input = alpha_input(tables_dict, table_name).T.sort_values(by="ref_code")
    df_alpha_input = pd.merge(
        df_alpha_input, metadata, left_index=True, right_on="ref_code"
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
    """Update the subset indicator with the number of unique source_mat_ids."""
    indicator.value = df["source_mat_id"].nunique()


def update_taxa_count_indicator(indicator, df):
    """Update the taxa count indicator with the number of unique taxa."""
    indicator.value = len(df)


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
    # Convert DF
    out = pd.pivot_table(
        df,
        index="ref_code",
        columns=taxon,
        values="abundance",
        fill_value=0,
    )

    # Normalize rows
    if kind == "beta":
        out = out.div(out.sum(axis=1), axis=0)

    assert df.ncbi_tax_id.nunique(), out.shape[1]
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
        pd.DataFrame: A pivot table with species abundances indexed by the key column and ref_code as columns.
    """
    key_column = get_key_column(table_name)
    print("Key column:", key_column)

    # select distinct ref_codes from the dataframe
    ref_codes = tables_dict[table_name]["ref_code"].unique()
    print("length of the ref_codes:", len(ref_codes))
    out = pd.pivot_table(
        tables_dict[table_name],
        values="abundance",
        index=[key_column],
        columns=["ref_code"],
        aggfunc="sum",
        fill_value=0,
    )
    print("table shape:", out.shape)
    return out


# Example usage
# alpha_input = diversity_input(df, king='alpha')
# beta_input = diversity_input(df, king='beta')
