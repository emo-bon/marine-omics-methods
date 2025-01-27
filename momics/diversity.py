import pandas as pd
import numpy as np
from typing import List, Dict

from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa


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

    # Create DataFrame with Shannon index and ref_code
    alpha_diversity_df = pd.DataFrame(
        {"ref_code": df["ref_code"], "Shannon": shannon_values}
    )

    # Merge with factors
    alpha_diversity_df = alpha_diversity_df.merge(factors, on="ref_code")

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


# version with merging metadata
# def beta_diversity_parametrized(df, metadata):
#     # beta diversity
#     df_beta_input = diversity_input(df, kind='beta', taxon="ncbi_tax_id")

#     beta = beta_diversity("braycurtis", df_beta_input)

#     # merge metadata
#     df_beta = pd.merge(beta.to_data_frame(), metadata, left_index=True, right_on='ref_code')
#     return df_beta


# helper functions
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
# Example tables: ['go', 'go_slim', 'ips', 'ko', 'pfam']
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
