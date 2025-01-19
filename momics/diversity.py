import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa


# this should go to stats, or use skbio method
def shannon_index(row):
    row = pd.to_numeric(row, errors="coerce")
    total_abundance = row.sum()
    if total_abundance == 0:
        return np.nan
    relative_abundance = row / total_abundance
    ln_relative_abundance = np.log(relative_abundance)
    ln_relative_abundance[relative_abundance == 0] = 0
    multi = relative_abundance * ln_relative_abundance * -1
    return multi.sum()  # Shannon entropy


def calculate_shannon_index(df):
    return df.apply(shannon_index, axis=1)


#######################
# diversity functions #
#######################
def calculate_alpha_diversity(df, factors):
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
def alpha_diversity_parametrized(table_list, table_name, metadata):
    df_alpha_input = alpha_input(table_list, table_name).T.sort_values(by="ref_code")
    df_alpha_input = pd.merge(df_alpha_input, metadata, left_index=True, right_on='ref_code')
    alpha = calculate_alpha_diversity(df_alpha_input, metadata)
    return alpha


def beta_diversity_parametrized(df, taxon, metric="braycurtis"):
    df_beta_input = diversity_input(df, kind='beta', taxon=taxon)
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
def diversity_input(df, kind='alpha', taxon="ncbi_tax_id"):
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
    if kind == 'beta':
        out = out.div(out.sum(axis=1), axis=0)

    assert df.ncbi_tax_id.nunique(), out.shape[1]
    return out


# Function to get the appropriate column based on the selected table
# Example tables: ['go', 'go_slim', 'ips', 'ko', 'pfam']
def get_key_column(table_name):
    if table_name in ["go", "go_slim"]:
        return "id"
    elif table_name == "ips":
        return "accession"
    elif table_name in ["ko", "pfam"]:
        return "entry"
    else:
        raise ValueError(f"Unknown table: {table_name}")
    

def alpha_input(table_list, table_name):
    key_column = get_key_column(table_name)
    print("Key column:", key_column)

    # select distinct ref_codes from the dataframe
    ref_codes = table_list[table_name]['ref_code'].unique()
    print('length of the ref_codes:', len(ref_codes))
    out = pd.pivot_table(table_list[table_name],
                         values='abundance',
                         index=[key_column],
                         columns=['ref_code'],
                         aggfunc='sum',
                         fill_value=0,
                         )
    print('table shape:', out.shape)
    return out
# Example usage
# alpha_input = diversity_input(df, king='alpha')
# beta_input = diversity_input(df, king='beta')
