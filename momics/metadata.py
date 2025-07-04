"""
Methods to manipulate, concatenate, merge and enrich metadata files from EMO-BON samplings.

Some of these methods work as temporary solution to bad or incomplete data validation of the metadata tables.

Hopefully, that will not be the case for ever.
"""

import os
import pandas as pd
from typing import Dict, List
from datetime import datetime
from mgo.udal import UDAL


def get_metadata(folder):
    # Load metadata
    sample_metadata = pd.read_csv(
        os.path.join(folder, "Batch1and2_combined_logsheets_2024-11-12.csv")
    )
    assert 'source_mat_id' in sample_metadata.columns, (
        "The sample metadata file does not contain the 'source_mat_id' column."
    )
    print(sample_metadata['source_mat_id'].unique())

    observatory_metadata = pd.read_csv(
        os.path.join(folder, "Observatory_combined_logsheets_validated.csv")
    )

    # Merge metadata
    full_metadata = pd.merge(
        sample_metadata,
        observatory_metadata,
        on=["obs_id", "env_package"],  # Matching conditions
        how="inner",  # Inner join
    )

    # Sort the merged dataframe by 'ref_code' column in ascending order
    full_metadata = full_metadata.sort_values(by="ref_code", ascending=True)

    # first convert some of the boolean cols
    full_metadata["failure"] = full_metadata["failure"].astype(str)
    # replace the 'nan' values with 'NA'
    full_metadata["failure"] = full_metadata["failure"].replace("nan", "NA")

    # adding replacement for the missing values for object type columns
    full_metadata = fill_na_for_object_columns(full_metadata)

    return full_metadata


def get_metadata_udal():
    """
    Load metadata from the UDAL API
    """
    udal = UDAL()

    sample_metadata = udal.execute("urn:embrc.eu:emobon:logsheets").data().reset_index()

    assert 'source_mat_id' in sample_metadata.columns, (
        "The sample metadata file does not contain the 'source_mat_id' column."
    )

    observatory_metadata = (
        udal.execute("urn:embrc.eu:emobon:observatories").data().set_index("obs_id")
    )

    # Merge metadata
    full_metadata = pd.merge(
        sample_metadata,
        observatory_metadata,
        on=["obs_id", "env_package"],  # Matching conditions
        how="inner",  # Inner join
    )

    # Sort the merged dataframe by 'ref_code' column in ascending order
    full_metadata = full_metadata.sort_values(by="ref_code", ascending=True)

    # first convert some of the boolean cols
    full_metadata["failure"] = full_metadata["failure"].astype(str)
    # replace the 'nan' values with 'NA'
    full_metadata["failure"] = full_metadata["failure"].replace("nan", "NA")

    # adding replacement for the missing values for object type columns
    full_metadata = fill_na_for_object_columns(full_metadata)

    return full_metadata


#####################
## Filter metadata ##
#####################
# Filter the metadata table based on the selections in box_granular
def filter_metadata_table(
    metadata_df: pd.DataFrame, selected_factors: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Filter the metadata DataFrame based on selected factors and their values.

    Args:
        metadata_df (pd.DataFrame): The metadata DataFrame to filter.
        selected_factors (Dict[str, List[str]]): A dictionary where keys are factor names and values are lists of selected values.
            If 'All' is in the list, that factor will not be filtered.
    Returns:
        pd.DataFrame: The filtered metadata DataFrame.
    """
    # Create a copy of the metadata DataFrame
    filtered_metadata = metadata_df.copy()
    # Apply filters for each selected factor
    for factor, selected_values in selected_factors.items():
        if "All" not in selected_values:
            filtered_metadata = filtered_metadata[
                filtered_metadata[factor].isin(selected_values)
            ]
    return filtered_metadata


## filter data according to the metadata
def filter_data(df: pd.DataFrame, filtered_metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Filter the DataFrame based on the filtered metadata.
    This function filters the DataFrame columns based on the 'ref_code' values in the filtered metadata.

    Args:
        df (pd.DataFrame): The DataFrame to filter.
        filtered_metadata (pd.DataFrame): The filtered metadata DataFrame.
    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    # filter columns names of df which are in the filtered metadata

    assert "source_mat_id" in filtered_metadata.columns, (
        "The filtered metadata does not contain the 'source_mat_id' column."
    )

    cols_to_keep = list(
        [
            col
            for col in df.columns.str.strip()
            if col in filtered_metadata["ref_code"].to_list()
        ]
    )

    filtered_df = df[cols_to_keep]
    return filtered_df


######################
## Enhance metadata ##
######################
def enhance_metadata(
    metadata: pd.DataFrame, df_validation: pd.DataFrame
) -> pd.DataFrame:
    """
    Enhance the metadata DataFrame by processing the 'collection_date' column and extracting the season.
    This function also filters the metadata based on the 'ref_code' values in the df_validation DataFrame.

    Args:
        metadata (pd.DataFrame): The metadata DataFrame to enhance.
        df_validation (pd.DataFrame): The DataFrame containing valid samples for filtering.
    Returns:
        pd.DataFrame: The enhanced metadata DataFrame.
    """
    metadata = process_collection_date(metadata)
    metadata = extract_season(metadata)

    # Filter the metadata on the 'ref_code' only for entries that are in df_valid
    metadata = metadata[metadata["ref_code"].isin(df_validation["ref_code"])]

    missing = df_validation[~df_validation["ref_code"].isin(metadata["ref_code"])]
    assert len(missing) == 0, "Missing samples in the metadata"
    assert len(metadata) == len(
        df_validation
    ), "Filtered metadata does not match the valid samples"

    # add column to identify properly the replicates
    metadata["replicate_info"] = (
        metadata["obs_id"].astype(str)
        + "_"
        + metadata["env_package"].astype(str)
        + "_"
        + metadata["collection_date"].astype(str)
        + "_"
        + metadata["size_frac"].astype(str)
    )

    return metadata


def process_collection_date(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Process the 'collection_date' column in the metadata DataFrame.
    This function converts the 'collection_date' column to datetime format,
    extracts the year, month, and day, and adds them as new columns.
    It also converts the month number to the month name (abbreviated).

    Args:
        metadata (pd.DataFrame): The metadata DataFrame containing the 'collection_date' column.

    Returns:
        pd.DataFrame: The updated metadata DataFrame with new columns for year, month, and day.
    """
    # Convert the 'collection_date' column to datetime
    metadata["collection_date"] = metadata["collection_date"].apply(
        lambda x: datetime.strptime(x, "%Y-%m-%d") if x is not None else None
    )
    # Extract the year from the 'collection_date' column
    metadata["year"] = metadata["collection_date"].apply(
        lambda x: x.year if x is not None else None
    )
    # Extract the month from the 'collection_date' column
    metadata["month"] = metadata["collection_date"].apply(
        lambda x: x.month if x is not None else None
    )
    # Convert month to month name
    metadata["month_name"] = metadata["month"].apply(
        lambda x: (
            datetime.strptime(str(x), "%m").strftime("%B")[:3]
            if x is not None
            else None
        )
    )
    # Extract the day from the 'collection_date' column
    metadata["day"] = metadata["collection_date"].apply(
        lambda x: x.day if x is not None else None
    )
    return metadata


def extract_season(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Add a 'season' column to the metadata DataFrame.
    This function determines the season based on the 'month' and 'day' columns
    and adds it as a new column to the DataFrame.

    Args:
        metadata (pd.DataFrame): The metadata DataFrame containing 'month' and 'day' columns.

    Returns:
        pd.DataFrame: The updated metadata DataFrame with a new 'season' column.
    """
    # Extract the season based on the month and day
    metadata["season"] = metadata.apply(extract_season_single, axis=1)
    return metadata


def extract_season_single(row):
    """
    Determine the season based on the month and day.
    This function is used as a helper function for the apply method."
    """
    if (
        (row["month"] == 3 and row["day"] >= 21)
        or (row["month"] == 4)
        or (row["month"] == 5)
        or (row["month"] == 6 and row["day"] < 21)
    ):
        return "Spring"
    elif (
        (row["month"] == 6 and row["day"] >= 21)
        or (row["month"] == 7)
        or (row["month"] == 8)
        or (row["month"] == 9 and row["day"] < 23)
    ):
        return "Summer"
    elif (
        (row["month"] == 9 and row["day"] >= 23)
        or (row["month"] == 10)
        or (row["month"] == 11)
        or (row["month"] == 12 and row["day"] < 21)
    ):
        return "Autumn"
    else:  # Winter
        return "Winter"


def fill_na_for_object_columns(df):
    """
    Fill NA values with 'NA' for object columns in the dataframe.

    Args:
        df (pd.DataFrame): The input dataframe.

    Returns:
        pd.DataFrame: The dataframe with NA values filled for object columns.
    """
    # Apply fillna only to object columns
    df[df.select_dtypes(include=["object"]).columns] = df.select_dtypes(
        include=["object"]
    ).apply(lambda col: col.fillna("NA"))
    return df
