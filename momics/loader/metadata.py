"""
Methods to manipulate, concatenate, merge and enrich metadata files from EMO-BON samplings.

Some of these methods work as temporary solution to bad or incomplete data validation of the metadata tables.

Hopefully, that will not be the case for ever.
"""

import os
import pandas as pd
from typing import Dict, List, Tuple
from datetime import datetime


######################
## Enhance metadata ##
######################
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
    metadata['collection_date'] = metadata['collection_date'].apply(
        lambda x: datetime.strptime(x, "%Y-%m-%d") if x is not None else None
    )
    # Extract the year from the 'collection_date' column
    metadata['year'] = metadata['collection_date'].apply(
        lambda x: x.year if x is not None else None
    )
    # Extract the month from the 'collection_date' column
    metadata['month'] = metadata['collection_date'].apply(
        lambda x: x.month if x is not None else None
    )
    # Convert month to month name
    metadata['month_name'] = metadata['month'].apply(
        lambda x: datetime.strptime(str(x), "%m").strftime("%B")[:3] if x is not None else None
    )
    # Extract the day from the 'collection_date' column
    metadata['day'] = metadata['collection_date'].apply(
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
    metadata['season'] = metadata.apply(extract_season_single, axis=1)
    return metadata


def extract_season_single(row):
    """
    Determine the season based on the month and day.
    This function is used as a helper function for the apply method."
    """
    if (row['month'] == 3 and row['day'] >= 21) or (row['month'] == 4) or (row['month'] == 5) or (row['month'] == 6 and row['day'] < 21):
        return 'Spring'
    elif (row['month'] == 6 and row['day'] >= 21) or (row['month'] == 7) or (row['month'] == 8) or (row['month'] == 9 and row['day'] < 23):
        return 'Summer'
    elif (row['month'] == 9 and row['day'] >= 23) or (row['month'] == 10) or (row['month'] == 11) or (row['month'] == 12 and row['day'] < 21):
        return 'Autumn'
    else:  # Winter
        return 'Winter'
    
