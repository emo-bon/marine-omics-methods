import os
import pandas as pd
from typing import Dict
from mgo.udal import UDAL


def load_parquets(folder: str) -> Dict[str, pd.DataFrame]:
    """
    Loads all .parquet files in a folder and stores them in a dictionary.

    The keys of the dictionary are the file names without the .parquet
    extension. If the filename contains more than one '.', only the last part
    of the name is included.

    Example:
        metagoflow_analyses.go_slim.parquet -> key = "go_slim"

    Args:
        folder (str): The path to the folder containing the .parquet files.

    Returns:
        dict: A dictionary containing the data frames of the .parquet files.
    """
    # Create an empty dictionary to store the data frames
    # In this disctionary the data tables will be stored as pandas data frames
    mgf_parquet_dfs = {}

    # Loop through the folder and load each .parquet file
    for file_name in os.listdir(folder):
        if file_name.endswith(".parquet"):
            file_path = os.path.join(folder, file_name)
            # Load the parquet file into a DataFrame
            df = pd.read_parquet(file_path)
            # Use the file name without extension as the dictionary key
            name = file_name.split(".")[-2].lower()
            mgf_parquet_dfs[name] = df
    return mgf_parquet_dfs


def load_parquets_udal():
    """
    Load parquet files into a dictionary by looping udal calls
    """
    udal = UDAL()

    parquets = {}
    for dataset in ["go", "go_slim", "ips", "ko", "pfam", "lsu", "ssu"]:
        parquets[dataset] = udal.execute(f"urn:embrc.eu:emobon:{dataset}").data()
    return parquets
