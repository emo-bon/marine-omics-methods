import os
import pandas as pd


def load_parquet_files(folder):
    """
    Loads all .parquet files in folder and stores them in a dictionary.
    The keys of the dictionary are the file names without the .parquet
    extension. If filename contains more '.' then only last part of the
    names is included.

    Example: metagoflow_analyses.go_slim.parquet -> key = "go_slim"

    Parameters
    ----------
    folder : str
        The path to the folder containing the .parquet files.

    Returns
    -------
    mgf_parquet_dfs : dict
        A dictionary containing the data frames of the .parquet files.

    """
    # Create an empty dictionary to store the data frames
    # In this disctionary the data tables will be stored as pandas data frames
    mgf_parquet_dfs = {}

    # Loop through the folder and load each .parquet file
    for file_name in os.listdir(folder):
        if file_name.endswith(".parquet"):
            file_path = os.path.join("parquet_files", file_name)
            # Load the parquet file into a DataFrame
            df = pd.read_parquet(file_path)
            # Use the file name without extension as the dictionary key
            name = file_name.split(".")[-2]
            mgf_parquet_dfs[name] = df
    return mgf_parquet_dfs
