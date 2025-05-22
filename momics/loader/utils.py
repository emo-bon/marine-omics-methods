import io
import pandas as pd
from typing import Dict, List


def bytes_to_df(data: bytes, sep: str = "\t") -> Dict[str, pd.DataFrame]:
    """
    Convert a dictionary of bytes to a dictionary of DataFrames.

    Args:
        data (Dict[str, bytes]): A dictionary where keys are filenames and values are byte strings.
        sep (str): The separator used in the data files. Default is tab ("\t").

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are filenames and values are DataFrames.
    """
    df = pd.read_csv(
        io.StringIO(
            str(
                data,
                "utf-8",
            )
        ),
        sep=sep,
    )
    return df
