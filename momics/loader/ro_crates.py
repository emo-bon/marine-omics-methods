import os
import requests
from typing import Dict


def get_rocrate_metadata_gh(sample_id: str) -> Dict:
    """
    Retrieves RO-Crate metadata from a GitHub repository.

    Args:
        sample_id (str): The ID of the sample.

    Returns:
        Dict: The metadata in JSON format.
    """
    url = f"https://api.github.com/repos/emo-bon/analysis-results-cluster-01-crate/contents/{sample_id}-ro-crate/ro-crate-metadata.json"
    req = requests.get(
        url,
        headers={
            "accept": "application/vnd.github.v3.raw",
        },
    )
    print("ro-crate-metadata.json request status", req.status_code)
    return req.json()


def get_rocrate_data(metadata_json: Dict, data_id: str):
    """
    Retrieves RO-Crate data file based on metadata.

    Args:
        metadata_json (Dict): The metadata in JSON format.
        data_id (str): The ID of the data file.
    Returns:
        str: The content of the data file.
    """
    raise NotImplementedError
