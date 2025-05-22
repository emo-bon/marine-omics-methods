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


def extract_data_by_name(metadata: Dict, name: str) -> Dict:
    """
    Extracts data from metadata based on the name.
    Args:
        metadata (Dict): The metadata in JSON format.
        name (str): The name of the data to extract.
    Returns:
        Dict: The extracted data.
    """
    for d in metadata["@graph"]:
        if "name" in d.keys() and d["name"] == name:
            return d
    return None


def extract_all_datafiles(metadata: Dict) -> list:
    """
    Extracts all data files from the metadata.
    Args:
        metadata (Dict): The metadata in JSON format.
    Returns:
        List[Dict]: A list of dictionaries containing data file information.
    """
    datafiles = []
    for d in metadata["@graph"]:
        if "name" in d.keys() and d["@type"] == "File":
            data_unit = {}
            data_unit["name"] = d["name"]
            try:
                # in MB
                data_unit["sizeMB"] = int(int(d["contentSize"]) / 1e6)
            except KeyError:
                data_unit["sizeMB"] = "unknown"

            try:
                data_unit["downloadUrl"] = d["downloadUrl"]
            except KeyError:
                data_unit["downloadUrl"] = "unknown"
            datafiles.append(data_unit)

    return datafiles
