import os
import requests
from typing import Dict
from dotenv import load_dotenv

load_dotenv()


def get_ro_crate_metadata_gh(sample_id: str, mode: str = "private") -> Dict:
    """
    Retrieves RO-Crate metadata from a GitHub repository.

    Args:
        sample_id (str): The ID of the sample.
        mode (str, optional): The mode of access, either "private" or "public". Defaults to "private".

    Returns:
        Dict: The metadata in JSON format.
    """
    token = os.getenv("GH_TOKEN")

    metdat = f"https://api.github.com/repos/emo-bon/metaGOflow-rocrates-dvc/contents/{sample_id}-ro-crate/ro-crate-metadata.json"
    req = requests.get(
        metdat,
        headers={
            "accept": "application/vnd.github.v3.raw",
            "authorization": f"token {token}",
        },
    )
    print("ro-crate-metadata.json request status", req.status_code)
    return req.json()


def get_ro_crate_data(metadata_json: Dict, data_id: str, mode: str = "private"):
    """
    Retrieves RO-Crate data file based on metadata.

    Args:
        metadata_json (Dict): The metadata in JSON format.
        data_id (str): The ID of the data file.
        mode (str, optional): The mode of access, either "private" or "public". Defaults to "private".

    Raises:
        NotImplementedError: This function is not yet implemented.
    """
    raise NotImplementedError
