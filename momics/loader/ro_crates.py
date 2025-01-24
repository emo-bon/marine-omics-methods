import os
import requests
from dotenv import load_dotenv
load_dotenv()


def get_ro_crate_metadata_gh(sample_id, mode='private'):
    # this return json
    token = os.getenv('GH_TOKEN')

    metdat = f"https://api.github.com/repos/emo-bon/metaGOflow-rocrates-dvc/contents/{sample_id}-ro-crate/ro-crate-metadata.json"
    req = requests.get(metdat, 
                    headers={
                        'accept': 'application/vnd.github.v3.raw',
                        'authorization': f'token {token}',
                        })
    print('ro-crate-metadata.json request status', req.status_code)
    return req.json()


def get_ro_crate_data(metadata_json, data_id, mode='private'):
    # this return data file
    return






