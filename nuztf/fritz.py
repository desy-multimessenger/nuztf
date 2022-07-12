import requests
import backoff
import os
from nuztf.credentials import load_credentials

# Fritz API URLs

API_BASEURL = "https://fritz.science"

fritz_token = load_credentials("fritz", token_based=True)


def fritz_api(method: str, endpoint_extension: str, data: dict = None):
    headers = {"Authorization": f"token {fritz_token}"}
    endpoint = os.path.join(API_BASEURL, endpoint_extension)
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=60,
)
def save_source_to_group(object_id: str, group_id: int):
    payload = {
        "objId": object_id,
        "inviteGroupIds": [group_id],
    }
    return fritz_api(
        method="POST", endpoint_extension="api/source_groups", data=payload
    )


def delete_source_from_group(object_id: str, group_id: int):
    payload = {
        "objId": object_id,
        "unsaveGroupIds": [group_id],
    }
    return fritz_api(
        method="POST", endpoint_extension="api/source_groups", data=payload
    )
