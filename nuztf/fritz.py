import requests
import os
from nuztf.credentials import load_credentials

# Fritz API URLs

API_BASEURL = "https://fritz.science"

fritz_token = load_credentials("fritz")


def fritz_api(method, endpoint_extension, data=None):
    headers = {"Authorization": f"token {fritz_token}"}
    endpoint = os.path.join(API_BASEURL, endpoint_extension)
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response
