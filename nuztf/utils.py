#!/usr/bin/env python3
# coding: utf-8

import re
import logging
import requests
import json
from json import JSONDecodeError
from collections import OrderedDict

from astropy.cosmology import FlatLambdaCDM
from nuztf.credentials import load_credentials
from requests.auth import HTTPBasicAuth

# same cosmology everywhere
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# test


def is_ztf_name(name: str):
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


def is_tns_name(name: str):
    """
    Checks if a string adheres to the TNS naming scheme
    """
    return re.match("^AT|SN(19|20)\d\d[a-z]{3,4}$", name)


def query_tns_by_name(name, logger=None):
    """
    Query the TNS API for a given name
    """
    if not is_tns_name(name):
        raise ValueError("String is not a TNS name")

    tns_api_token = load_credentials("tns_api_token", token_based=True)

    name = name[2:]

    tns_bot_id = "115364"
    tns_bot_name = "ZTF_DESY"

    if logger is None:
        logger = logging.getLogger(__name__)

    # logger = logging.getLogger(__name__)
    logger.info(f"Obtaining TNS information for {name}.")

    queryurl_tns = "https://www.wis-tns.org/api/get/object"

    tns_marker = (
        'tns_marker{"tns_id": "'
        + str(tns_bot_id)
        + '", "type": "bot", "name": "'
        + tns_bot_name
        + '"}'
    )
    headers = {"User-Agent": tns_marker}

    get_obj = [("objname", name), ("objid", ""), ("photometry", "1"), ("spectra", "0")]
    json_file = OrderedDict(get_obj)

    get_data = {"api_key": tns_api_token, "data": json.dumps(json_file)}

    response = requests.post(queryurl_tns, headers=headers, data=get_data)

    try:
        res = response.json()
        if "name" in res["data"]["reply"].keys():
            if "110" in res["data"]["reply"]["name"].keys():
                return None
        else:
            return res

    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException
