#!/usr/bin/env python3
# coding: utf-8

import json
import logging
import math
import re
from collections import OrderedDict, defaultdict
from json import JSONDecodeError

import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from nuztf.credentials import load_credentials
from requests.auth import HTTPBasicAuth

# same cosmology everywhere
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def is_icecube_name(name) -> bool:
    """
    Checks if a string adheres to the IceCube naming scheme
    (e.g. IC201021B)
    """
    if re.match(
        r"^IC((\d{2}((0[13578]|1[02])(0[1-9]|[12]\d|3[01])|(0[13456789]|1[012])(0[1-9]|[12]\d|30)|02(0[1-9]|1\d|2[0-8])))|([02468][048]|[13579][26])0229)[a-zA-Z]$",
        name,
    ):
        match = True
    else:
        match = False
    return match


def is_ligo_name(name) -> bool:
    """
    Checks if a string adheres to the LVT naming scheme
    (e.g. S190814bv)
    """
    if re.match(
        r"^(S|GW|MS)((\d{2}((0[13578]|1[02])(0[1-9]|[12]\d|3[01])|(0[13456789]|1[012])(0[1-9]|[12]\d|30)|02(0[1-9]|1\d|2[0-8])))|([02468][048]|[13579][26])0229)[a-z]{0,4}$",
        name,
    ):
        match = True
    else:
        match = False
    return match


def is_ztf_name(name: str) -> bool:
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    if re.match(r"^ZTF[1-2]\d[a-z]{7}$", name):
        match = True
    else:
        match = False
    return match


def is_tns_name(name: str) -> bool:
    """
    Checks if a string adheres to the TNS naming scheme
    """
    if re.match(r"^AT|SN(19|20)\d\d[a-z]{3,4}$", name):
        matches = True
    else:
        matches = False
    return matches


def reformat_downloaded_results(photopoints: list, ztf_id: str) -> dict:
    """
    Massage the TransientView photopoint output so it matches what the archive DB returns
    """

    # Find the index of the latest detection
    latest_jd_index = np.argmax(np.asarray([pp["body"]["jd"] for pp in photopoints]))

    resdict = {
        "candid": photopoints[latest_jd_index]["id"],
        "objectId": ztf_id,
        "schemavsn": "3.3",
        "publisher": "Ampel",
        "candidate": photopoints[latest_jd_index]["body"],
    }

    prv_candidates = []

    for i, pp in enumerate(photopoints):
        if i != latest_jd_index:  # dict is readonly, we can't pop
            prv_candidates.append(pp["body"])

    resdict["prv_candidates"] = prv_candidates

    return resdict


def deres(nside, ipix, min_nside=1):
    """
    Originally from Ampel-ZTF-archive/ampel/ztf/archive/server (by JvS)

    Decompose a set of (nested) HEALpix indices into sets of complete superpixels at lower resolutions.
    :param nside: nside of given indices
    :param ipix: pixel indices
    :min_nside: minimum nside of complete pixels
    """
    remaining_pixels = set(ipix)
    decomposed = defaultdict(list)
    for log2_nside in range(int(math.log2(min_nside)), int(math.log2(nside)) + 1):
        super_nside = 2**log2_nside
        # number of base_nside pixels per nside superpixel
        scale = (nside // super_nside) ** 2
        # sort remaining base_nside pixels by superpixel
        by_superpixel = defaultdict(list)
        for pix in remaining_pixels:
            by_superpixel[pix // scale].append(pix)
        # represent sets of pixels that fill a superpixel
        # as a single superpixel, and remove from the working set
        for superpix, members in by_superpixel.items():
            if len(members) == scale:
                decomposed[super_nside].append(superpix)
                remaining_pixels.difference_update(members)

    decomposed_dict = dict(decomposed)

    healpix_regions = [
        {"nside": nside, "pixels": pixels} for nside, pixels in decomposed_dict.items()
    ]

    return healpix_regions


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
