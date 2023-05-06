#!/usr/bin/env python3

import gzip
import io
import logging
from base64 import b64encode
from json import JSONDecodeError

import backoff
import numpy as np
import requests
from ampel.util.json import load
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from astropy.io import fits  # type: ignore
from astropy.time import Time  # type: ignore
from requests.auth import HTTPBasicAuth

from nuztf import utils
from nuztf.credentials import load_credentials

API_BASEURL = "https://ampel.zeuthen.desy.de"
API_ZTF_ARCHIVE_URL = API_BASEURL + "/api/ztf/archive/v3"
API_CATALOGMATCH_URL = API_BASEURL + "/api/catalogmatch"

_, ampel_api_archive_token = load_credentials("ampel_api_archive_token")


def merge_alerts(alert_list: list) -> list:
    """ """
    merged_list = []
    keys = list(set([x["objectId"] for x in alert_list]))

    for objectid in keys:
        alerts = [x for x in alert_list if x["objectId"] == objectid]
        if len(alerts) == 1:
            merged_list.append(alerts[0])
        else:
            jds = [x["candidate"]["jd"] for x in alerts]
            order = [jds.index(x) for x in sorted(jds)[::-1]]
            latest = alerts[jds.index(max(jds))]
            latest["candidate"]["jdstarthist"] = min(
                [x["candidate"]["jdstarthist"] for x in alerts]
            )

            for index in order[1:]:
                x = alerts[index]

                # Merge previous detections

                for prv in x["prv_candidates"] + [x["candidate"]]:
                    if prv not in latest["prv_candidates"]:
                        latest["prv_candidates"] = [prv] + latest["prv_candidates"]

            merged_list.append(latest)

    return merged_list


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_cone(
    ra: float,
    dec: float,
    radius: float,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    logger=None,
) -> list:
    """Function to query ampel via a cone search"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    queryurl_conesearch = (
        API_ZTF_ARCHIVE_URL + f"/alerts/cone_search?ra={ra}&dec={dec}&"
        f"radius={radius}&jd_start={t_min_jd}&"
        f"jd_end={t_max_jd}&with_history={hist}&"
        f"with_cutouts={cutouts}&chunk_size={chunk_size}"
    )

    logger.debug(queryurl_conesearch)

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_conesearch,
        headers=headers,
    )

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_timerange(
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    logger=None,
) -> list:
    """Function to query ampel via a time-range search"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    queryurl_timerange = (
        API_ZTF_ARCHIVE_URL + f"/alerts/time_range?jd_start={t_min_jd}&"
        f"jd_end={t_max_jd}&with_history={hist}&"
        f"with_cutouts={cutouts}&chunk_size={chunk_size}"
    )

    logger.debug(queryurl_timerange)

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_timerange,
        headers=headers,
    )

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res


def ensure_cutouts(alert: list, logger=None):
    """Make sure alert contains cutouts (if not, query them from AMPEL API"""

    if logger is None:
        logger = logging.getLogger(__name__)

    candid = alert[0]["candid"]
    ztf_id = alert[0]["objectId"]

    if "cutoutScience" in alert[0].keys():
        if "stampData" in alert[0]["cutoutScience"].keys():
            logger.debug("Alert already contains cutouts.")

            return alert

    logger.debug(f"{ztf_id}: Querying API for cutouts.")

    final_cutouts = {}

    cutouts = ampel_api_cutout(candid)

    if "detail" in cutouts.keys():
        if cutouts["detail"] == "Not Found":
            for k in [
                "science",
                "difference",
                "template",
                "Cutoutscience",
                "Cutoutdifference",
                "Cutouttemplate",
            ]:
                final_cutouts[f"cutout{k.title()}"] = {
                    "stampData": create_empty_cutout()
                }
    else:
        for k in cutouts:
            final_cutouts[f"cutout{k.title()}"] = {"stampData": cutouts[k]}

    alert[0] = {**alert[0], **final_cutouts}

    logger.debug(f"{ztf_id}: Added cutouts.")

    return alert


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_name(
    ztf_name: str,
    with_history: bool = True,
    with_cutouts: bool = False,
    limit: int = 999999,
    logger=None,
) -> list:
    """Function to query ampel via name"""
    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    queryurl_ztf_name = (
        API_ZTF_ARCHIVE_URL
        + f"/object/{ztf_name}/alerts?with_history={hist}&with_cutouts={cutouts}&limit={limit}"
    )

    logger.debug(queryurl_ztf_name)

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_ztf_name,
        headers=headers,
    )

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        query_res = [i for i in response.json()]
        query_res = merge_alerts(query_res)

    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_lightcurve(
    ztf_name: str,
    t_min_jd=Time("2017-01-01T00:00:00.0", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    program_id: int = None,
    logger=None,
) -> list:
    """
    Function to query ampel via name, returns a virtual alert
    constructed by AMPEL containing ALL photopoints and upper limits

    """

    if logger is None:
        logger = logging.getLogger(__name__)

    if program_id is None:
        queryurl_lightcurve = (
            API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/photopoints?jd_start={t_min_jd}&"
            f"jd_end={t_max_jd}"
        )
    else:
        queryurl_lightcurve = (
            API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/photopoints?jd_start={t_min_jd}&"
            f"jd_end={t_max_jd}&programid={program_id}"
        )

    logger.debug(queryurl_lightcurve)

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_lightcurve,
        headers=headers,
    )

    if response.status_code == 503:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    try:
        query_res = [response.json()]

    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_healpix(
    ipix: int,
    nside: int = 64,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    logger=None,
) -> list:
    """Function to query ampel based on a healpix pixel-index (nside is the pixelization degree)"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    queryurl_healpix = (
        API_ZTF_ARCHIVE_URL
        + f"/alerts/healpix?nside={nside}&ipix={ipix}&jd_start={t_min_jd}&jd_end={t_max_jd}&with_history={hist}&with_cutouts={cutouts}&chunk_size={chunk_size}"
    )

    logger.debug(queryurl_healpix)

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_healpix,
        headers=headers,
    )

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=1200,
)
def ampel_api_acknowledge_chunk(resume_token: str, chunk_id: int, logger=None):
    """
    After receiving a chunk, acknowledge that we got it
    (otherwise large alert queries will start looping)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    endpoint = (
        API_ZTF_ARCHIVE_URL + f"/stream/{resume_token}/chunk/{chunk_id}/acknowledge"
    )

    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {ampel_api_archive_token}",
    }

    payload = {"resume_token": resume_token, "chunk_id": chunk_id}

    logger.debug(f"Acknowledging:\n{payload}")

    response = requests.post(url=endpoint, json=payload, headers=headers)


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=1200,
)
def ampel_api_skymap(
    pixels: list,
    nside: int = 64,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    resume_token: str = None,
    warn_exceeding_chunk: bool = True,
    program_id: int = None,
    logger=None,
) -> tuple:
    """
    Function to query ampel based on a list of healpix pixels (nside is the pixelization degree)
    """

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    # First, we create a json body to post
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {ampel_api_archive_token}",
    }
    if "v2" in API_ZTF_ARCHIVE_URL:
        lt = "lt"
        gt = "gt"
    else:
        lt = "$lt"
        gt = "$gt"

    # Now we reduce the query size
    regions = utils.deres(nside=nside, ipix=pixels)

    query = {
        "regions": regions,
        "jd": {
            lt: t_max_jd,
            gt: t_min_jd,
        },
        "latest": "false",
        "with_history": hist,
        "with_cutouts": cutouts,
        "chunk_size": chunk_size,
    }

    if resume_token:
        query["resume_token"] = resume_token

    if program_id is not None:
        query["programid"] = program_id

    queryurl_skymap = API_ZTF_ARCHIVE_URL + f"/alerts/healpix/skymap"

    logger.debug(f"Query url:\n{queryurl_skymap}")
    logger.debug(f"Query:\n{query}")

    response = requests.post(url=queryurl_skymap, json=query, headers=headers)

    logger.debug(response)
    logger.debug(response.status_code)

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        res_json = response.json()
        remaining_chunks = res_json["remaining"]["chunks"]
        logger.debug(f"Remaining chunks: {remaining_chunks}")
        chunk_id = res_json.get("chunk", None)
        resume_token = response.json().get("resume_token", None)
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size and warn_exceeding_chunk:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res, resume_token, chunk_id, remaining_chunks


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_cutout(candid: int, logger=None):
    """Function to query ampel for cutouts by candidate ID"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if "v2" in API_ZTF_ARCHIVE_URL:
        queryurl_cutouts = API_ZTF_ARCHIVE_URL + f"/cutouts/{candid}"
    else:
        queryurl_cutouts = API_ZTF_ARCHIVE_URL + f"/alert/{candid}/cutouts"

    headers = {"Authorization": f"Bearer {ampel_api_archive_token}"}

    response = requests.get(
        queryurl_cutouts,
        headers=headers,
    )

    logger.debug(queryurl_cutouts)

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        cutouts = response.json()
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    return cutouts


def create_empty_cutout():
    """Function to reate an empty image for missing cutouts"""
    npix = 63

    blank = np.ones((npix, npix))

    for i in range(npix):
        c = abs(npix / 2 - i) / (0.5 * npix)
        blank[i - 1][i - 1] = c
        blank[i - 1][npix - i - 1] = c

    hdu = fits.PrimaryHDU(blank)
    hdul = fits.HDUList([hdu])
    comp = io.BytesIO()
    hdul.writeto(comp)
    blank_compressed = gzip.compress(comp.getvalue())
    blank_compressed = b64encode(blank_compressed)

    return blank_compressed


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_catalog(
    catalog: str,
    catalog_type: str,
    ra_deg: float,
    dec_deg: float,
    search_radius_arcsec: float = 10,
    search_type: str = "all",
    logger=None,
):
    """
    Method for querying catalogs via the Ampel API
    'catalog' must be the name of a supported catalog, e.g.
    SDSS_spec, PS1, NEDz_extcats...
    For a full list of catalogs, confer
    https://ampel.zeuthen.desy.de/api/catalogmatch/catalogs

    """
    assert catalog_type in ["extcats", "catsHTM"]
    assert search_type in ["all", "nearest"]

    if logger is None:
        logger = logging.getLogger(__name__)

    queryurl_catalogmatch = API_CATALOGMATCH_URL + "/cone_search/" + search_type

    # First, we create a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    query = {
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "catalogs": [
            {"name": catalog, "rs_arcsec": search_radius_arcsec, "use": catalog_type}
        ],
    }

    logger.debug(queryurl_catalogmatch)
    logger.debug(query)

    response = requests.post(url=queryurl_catalogmatch, json=query, headers=headers)

    if response.status_code == 503:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    try:
        res = response.json()[0]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    return res


def get_preprocessed_results(file_basename: str, logger=None) -> None | list:
    """
    Access the DESY Cloud to look if there are precomputed results from an AMPEL run there
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    desy_cloud_token = load_credentials("desy_cloud_token", token_based=True)

    filename = file_basename + ".json.gz"

    res = requests.get(
        f"https://syncandshare.desy.de/public.php/webdav/{filename}",
        headers={"X-Requested-With": "XMLHttpRequest"},
        auth=(desy_cloud_token, "bla"),
    )

    if res.status_code != 200:
        logger.warning(
            "\n\n-------------------- !! -------------------\nSomething went wrong with your query.\nCheck your credentials and make sure Ampel\nhas run correctly at Desy.\n-------------------- !! -------------------\n\n"
        )
        return None

    with open(f"{filename}", "wb") as f:
        f.write(res.content)

    res = []
    with gzip.open(filename, "rb") as f_in:
        data = load(f_in)
        for t in data:
            ztf_id = ZTFIdMapper.to_ext_id(t.stock.get("stock"))
            pp = t.get_photopoints()
            pp_reformatted = utils.reformat_downloaded_results(
                photopoints=pp, ztf_id=ztf_id
            )
            res.append(pp_reformatted)

    return res
