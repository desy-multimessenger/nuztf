import io, logging, gzip
import requests
import backoff
from base64 import b64decode
from json import JSONDecodeError
import numpy as np

from astropy.time import Time
from astropy.io import fits
from requests.auth import HTTPBasicAuth

from nuztf.credentials import load_credentials

# AMPEL API URLs

API_BASEURL = "https://ampel.zeuthen.desy.de"
API_ZTF_ARCHIVE_URL = API_BASEURL + "/api/ztf/archive/v2"
API_CATALOGMATCH_URL = API_BASEURL + "/api/catalogmatch"
API_CUTOUT_URL = API_BASEURL + "/api/ztf/archive/v2/cutouts"

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

    return query_res


def add_cutouts(alert: list):
    candid = alert[0]["candid"]
    cutouts = ampel_api_cutout(candid)

    final_cutouts = {}

    if "detail" in cutouts.keys():
        if cutouts["detail"] == "Not Found":
            for k in ["science", "difference", "template"]:
                final_cutouts[f"cutout{k.title()}"] = {"data": create_empty_cutout()}
    else:
        for k in cutouts:
            final_cutouts[f"cutout{k.title()}"] = {
                "data": b64decode(cutouts[k]),
            }

    alert[0].update({"cutouts": final_cutouts})


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_name(
    ztf_name: str, with_history: bool = True, with_cutouts: bool = False, logger=None
) -> list:
    """Function to query ampel via name"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    queryurl_ztf_name = (
        API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/alerts?with_history={hist}"
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

    if with_cutouts:
        add_cutouts(query_res)

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_healpix(
    ipix: int,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    logger=None,
) -> list:
    """Function to query ampel based on a healpix pixel-index (gird has nside=64)"""

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
        + f"/alerts/healpix?ipix={ipix}&jd_start={t_min_jd}&jd_end={t_max_jd}&with_history={hist}&with_cutouts={cutouts}&chunk_size={chunk_size}"
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

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_cutout(candid: int, logger=None):
    """Function to query ampel for cutouts by candidate ID"""

    if logger is None:
        logger = logging.getLogger(__name__)

    queryurl_cutouts = API_CUTOUT_URL + f"/{candid}"

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
    searchradius_arcsec: float = 10,
    searchtype: str = "all",
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
    assert searchtype in ["all", "nearest"]

    if logger is None:
        logger = logging.getLogger(__name__)

    queryurl_catalogmatch = API_CATALOGMATCH_URL + "/cone_search/" + searchtype

    # First, we create a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    query = {
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "catalogs": [
            {"name": catalog, "rs_arcsec": searchradius_arcsec, "use": catalog_type}
        ],
    }

    logger.debug(queryurl_catalogmatch)
    logger.debug(query)

    response = requests.post(url=queryurl_catalogmatch, json=query, headers=headers)

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        res = response.json()[0]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    return res
