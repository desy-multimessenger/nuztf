import logging
from json import JSONDecodeError

import backoff
import requests
from astropy.io import fits  # type: ignore
from astropy.time import Time  # type: ignore

from nuztf.ampel.urls import API_CATALOGMATCH_URL


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
    if not catalog_type in ["extcats", "catsHTM"]:
        raise ValueError(
            f"Expected parameter catalog_type in ['extcats', 'catsHTM'], got {catalog_type}"
        )
    if not search_type in ["all", "nearest"]:
        raise ValueError(
            f"Expected parameter catalog_type in ['all', 'nearest'], got {search_type}"
        )

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
