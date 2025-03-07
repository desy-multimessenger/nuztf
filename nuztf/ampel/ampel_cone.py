import logging
from json import JSONDecodeError

import backoff
import requests
from astropy.io import fits  # type: ignore
from astropy.time import Time  # type: ignore

from nuztf.ampel.urls import API_ZTF_ARCHIVE_URL
from nuztf.ampel.utils import get_ampel_token


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
    """Query ampel via a cone search"""

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

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

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
