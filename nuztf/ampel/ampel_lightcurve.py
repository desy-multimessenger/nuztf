import logging
from json import JSONDecodeError

import backoff
import requests
from astropy.io import fits  # type: ignore
from astropy.time import Time  # type: ignore

from nuztf.ampel.urls import API_ZTF_ARCHIVE_URL
from nuztf.ampel.utils import get_ampel_token, merge_alerts


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

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

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

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

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
def ampel_api_alerts(
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
            API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/alerts?jd_start={t_min_jd}&"
            f"jd_end={t_max_jd}&with_history=true"
        )
    else:
        queryurl_lightcurve = (
            API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/alerts?jd_start={t_min_jd}&"
            f"jd_end={t_max_jd}&programid={program_id}"
        )

    logger.debug(queryurl_lightcurve)

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

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
