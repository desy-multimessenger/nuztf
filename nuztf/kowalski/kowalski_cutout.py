"""
This module contains functions to query Kowalski for cutouts
"""

import logging
from base64 import b64encode

from penquins import Kowalski

from nuztf.kowalski.config import get_kowalski

logger = logging.getLogger(__name__)


def kowalski_api_cutout(
    candid: int,
    kowalski: Kowalski | None = None,
):
    """
    Download alert data from Kowalski

    :param candid: Candidate ID
    :param kowalski: Kowalski object
    :return: Alert data
    """
    if kowalski is None:
        kowalski = get_kowalski()

    query_config = {
        "query_type": "find",
        "query": {
            "catalog": "ZTF_alerts",
            "filter": {
                "candid": {"$eq": int(candid)},
            },
            "projection": {
                "cutoutScience": 1,
                "cutoutTemplate": 1,
                "cutoutDifference": 1,
            },
        },
    }

    query_result = kowalski.query(query_config)

    if "data" in query_result:
        alerts = query_result["data"]
    else:
        alerts = query_result.get("default").get("data")

    cutouts = alerts[-1]

    return cutouts


def ensure_kowalski_cutouts(alert: list[dict]) -> list[dict]:
    """
    Make sure alert contains cutouts (if not, query them from Kowalski API)

    :param alert: Alert data
    :return: Alert data with cutouts
    """
    candid = alert[0]["candid"]
    ztf_id = alert[0]["objectId"]

    if "cutoutScience" in alert[0].keys():
        if "stampData" in alert[0]["cutoutScience"].keys():
            logger.debug("Alert already contains cutouts.")
            return alert

    logger.debug(f"{ztf_id}: Querying API for cutouts.")

    cutouts = kowalski_api_cutout(candid)

    for k in ["Science", "Difference", "Template"]:
        key = f"cutout{k}"
        cutouts[key]["stampData"] = b64encode((cutouts[key]["stampData"]))
        alert[0][key] = cutouts[key]

    logger.debug(f"{ztf_id}: Added cutouts.")

    return alert
