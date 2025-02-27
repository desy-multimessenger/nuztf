import gzip
import io
import logging
from base64 import b64encode
from json import JSONDecodeError

import backoff
import numpy as np
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
def ampel_api_cutout(candid: int, logger=None):
    """Function to query ampel for cutouts by candidate ID"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if "v2" in API_ZTF_ARCHIVE_URL:
        queryurl_cutouts = API_ZTF_ARCHIVE_URL + f"/cutouts/{candid}"
    else:
        queryurl_cutouts = API_ZTF_ARCHIVE_URL + f"/alert/{candid}/cutouts"

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

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


def ensure_ampel_cutouts(alert: list, logger=None):
    """Make sure alert contains cutouts (if not, query them from AMPEL API)"""

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
