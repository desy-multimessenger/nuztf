import logging
import time
from pathlib import Path

import backoff
import numpy as np
import pyvo.dal
import requests
from astropy.time import Time
from pyvo.auth import authsession, securitymethods

from nuztf.observations.depot import data
from nuztf.observations.shared import coverage_dir, get_date, partial_flag

logger = logging.getLogger(__name__)


# IPAC TAP login
headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "text/plain"}
IPAC_TAP_URL = "https://irsa.ipac.caltech.edu"

# The following functions are used to find the paths to the cached coverage files


def coverage_tap_path(jd: float) -> Path:
    """
    Find the path to the TAP coverage file for a given JD

    :param jd: JD
    :return: Output path
    """
    if (Time.now().jd - jd) < 1:
        partial_ext = partial_flag
    else:
        partial_ext = ""

    output_path = coverage_dir.joinpath(f"{get_date(jd)}{partial_ext}_tap.json")

    return output_path


# The following functions are used to write the coverage to the cache


@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=60)
def write_coverage_tap(jds: [float]):
    """
    Write the coverage for a list of JDs to the cache using TAP

    (JDs must be half-integer)

    :param jds: Half-integer JDs
    :return: None
    """
    with requests.Session() as session:
        # Create a session and do the login.
        # The cookie will end up in the cookie jar of the session.
        response = session.post(IPAC_TAP_URL, data=data, headers=headers)
        response.raise_for_status()
        auth = authsession.AuthSession()
        auth.credentials.set(securitymethods.ANONYMOUS, session)
        client = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP", auth)

        for jd in jds:
            assert jd - int(jd) == 0.5, "JD must be a half-integer"

            obstable = client.search(
                f"""
            SELECT expid,obsjd,fid,field,exptime,rcid,maglimit,infobits,ipac_gid
            FROM ztf.ztf_current_meta_sci WHERE (obsjd BETWEEN {jd} AND {jd + 1})
            """
            ).to_table()
            names = (
                "expid",
                "exptime",
                "field",
                "fid",
                "rcid",
                "maglimit",
                "infobits",
                "ipac_gid",
            )
            renames = (
                "exposure_id",
                "exposure_time",
                "field_id",
                "filter_id",
                "qid",
                "maglim",
                "status",
                "programid",
            )
            obstable.rename_columns(names, renames)

            obs = obstable.to_pandas()
            obs["nalertpackets"] = np.nan

            output_path = coverage_tap_path(jd)

            obs.to_json(output_path)
            time.sleep(1)

            logger.warning(
                f"Coverage for JD {jd} written to {output_path} using TAP, "
                f"but this will only include programid=1"
            )
