import glob
import logging
import os
import time
import warnings
from glob import glob
from pathlib import Path
from typing import Optional

import backoff
import numpy as np
import pandas as pd
import pyvo.dal
import requests
from astropy import units as u
from astropy.time import Time
from pyvo.auth import authsession, securitymethods
from ztfquery import skyvision
from ztfquery.fields import get_fields_containing_target
from ztfquery.io import LOCALSOURCE

from nuztf import credentials
from nuztf.fritz import fritz_api
from nuztf.observations_depot import (
    MNS,
    coverage_depot_path,
    coverage_dir,
    get_date,
    partial_flag,
)

logger = logging.getLogger(__name__)

username, password = credentials.load_credentials("irsa")

# Gather login information
data = {"username": username, "password": password}

headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

login_url = "https://irsa.ipac.caltech.edu"


def get_coverage(jds: [int]) -> Optional[pd.DataFrame]:
    # Clear any logs flagged as partial/incomplete

    cache_files = glob(f"{coverage_dir}/*.csv")
    partial_logs = [x for x in cache_files if partial_flag in x]

    if len(partial_logs) > 0:
        logger.debug(f"Removing the following partial logs: {partial_logs}")
        for partial_log in partial_logs:
            os.remove(partial_log)

    # Only write missing logs

    missing_logs = []

    for jd in jds:
        if not os.path.exists(coverage_path(jd)):
            missing_logs.append(jd)
        else:
            df = pd.read_csv(coverage_path(jd))
            if len(df) == 0:
                missing_logs.append(jd)

    if len(missing_logs) > 0:
        logger.debug(
            f"Some logs were missing from the cache. Querying for the following JDs: {missing_logs}"
        )
        write_coverage(missing_logs)

    # Load logs from cache

    results = []

    for jd in jds:
        res = pd.read_csv(coverage_path(jd))
        results.append(res)

    if results:
        result_df = pd.concat(results)
        return result_df
    else:
        return None


def get_obs_summary(t_min, t_max=None, max_days: int = None):
    """
    Get observation summary from Skyvision (or IRSA if Skyvision fails)
    """
    now = Time.now()

    if t_max and max_days:
        raise ValueError("Choose either t_max or max_days, not both")

    if t_max is None:
        if max_days is None:
            t_max = now
        else:
            t_max = t_min + (max_days * u.day)

    if t_max > now:
        t_max = now

    logger.info("Getting observation logs from skyvision.")
    mns = get_obs_summary_skyvision(t_min=t_min, t_max=t_max)

    # Useless for now as long as Fritz also depends on TAP

    # if len(mns.data) == 0:
    #     logger.debug("Empty observation log, try Fritz instead.")
    #     mns = get_obs_summary_fritz(t_min=t_min, t_max=t_max)

    if len(mns.data) == 0:
        logger.debug("Empty observation log, try IRSA instead.")
        mns = get_obs_summary_irsa(t_min=t_min, t_max=t_max)

    logger.debug(f"Found {len(mns.data)} observations in total.")

    return mns


def get_obs_summary_irsa(t_min, t_max):
    """
    Get observation summary from IRSA
    """
    jds = np.arange(int(t_min.jd), int(t_max.jd) + 1)

    logger.debug("Getting coverage")

    df = get_coverage(jds)
    mns = MNS(df)

    mns.data.query(f"obsjd >= {t_min.jd} and obsjd <= {t_max.jd}", inplace=True)
    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    logger.debug("Done")

    return mns


def get_obs_summary_skyvision(t_min, t_max):
    """
    Get observation summary from Skyvision
    """
    t_min_date = t_min.to_value("iso", subfmt="date")
    t_max_date = t_max.to_value("iso", subfmt="date")

    # ztfquery saves nightly observations in a cache, and does not redownload them.
    # If the nightly log was not complete, it will never be updated.
    # Here we simply clear the cache and cleanly re-download everything.

    logger.debug(
        f"Skyvision: Obtaining nightly observation logs from {t_min_date} to {t_max_date}"
    )

    skyvision_log = os.path.join(LOCALSOURCE, "skyvision")

    for filename in os.listdir(skyvision_log):
        if ".csv" in filename:
            path = os.path.join(skyvision_log, filename)
            os.remove(path)

    mns = skyvision.CompletedLog.from_daterange(
        t_min_date, end=t_max_date, verbose=False
    )

    mns.data["obsjd"] = Time(list(mns.data.datetime.values), format="isot").jd

    mns.data.query(f"obsjd >= {t_min.jd} and obsjd <= {t_max.jd}", inplace=True)

    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    logger.debug(f"Found {len(mns.data)} observations in total.")

    return mns


def get_obs_summary_fritz(t_min, t_max) -> MNS | None:
    """
    Get observation summary from Fritz
    """

    res = fritz_api(method="GET", endpoint_extension="api/observation", data=data)
    res = res.json().get("data", {}).get("observations")

    t_min_date = t_min.to_value("iso").split(".")[0]
    t_max_date = t_max.to_value("iso").split(".")[0]

    logger.debug(
        f"Fritz: Obtaining nightly observation logs from {t_min_date} to {t_max_date}"
    )

    params = {
        "telescopeName": "Palomar 1.2m Oschin",
        "startDate": t_min_date,
        "endDate": t_max_date,
    }

    res = fritz_api(method="GET", endpoint_extension="api/observation", data=params)
    res = res.json().get("data", {}).get("observations")

    if res is None:
        return res

    resdict = {}
    for i, entry in enumerate(res):
        resdict.update(
            {
                i: {
                    "datetime": entry["obstime"],
                    "date": Time(entry["obstime"]).to_value("iso", subfmt="date"),
                    "exptime": entry["exposure_time"],
                    "fid": entry["id"],
                    "field": entry["field"]["field_id"],
                    "obsjd": Time(entry["obstime"]).to_value("jd"),
                    "maglim": entry["limmag"],
                }
            }
        )

    df = pd.DataFrame.from_dict(resdict, orient="index")

    mns = MNS(df)

    logger.debug(f"Found {len(mns.data)} observations in total.")

    return mns
