import json
import logging
import os
from glob import glob

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.time import Time
from requests.auth import HTTPBasicAuth
from requests.exceptions import HTTPError
from tqdm import tqdm

from nuztf import credentials
from nuztf.observations import MNS, coverage_dir, partial_flag

logger = logging.getLogger(__name__)

username, password = credentials.load_credentials("ipacdepot")


class NoDepotEntry(Exception):
    """
    No entry in the depot for a given date
    """


def download_depot_log(date):
    """
    Download the depot log for a given date

    :param date: date in YYYYMMDD format
    :return: json log
    """
    url = f"https://ztfweb.ipac.caltech.edu/ztf/depot/{date}/ztf_recentproc_{date}.json"
    response = requests.get(url, auth=HTTPBasicAuth(username, password))
    if response.status_code == 404:
        raise NoDepotEntry(f"No depot entry for {date}")
    response.raise_for_status()
    return response.json()


def coverage_depot_cache(jd: float) -> str:
    """
    Return the path to the cached coverage file for a given JD

    :param jd: JD
    :return: path to cached coverage file
    """
    if (Time.now().jd - jd) < 1:
        partial_ext = partial_flag
    else:
        partial_ext = ""

    return os.path.join(coverage_dir, f"{jd}{partial_ext}.json")


def write_depot_coverage(jds: list[int]):
    """
    Write the depot coverage for a list of JDs to the cache

    :param jds: JDs
    :return: None
    """
    for jd in tqdm(jds):
        try:
            date = str(Time(jd, format="jd").isot).split("T")[0].replace("-", "")
            log = download_depot_log(date)
        except NoDepotEntry:
            log = {}
        path = coverage_depot_cache(jd)
        with open(path, "w") as f:
            json.dump(log, f)


def get_coverage_depot(jds: [int]) -> pd.DataFrame | None:
    """
    Get a dataframe of the depot coverage for a list of JDs

    :param jds: JDs
    :return: Coverage dataframe
    """
    # Clear any logs flagged as partial/incomplete

    cache_files = glob(f"{coverage_dir}/*.json")
    partial_logs = [x for x in cache_files if partial_flag in x]

    if len(partial_logs) > 0:
        logger.debug(f"Removing the following partial logs: {partial_logs}")
        for partial_log in partial_logs:
            os.remove(partial_log)

    # Only write missing logs

    missing_logs = []

    for jd in jds:
        if not os.path.exists(coverage_depot_cache(jd)):
            missing_logs.append(jd)
        else:
            df = pd.read_json(coverage_depot_cache(jd))
            if len(df) == 0:
                missing_logs.append(jd)

    if len(missing_logs) > 0:
        logger.debug(
            f"Some logs were missing from the cache. "
            f"Querying for the following JDs: {missing_logs}"
        )
        write_depot_coverage(missing_logs)

    # Load logs from cache

    results = []

    for jd in tqdm(jds):
        res = pd.read_json(coverage_depot_cache(jd))
        results.append(res)

    if results:
        result_df = pd.concat(results, ignore_index=True)
        return result_df
    else:
        return None


def get_obs_summary_depot(t_min: Time, t_max: Time) -> MNS | None:
    """
    Get observation summary from depot
    """

    jds = np.arange(t_min.jd, t_max.jd + 1)

    res = get_coverage_depot(jds)

    if len(res) == 0:
        return None

    res["date"] = Time(res["obsjd"].to_numpy(), format="jd").isot

    mns = MNS(df=res)

    mns.data.query(f"obsjd >= {t_min.jd} and obsjd <= {t_max.jd}", inplace=True)

    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    return mns


def get_obs_summary(t_min, t_max=None, max_days: float = None) -> MNS | None:
    """
    Get observation summary from IPAC depot
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

    logger.info("Getting observation logs from IPAC depot.")
    mns = get_obs_summary_depot(t_min=t_min, t_max=t_max)

    if mns is not None:
        logger.debug(
            f"Found {len(set(mns.data['exposure_id']))} observations in total."
        )
    else:
        logger.debug("Found no observations on IPAC depot.")

    return mns
