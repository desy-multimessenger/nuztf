import json
import logging
import os
import time
from glob import glob
from pathlib import Path

import backoff
import numpy as np
import pandas as pd
import pyvo.dal
import requests
from astropy import units as u
from astropy.time import Time
from pyvo.auth import authsession, securitymethods
from requests.auth import HTTPBasicAuth
from requests.exceptions import HTTPError
from tqdm import tqdm
from ztfquery import skyvision
from ztfquery.io import LOCALSOURCE

from nuztf import credentials

coverage_dir = Path(LOCALSOURCE).joinpath("all_obs")
coverage_dir.mkdir(exist_ok=True)

partial_flag = "_PARTIAL"

logger = logging.getLogger(__name__)

# IPAC TAP login
username, password = credentials.load_credentials("ipacdepot")
data = {"username": username, "password": password}
headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "text/plain"}
IPAC_TAP_URL = "https://irsa.ipac.caltech.edu"


class MNS:
    def __init__(self, df):
        self.data = df


class NoDepotEntry(Exception):
    """
    No entry in the depot for a given date
    """


def get_date(jd: float) -> str:
    """
    Get the date in YYYYMMDD format from a JD

    :param jd: JD
    :return: String date
    """
    return str(Time(jd, format="jd").isot).split("T")[0].replace("-", "")


def download_depot_log(date):
    """
    Download the depot log for a given date

    :param date: date in YYYYMMDD format
    :return: json log
    """
    url = f"https://ztfweb.ipac.caltech.edu/ztf/depot/{date}/ztf_recentproc_{date}.json"
    response = requests.get(url, auth=HTTPBasicAuth(username, password))
    if response.status_code == 404:
        raise NoDepotEntry(f"No depot entry for {date} at url {url}")
    response.raise_for_status()
    return response.json()


# The following functions are used to find the paths to the cached coverage files


def coverage_depot_path(jd: float) -> Path:
    """
    Return the path to the cached coverage file for a given JD

    :param jd: JD
    :return: path to cached coverage file
    """

    date = get_date(jd)

    if (Time.now().jd - jd) < 1:
        partial_ext = partial_flag
    else:
        partial_ext = ""

    return coverage_dir.joinpath(f"{date}{partial_ext}.json")


def coverage_skyvision_path(jd: float) -> Path:
    """
    Find the path to the Skyvision coverage file for a given JD

    :param jd: JD
    :return: Output path
    """
    if (Time.now().jd - jd) < 1:
        partial_ext = partial_flag
    else:
        partial_ext = ""

    output_path = coverage_dir.joinpath(f"{get_date(jd)}{partial_ext}_skyvision.json")

    return output_path


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


def write_coverage_depot(jds: list[float]):
    """
    Write the depot coverage for a list of JDs to the cache

    Requires a valid IPAC depot login

    :param jds: JDs
    :return: None
    """
    for jd in tqdm(jds):
        date = get_date(jd)
        try:
            log = download_depot_log(date)
        except NoDepotEntry as exc:
            log = {}
            logger.warning(f"No depot entry for {date}: {exc}")

        path = coverage_depot_path(jd)
        with open(path, "w") as f:
            json.dump(log, f)


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


def write_coverage_skyvision(jds: list[float]):
    """
    Write the Skyvision coverage for a list of JDs to the cache

    :param jds: JDs
    :return: None
    """
    for jd in tqdm(jds):
        date = get_date(jd)

        assert jd - int(jd) == 0.5, "JD must be a half-integer"

        res = skyvision.get_log(f"{date[:4]}-{date[4:6]}-{date[6:8]}", verbose=False)

        path = coverage_skyvision_path(jd)

        if res is not None:
            # The skyvision log has a bug where some entries have a FieldID of "NONE"
            mask = (
                pd.notnull(res["FieldID"]) & res["FieldID"].astype(str).str.isnumeric()
            )
            res = res[mask]
            jds = [
                Time(f'{row["UT Date"]}T{row["UT Time"]}').jd
                for _, row in res.iterrows()
            ]
            res["obsjd"] = jds
            res["status"] = (res["Observation Status"] == "FAILED").astype(int)
            res["filter_id"] = res["Filter"].apply(
                lambda x: 1 if x == "FILTER_ZTF_G" else 2 if x == "FILTER_ZTF_R" else 3
            )
            res["maglim"] = 20.5
            res["field_id"] = res["FieldID"].astype(int)
            res["exposure_time"] = res["Exptime"]
            res = res[
                ["obsjd", "filter_id", "field_id", "exposure_time", "maglim", "status"]
            ]

            new_res = []
            for _, row in res.iterrows():
                new = row.to_dict()
                new_res += [dict(qid=int(i), **new) for i in range(64)]
            pd.DataFrame(new_res).to_json(path)

        else:
            with open(path, "w") as f:
                json.dump({}, f)


# The following functions are used to get the coverage from the cache


def get_coverage(jds: [int]) -> pd.DataFrame | None:
    """
    Get a dataframe of the coverage for a list of JDs

    Will use the cache if available, otherwise will query the depot, and lastly TAP

    :param jds: JDs
    :return: Coverage dataframe
    """
    # Clear any logs flagged as partial/incomplete

    cache_files = coverage_dir.glob("*.json")
    partial_logs = [x for x in cache_files if partial_flag in str(x)]

    if len(partial_logs) > 0:
        logger.debug(f"Removing the following partial logs: {partial_logs}")
        for partial_log in partial_logs:
            partial_log.unlink()

    # Only write missing logs

    missing_logs = []

    for jd in jds:
        depot_path = coverage_depot_path(jd)
        if not depot_path.exists():
            missing_logs.append(jd)
        else:
            df = pd.read_json(coverage_depot_path(jd))
            if len(df) == 0:
                missing_logs.append(jd)

    if len(missing_logs) > 0:
        logger.info(
            f"Some logs were missing from the cache. "
            f"Querying for the following JDs in depot: {missing_logs}"
        )
        write_coverage_depot(missing_logs)

    # Try skyvision for missing logs
    still_missing_logs = []
    for jd in missing_logs:
        skyvision_path = coverage_skyvision_path(jd)
        if not skyvision_path.exists():
            still_missing_logs.append(jd)
        else:
            df = pd.read_json(coverage_skyvision_path(jd))
            if len(df) == 0:
                still_missing_logs.append(jd)

    if len(still_missing_logs) > 0:
        logger.info(
            f"Some logs were still missing from the cache. "
            f"Querying for the following JDs from skyvision: {still_missing_logs}"
        )
        write_coverage_skyvision(still_missing_logs)

    # Try TAP for missing logs

    completely_missing_logs = []

    for jd in still_missing_logs:
        depot_path = coverage_depot_path(jd)
        if not depot_path.exists():
            completely_missing_logs.append(jd)
        else:
            df = pd.read_json(coverage_depot_path(jd))
            if len(df) == 0:
                completely_missing_logs.append(jd)

    if len(completely_missing_logs) > 0:
        logger.info(
            f"Some logs were still missing from the cache. "
            f"Querying for the following JDs from TAP: {completely_missing_logs}"
        )
        write_coverage_tap(completely_missing_logs)

    # Load logs from cache

    results = []

    for jd in tqdm(jds):
        res = pd.read_json(coverage_depot_path(jd))
        if len(res) > 0:
            results.append(res)
        else:
            res = pd.read_json(coverage_skyvision_path(jd))
            if len(res) > 0:
                results.append(res)
            else:
                res = pd.read_json(coverage_tap_path(jd))
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

    jds = np.arange(int(t_min.jd) - 0.5, int(t_max.jd) + 1.5)

    res = get_coverage(jds)

    if len(res) == 0:
        return None

    res["date"] = Time(res["obsjd"].to_numpy(), format="jd").isot

    mns = MNS(df=res)

    mns.data.query(f"obsjd >= {t_min.jd} and obsjd <= {t_max.jd}", inplace=True)

    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

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
        logger.info(f"Found {len(set(mns.data))} observations in total.")
    else:
        logger.warning("Found no observations on IPAC depot or TAP.")

    return mns
