import glob

import pyvo.dal
import time
from typing import Optional
from nuztf import credentials
from pyvo.auth import securitymethods, authsession
import requests
import pandas as pd
import logging
from glob import glob
import backoff

import os
import warnings
import numpy as np
from astropy.time import Time
from astropy import units as u
from ztfquery import skyvision
from ztfquery.io import LOCALSOURCE
from ztfquery.fields import get_fields_containing_target

logger = logging.getLogger(__name__)

username, password = credentials.load_credentials("irsa")

# Gather login information
data = {"username": username, "password": password}

headers = {"Content-Type": "application/x-www-form-urlencoded", "Accept": "text/plain"}

login_url = "https://irsa.ipac.caltech.edu"

coverage_dir = os.path.join(LOCALSOURCE, "all_obs")

try:
    os.makedirs(coverage_dir)
except OSError:
    pass

partial_flag = "_PARTIAL"


def coverage_path(jd: float) -> str:

    if (Time.now().jd - jd) < 1:
        partial_ext = partial_flag
    else:
        partial_ext = ""

    return os.path.join(coverage_dir, f"{jd}{partial_ext}.csv")


class MNS:
    def __init__(self, df):
        self.data = df


@backoff.on_exception(backoff.expo, requests.exceptions.RequestException, max_time=60)
def write_coverage(jds: [int]):
    with requests.Session() as session:
        # Create a session and do the login.
        # The cookie will end up in the cookie jar of the session.
        response = session.post(login_url, data=data, headers=headers)
        response.raise_for_status()
        auth = authsession.AuthSession()
        auth.credentials.set(securitymethods.ANONYMOUS, session)
        client = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP", auth)

        for jd in jds:
            obstable = client.search(
                f"""
            SELECT field,rcid,fid,expid,obsjd,exptime,maglimit,ipac_gid,seeing
            FROM ztf.ztf_current_meta_sci WHERE (obsjd BETWEEN {jd} AND {jd + 1})
            """
            ).to_table()
            names = ("ipac_gid",)
            renames = ("programid",)
            obstable.rename_columns(names, renames)

            obs_grouped_by_exp = obstable.to_pandas().groupby("expid")

            output_path = coverage_path(jd)

            with open(output_path, "w") as fid:

                fid.write(
                    "obsid,field,obsjd,datetime,seeing,limmag,exposure_time,fid,processed_fraction\n"
                )

                for group_name, df_group in obs_grouped_by_exp:
                    processed_fraction = len(df_group["field"]) / 64.0

                    t = Time(df_group["obsjd"].iloc[0], format="jd")
                    t.format = "isot"

                    line = f'{int(df_group["expid"].iloc[0])},{int(df_group["field"].iloc[0])},{t.jd},{t},{df_group["seeing"].median():.10f},{df_group["maglimit"].median():.10f},{df_group["exptime"].iloc[0]},{int(df_group["fid"].iloc[0])},{processed_fraction} \n'
                    fid.write(line)

            time.sleep(1)


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

    logger.info("Getting observation logs from skyvision.")
    mns = get_obs_summary_skyvision(t_min, t_max, max_days=max_days)

    if len(mns.data) == 0:
        logger.debug("Empty observation log, try IRSA instead.")
        mns = get_obs_summary_irsa(t_min, t_max, max_days=max_days)

    logger.debug(f"Found {len(mns.data)} observations in total.")

    return mns


def get_obs_summary_irsa(t_min, t_max=None, max_days: int = None):
    """
    Get observation summary from IRSA
    """

    if t_max is None:
        if max_days is None:
            t_max = Time.now()
        else:
            t_max = t_min + (max_days * u.day)

    jds = np.arange(int(t_min.jd), int(t_max.jd) + 1)

    logger.debug("Getting coverage")

    df = get_coverage(jds)

    mns = MNS(df)

    mns.data.query(f"obsjd >= {t_min.jd} and obsjd <= {t_max.jd}", inplace=True)
    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    logger.debug("Done")

    return mns


def get_obs_summary_skyvision(t_min, t_max=None, max_days: int = None):
    """
    Get observation summary from Skyvision
    """

    t_min_jd = t_min.jd
    t_min_date = t_min.to_value("iso", subfmt="date")

    if t_max is None:
        if max_days is None:
            t_max = Time.now()
        else:
            t_max = t_min + max_days * u.day

    elif t_max is not None and max_days is not None:
        raise ValueError("Choose either t_max or max_days, not both")

    t_max_jd = t_max.jd
    t_max_date = t_max.to_value("iso", subfmt="date")

    # ztfquery saves nightly observations in a cache, and does not redownload them.
    # If the nightly log was not complete, it will never be updated.
    # Here we simply clear the cache and cleanly re-download everything.

    logger.debug(
        f"Obtaining nightly observation logs from {t_min_date} to {t_max_date}"
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

    mns.data.query(f"obsjd >= {t_min_jd} and obsjd <= {t_max_jd}", inplace=True)

    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    logger.debug(f"Found {len(mns.data)} observations in total.")

    return mns


def get_most_recent_obs(ra: float, dec: float, lookback_weeks_max: int = 12):
    """ """

    fields = get_fields_containing_target(ra, dec)._data

    logger.info(f"Target in fields {fields}")

    mask = 0.0

    t_max = Time.now()

    lookback_weeks = 1

    while np.sum(mask) < 1 and lookback_weeks <= lookback_weeks_max:

        t_min = t_max - 1 * u.week

        if lookback_weeks > 1:
            logger.debug(
                f"Searching for most recent obs during the last {lookback_weeks_max} weeks. Looking back {lookback_weeks} weeks."
            )
        else:
            logger.debug(
                f"Searching for most recent obs during the last {lookback_weeks_max} weeks. Looking back {lookback_weeks} week."
            )

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)

            mns = get_obs_summary(t_min=t_min, t_max=t_max)

        mask = np.array([x in fields for x in mns.data["field"]])

        t_max = t_min

        lookback_weeks += 1

    if np.sum(mask) >= 1:
        index = list(mns.data["datetime"]).index(max(mns.data["datetime"][mask]))
        mro = mns.data.iloc[index]
        return mro

    else:
        logger.warning(
            f"Found no observation during the last {lookback_weeks_max} weeks."
        )
        return None
