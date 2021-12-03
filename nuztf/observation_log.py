#!/usr/bin/env python3
# coding: utf-8

import logging, os, warnings
import numpy as np

from astropy.time import Time
from astropy import units as u
from ztfquery import skyvision
from ztfquery.io import LOCALSOURCE
from ztfquery.fields import get_fields_containing_target

logger = logging.getLogger(__name__)


def get_obs_summary(t_min, t_max=None, max_days: int = None, logger=None):
    """ """
    if logger is None:
        logger = logging.getLogger(__name__)

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


def get_most_recent_obs(
    ra: float, dec: float, lookback_weeks_max: int = 12, logger=None
):
    """ """
    if logger is None:
        logger = logging.getLogger(__name__)

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

            mns = get_obs_summary(t_min=t_min, t_max=t_max, logger=logger)

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
