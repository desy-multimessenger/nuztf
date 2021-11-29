#!/usr/bin/env python3

import datetime, logging, os
import numpy as np
from astropy.time import Time
from astropy import units as u
from ztfquery import skyvision
from ztfquery.io import LOCALSOURCE
from ztfquery.fields import get_fields_containing_target

logger = logging.getLogger(__name__)


def get_obs_summary(t_min, max_days: int = None, logger=None):
    """ """

    if logger is None:
        logger = logging.getLogger(__name__)

    mns_time = str(t_min).split("T")[0].replace("-", "")
    now = datetime.datetime.now()

    if max_days is not None:
        date_1 = datetime.datetime.strptime(mns_time, "%Y%m%d")
        end_date = date_1 + datetime.timedelta(days=max_days)
        end_date = end_date.strftime("%Y-%m-%d")
    else:
        end_date = now.strftime("%Y-%m-%d")

    start_date_jd = t_min.jd
    start_date_jd = Time(start_date_jd, format="jd").jd

    # ztfquery saves nightly observations in a cache, and does not redownload them.
    # If the nightly log was not complete, it will never be updated.
    # Here we simply clear the cache and cleanly re-download everything.

    logger.debug(f"Obtaining nightly observation logs from {mns_time} to {end_date}")

    skyvision_log = os.path.join(LOCALSOURCE, "skyvision")

    for filename in os.listdir(skyvision_log):
        if ".csv" in filename:
            path = os.path.join(skyvision_log, filename)
            os.remove(path)

    mns = skyvision.CompletedLog.from_daterange(mns_time, end=end_date, verbose=False)

    mns.data["obsjd"] = Time(list(mns.data.datetime.values), format="isot").jd

    mns.data.query(f"obsjd > {start_date_jd}", inplace=True)

    mns.data.reset_index(inplace=True)
    mns.data.drop(columns=["index"], inplace=True)

    return mns


def get_most_recent_obs(ra, dec):
    fields = get_fields_containing_target(ra, dec)._data

    logger.info(f"Target in fields {fields}")

    mask = 0.0
    day_range = 5.0

    while np.sum(mask) < 1:
        mns = get_obs_summary(Time.now() - 30.0 * u.day)

        mask = np.array([x in fields for x in mns.data["field"]])

        day_range += 5.0

    index = list(mns.data["datetime"]).index(max(mns.data["datetime"][mask]))
    mro = mns.data.iloc[index]

    return mro
