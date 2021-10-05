#!/usr/bin/env python3
# License: BSD-3-Clause

import os, time
from astropy.time import Time
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def lightcurve_from_alert(
		alert: dict,
		figsize: list=[6.47, 4],
		title: str=None,
		include_ulims: bool=False,
		mag_range: list=None,
		logger=None,
	):
    """ plot AMPEL alerts as lightcurve """


    BAND_NAMES = {1: "ZTF g", 2: "ZTF r", 3: "ZTF i"}
    BAND_COLORS = {1: "green", 2: "red", 3: "orange"}

    name = alert[0]["objectId"]
    candid = alert[0]["candidate"]
    prv_candid = alert[0]["prv_candidates"]

    if logger is not None:
        logger.debug(f"Plotting {name}")

    df = pd.DataFrame(candid, index=[0])
    df_ulims = pd.DataFrame()

    i = 0
    for prv in prv_candid:
        if "magpsf" in prv.keys() and "isdiffpos" in prv.keys():
            i+=1
            ser = pd.Series(prv, name=i)
            df = df.append(ser)
        else:
            df_ulims = df_ulims.append(prv, ignore_index=True)
            i+=1

    df["mjd"] = df["jd"] - 2400000.5
    df_ulims["mjd"] = df_ulims["jd"] - 2400000.5

    def t0_dist(obsmjd):
        """ """
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return obsmjd - t0

    def t0_to_mjd(dist_to_t0):
        """ """
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return t0 + dist_to_t0

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    fig.subplots_adjust(top=0.8)
    ax2 = ax.secondary_xaxis("top", functions=(t0_dist, t0_to_mjd))

    # Get time now as UTC time
    ts = time.time()
    utc_now = datetime.utcfromtimestamp(ts)
    utc_string = utc_now.strftime("%Y-%m-%d")
    ax2.set_xlabel(f"Days from {utc_string}")

    if title is None:
        fig.suptitle(f"{name}", fontweight="bold")
    else:
        fig.suptitle(title, fontweight="bold")
    ax.grid(b=True, axis="both")
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude [AB]")

    if mag_range is None:
        ax.set_ylim([23, 15])
    else:
        ax.set_ylim([np.max(mag_range), np.min(mag_range)])

    for fid in BAND_NAMES.keys():
        df_temp = df.query("fid == @fid")
        ax.errorbar(
            df_temp["mjd"],
            df_temp["magpsf"],
            df_temp["sigmapsf"],
            color=BAND_COLORS[fid],
            fmt=".",
            label=BAND_NAMES[fid],
            mec="black",
            mew=0.5,
        )

        if include_ulims:
            df_temp2 = df_ulims.query("fid == @fid")
            ax.scatter(
                df_temp2["mjd"],
                df_temp2["diffmaglim"],
                c=BAND_COLORS[fid],
                marker="v",
                s=1.3,
                alpha=0.5,
            )

    plt.tight_layout()

    return fig, [ax, ax2]