#!/usr/bin/env python3
# License: BSD-3-Clause

import os, time
from astropy.time import Time
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

# For absolute magnitude calculation
GENERIC_COSMOLOGY = FlatLambdaCDM(H0=70, Om0=0.3)


def lightcurve_from_alert(
        alert: dict,
        figsize: list=[6.47, 4],
        title: str=None,
        include_ulims: bool=False,
        mag_range: list=None,
        z: float=None,
        logger=None,
    ): 
    """ plot AMPEL alerts as lightcurve """

    if np.isnan(z):
        z = None

    # ZTF color and naming scheme
    BAND_NAMES = {1: "ZTF g", 2: "ZTF r", 3: "ZTF i"}
    BAND_COLORS = {1: "green", 2: "red", 3: "orange"}

    name = alert[0]["objectId"]
    candid = alert[0]["candidate"]
    prv_candid = alert[0]["prv_candidates"]

    if logger is not None:
        logger.debug(f"Plotting {name}")

    df = pd.DataFrame(candid, index=[0])
    df_ulims = pd.DataFrame()

    # Filter out images with negative difference flux
    i = 0
    for prv in prv_candid:
        # Go through the alert history
        if "magpsf" in prv.keys() and "isdiffpos" in prv.keys():
            i+=1
            ser = pd.Series(prv, name=i)
            df = df.append(ser)
        else:
            df_ulims = df_ulims.append(prv, ignore_index=True)
            i+=1

    df["mjd"] = df["jd"] - 2400000.5
    df_ulims["mjd"] = df_ulims["jd"] - 2400000.5

    # Helper functions for the axis conversion (from MJD to days from today)
    def t0_dist(obsmjd):
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return obsmjd - t0

    def t0_to_mjd(dist_to_t0):
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return t0 + dist_to_t0

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    fig.subplots_adjust(top=0.8)
    ax2 = ax.secondary_xaxis("top", functions=(t0_dist, t0_to_mjd))

    # If redshift is given, calculate absolute magnitude via luminosity distance
    # and plot as right axis
    if z is not None:

        dist_l = GENERIC_COSMOLOGY.luminosity_distance(z).to(u.pc).value

        def mag_to_absmag(mag):
            absmag = mag - 5 * (np.log10(dist_l) - 1)
            return absmag

        def absmag_to_mag(absmag):
            mag = absmag + 5 * (np.log10(dist_l) - 1)
            return mag

        ax3 = ax.secondary_yaxis("right", functions=(mag_to_absmag, absmag_to_mag))
        ax3.set_ylabel(f"Absolute Magnitude [AB]")

    # Get time now as UTC time
    ts = time.time()
    utc_now = datetime.utcfromtimestamp(ts)
    utc_string = utc_now.strftime("%Y-%m-%d")
    ax2.set_xlabel(f"Days from {utc_string}")

    # Give the figure a title
    if title is None:
        fig.suptitle(f"{name}", fontweight="bold")
    else:
        fig.suptitle(title, fontweight="bold")

    # grid line every 100 days
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.grid(b=True, axis="both", alpha=0.5)
    ax.set_xlabel("MJD")
    ax.set_ylabel("Magnitude [AB]")

    if mag_range is None:
        ax.set_ylim([23, 15])
    else:
        ax.set_ylim([np.max(mag_range), np.min(mag_range)])

    for fid in BAND_NAMES.keys():

        # Plot datapoints
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

        # Plot upper limits
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

    if z is not None:
        axes = [ax, ax2, ax3]
    else:
        axes = [ax, ax2]

    return fig, axes
