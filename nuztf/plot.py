#!/usr/bin/env python3
# License: BSD-3-Clause

import os, time, gzip, io
from datetime import datetime
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.colors import Normalize
from base64 import b64decode

from astropy.time import Time
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.io import fits
from astropy import visualization
from ztfquery.utils.stamps import get_ps_stamp

from nuztf.cat_match import get_cross_match_info
from nuztf.ampel_api import ensure_cutouts

# For absolute magnitude calculation
GENERIC_COSMOLOGY = FlatLambdaCDM(H0=70, Om0=0.3)


def lightcurve_from_alert(
    alert: list,
    # figsize: list=[6.47, 4],
    figsize: list = [8, 5],
    title: str = None,
    include_ulims: bool = True,
    include_cutouts: bool = True,
    include_crossmatch: bool = True,
    mag_range: list = None,
    z: float = None,
    legend: bool = False,
    grid_interval: int = None,
    logger=None,
):
    """plot AMPEL alerts as lightcurve"""

    if logger is None:
        import logging

        logger = logging.getLogger(__name__)
    else:
        logger = logger

    if z is not None:
        if np.isnan(z):
            z = None
            logger.debug("Redshift is nan, will be ignored")

    # ZTF color and naming scheme
    BAND_NAMES = {1: "ZTF g", 2: "ZTF r", 3: "ZTF i"}
    BAND_COLORS = {1: "green", 2: "red", 3: "orange"}

    name = alert[0]["objectId"]
    candidate = alert[0]["candidate"]
    prv_candid = alert[0]["prv_candidates"]

    if include_cutouts:
        if "cutoutScience" in alert[0].keys():
            if "stampData" in alert[0]["cutoutScience"].keys():
                logger.debug(f"{name}: Cutouts are present.")
            else:
                logger.debug(f"{name}: Cutouts are missing data. Will obtain them")
                alert = ensure_cutouts(alert, logger=logger)
        else:
            logger.debug(
                "The alert dictionary does not contain cutouts. Will obtain them."
            )
            alert = ensure_cutouts(alert, logger=logger)

    logger.debug(f"Plotting {name}")

    df = pd.DataFrame(candidate, index=[0])
    df_ulims = pd.DataFrame()

    # Filter out images with negative difference flux
    i = 0
    for prv in prv_candid:
        # Go through the alert history
        if "magpsf" in prv.keys() and "isdiffpos" in prv.keys():
            i += 1
            ser = pd.Series(prv, name=i)
            df = df.append(ser)
        else:
            df_ulims = df_ulims.append(prv, ignore_index=True)
            i += 1

    df["mjd"] = df["jd"] - 2400000.5
    df_ulims["mjd"] = df_ulims["jd"] - 2400000.5

    # Helper functions for the axis conversion (from MJD to days from today)
    def t0_dist(obsmjd):
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return obsmjd - t0

    def t0_to_mjd(dist_to_t0):
        t0 = Time(time.time(), format="unix", scale="utc").mjd
        return t0 + dist_to_t0

    def mjd_to_date(mjd):
        # return mjd + 20000
        return Time(mjd, format="mjd").mjd

    def date_to_mjd(date):
        return Time(date, format="mjd").mjd

    fig = plt.figure(figsize=figsize)

    if include_cutouts:
        lc_ax1 = fig.add_subplot(5, 4, (9, 19))
        cutoutsci = fig.add_subplot(5, 4, (1, 5))
        cutouttemp = fig.add_subplot(5, 4, (2, 6))
        cutoutdiff = fig.add_subplot(5, 4, (3, 7))
        cutoutps1 = fig.add_subplot(5, 4, (4, 8))
    else:
        lc_ax1 = fig.add_subplot(1, 1, 1)
        fig.subplots_adjust(top=0.8, bottom=0.15)

    plt.subplots_adjust(wspace=0.4, hspace=1.8)

    if include_cutouts:
        for cutout_, ax_, type_ in zip(
            [alert[0], alert[0], alert[0]],
            [cutoutsci, cutouttemp, cutoutdiff],
            ["Science", "Template", "Difference"],
        ):
            create_stamp_plot(alert=cutout_, ax=ax_, type=type_)

        img = get_ps_stamp(
            candidate["ra"], candidate["dec"], size=240, color=["y", "g", "i"]
        )
        cutoutps1.imshow(np.asarray(img))
        cutoutps1.set_title("PS1", fontdict={"fontsize": "small"})
        cutoutps1.set_xticks([])
        cutoutps1.set_yticks([])

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

        lc_ax3 = lc_ax1.secondary_yaxis(
            "right", functions=(mag_to_absmag, absmag_to_mag)
        )

        if not include_cutouts:
            lc_ax3.set_ylabel(f"Absolute Magnitude [AB]")

    # Give the figure a title
    if not include_cutouts:
        if title is None:
            fig.suptitle(f"{name}", fontweight="bold")
        else:
            fig.suptitle(title, fontweight="bold")

    if grid_interval is not None:
        lc_ax1.xaxis.set_major_locator(MultipleLocator(grid_interval))

    lc_ax1.grid(b=True, axis="both", alpha=0.5)
    lc_ax1.set_ylabel("Magnitude [AB]")

    if not include_cutouts:
        lc_ax1.set_xlabel("MJD")

    # Determine magnitude limits
    if mag_range is None:
        max_mag = np.max(df.magpsf.values) + 0.3
        min_mag = np.min(df.magpsf.values) - 0.3
        lc_ax1.set_ylim([max_mag, min_mag])
    else:
        lc_ax1.set_ylim([np.max(mag_range), np.min(mag_range)])

    for fid in BAND_NAMES.keys():

        # Plot older datapoints
        df_temp = df.iloc[1:].query("fid == @fid")
        lc_ax1.errorbar(
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
            lc_ax1.scatter(
                df_temp2["mjd"],
                df_temp2["diffmaglim"],
                c=BAND_COLORS[fid],
                marker="v",
                s=1.3,
                alpha=0.5,
            )

    # Plot datapoint from alert
    df_temp = df.iloc[0]
    fid = df_temp["fid"]
    lc_ax1.errorbar(
        df_temp["mjd"],
        df_temp["magpsf"],
        df_temp["sigmapsf"],
        color=BAND_COLORS[fid],
        fmt=".",
        label=BAND_NAMES[fid],
        mec="black",
        mew=0.5,
        markersize=12,
    )

    if legend:
        plt.legend()

    # Now we create an infobox
    if include_cutouts:
        info = []

        info.append(name)
        info.append("------------------------")
        info.append(f"RA: {candidate['ra']:.8f}")
        info.append(f"Dec: {candidate['dec']:.8f}")
        if "drb" in candidate.keys():
            info.append(f"drb: {candidate['drb']:.3f}")
        else:
            info.append(f"rb: {candidate['rb']:.3f}")
        info.append("------------------------")

        for entry in ["sgscore1", "distpsnr1", "srmag1"]:
            info.append(f"{entry[:-1]}: {candidate[entry]:.3f}")
            # for k in [k for k in candidate.keys() if kk in k]:
            #     info.append(f"{k}: {candidate.get(k):.3f}")

        fig.text(0.77, 0.55, "\n".join(info), va="top", fontsize="medium", alpha=0.5)

    if include_crossmatch:
        xmatch_info = get_cross_match_info(alert[0])
        if include_cutouts:
            ypos = 0.975
        else:
            ypos = 0.035

        fig.text(
            0.5,
            ypos,
            xmatch_info,
            va="top",
            ha="center",
            fontsize="medium",
            alpha=0.5,
        )

    # Ugly hack because secondary_axis does not work with astropy.time.Time datetime conversion
    mjd_min = np.min(df.mjd.values)
    mjd_max = np.max(df.mjd.values)
    length = mjd_max - mjd_min

    lc_ax1.set_xlim([mjd_min - (length / 20), mjd_max + (length / 20)])

    lc_ax2 = lc_ax1.twiny()

    datetimes = [Time(x, format="mjd").datetime for x in [mjd_min, mjd_max]]

    lc_ax2.scatter(
        [Time(x, format="mjd").datetime for x in [mjd_min, mjd_max]], [20, 20], alpha=0
    )
    lc_ax2.tick_params(axis="both", which="major", labelsize=6, rotation=45)
    lc_ax1.tick_params(axis="x", which="major", labelsize=6, rotation=45)
    lc_ax1.ticklabel_format(axis="x", style="plain")
    lc_ax1.tick_params(axis="y", which="major", labelsize=9)

    if z is not None:
        lc_ax3.tick_params(axis="both", which="major", labelsize=9)

    if z is not None:
        axes = [lc_ax1, lc_ax2, lc_ax3]
    else:
        axes = [lc_ax1, lc_ax2]

    return fig, axes


def create_stamp_plot(alert: dict, ax, type: str):
    """Helper function to create cutout subplot"""

    with gzip.open(
        io.BytesIO(b64decode(alert[f"cutout{type}"]["stampData"])), "rb"
    ) as f:
        data = fits.open(io.BytesIO(f.read()), ignore_missing_simple=True)[0].data
    vmin, vmax = np.percentile(data[data == data], [0, 100])
    data_ = visualization.AsinhStretch()((data - vmin) / (vmax - vmin))
    ax.imshow(
        data_,
        norm=Normalize(*np.percentile(data_[data_ == data_], [0.5, 99.5])),
        aspect="auto",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(type, fontdict={"fontsize": "small"})
