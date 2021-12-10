#!/usr/bin/env python3
# coding: utf-8

import os
import logging

import matplotlib.pyplot as plt
import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astroquery.exceptions import RemoteServiceError
from astroquery.ned import Ned
from ztfquery.lightcurve import LCQuery
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM

from nuztf.style import plot_dir, big_fontsize, base_width, base_height, dpi
from nuztf.observation_log import get_most_recent_obs
from nuztf.parse_nu_gcn import find_gcn_no, parse_gcn_circular

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

logger = logging.getLogger(__name__)


def format_date(t, atel=True):
    t.format = "fits"

    if atel:
        frac_days = f"{t.mjd - int(t.mjd):.2f}"[1:]
        t.out_subfmt = "date"

        dt = "".join([t.value, frac_days])
    else:
        dt = t.value

    return dt


def plot_irsa_lightcurve(
    source_name: str,
    nu_name: list = None,
    source_coords: list = None,
    source_redshift: float = None,
    plot_mag: bool = False,
    atel: bool = True,
    plot_folder: str = plot_dir,
    extra_folder: str = None,
    logger=None,
    check_obs=True,
    check_obs_lookback_weeks=4,
):
    plot_title = source_name

    if logger is None:
        logger = logging.getLogger(__name__)
    else:
        logger = logger

    if np.logical_and("ZTF" in source_name, source_coords is not None):

        c = SkyCoord(source_coords[0], source_coords[1], unit=u.deg, frame='icrs')

        r = 0.5 * u.arcsecond

        result_table = Ned.query_region(c, radius=r)

        if len(result_table) == 1:
            if "ZTF" in plot_title:
                plot_title += f' ({result_table["Object Name"][0]})'

            source_coords = [result_table["RA"][0], result_table["DEC"][0]]

            if str(result_table["Redshift"][0]) != "--":
                source_redshift = result_table["Redshift"]

            logger.info(f"Found likely match to {source_name}"
                        f"(type = '{result_table['Type'][0]}'. "
                        f"distance = {result_table['Separation'][0]} arcsec')")
        elif len(result_table) > 1:
            logger.warning(f"Found multiple possible cross-matches: {result_table['Object Name']}")
        else:
            logger.info("No NED crossmatch found.")

    if source_coords is None:

        result_table = Ned.query_object(source_name)
        if "ZTF" in plot_title:
            plot_title += f' ({result_table["Object Name"][0]})'
        source_coords = [result_table["RA"][0], result_table["DEC"][0]]

        if str(result_table["Redshift"][0]) != "--":
            source_redshift = result_table["Redshift"]

        # sc = SkyCoord.from_name(source_name)
        logger.info(
            f"Using Astropy NED query result for name {source_name} ({source_coords})"
        )

    df = LCQuery.from_position(source_coords[0], source_coords[1], 1.0).data

    data = Table.from_pandas(df)

    logger.info(f"There are a total of {len(data)} detections for {source_name}")

    plt.figure(figsize=(base_width, base_height), dpi=dpi)

    ax2 = plt.subplot(111)

    ax = ax2.twiny()

    # If you have a redshift, you can add a second y axis!

    if source_redshift is None:
        logger.info("Querying NED to check for a redshift")
        try:
            result_table = Ned.query_object(source_name)
            if len(result_table["Redshift"]) == 1:

                if str(result_table["Redshift"][0]) == "--":
                    raise RemoteServiceError

                source_redshift = result_table["Redshift"][0]
                logger.info(f"Found a redshift of {source_redshift}")
            elif len(result_table["Redshift"]) > 1:
                logger.warning(f"Found multiple redshifts: {result_table}")
            else:
                raise RemoteServiceError
        except RemoteServiceError:
            logger.info("No redshift found")

    if source_redshift is not None:

        ax1b = ax.twinx()

        redshift = 1.0 + source_redshift

        if plot_mag:
            dist_mod = 5 * (
                np.log10(cosmo.luminosity_distance(z=(redshift - 1)).to(u.pc).value)
                - 1.0
            )
        else:
            conversion_factor = (
                4
                * np.pi
                * cosmo.luminosity_distance(z=(redshift - 1)).to(u.cm) ** 2.0
                / (redshift)
            )

    cmap = {"zg": "g", "zr": "r", "zi": "orange"}

    wl = {
        "zg": 472.27,
        "zr": 633.96,
        "zi": 788.61,
    }

    markersize = 2.0

    latest_index = list(data["mjd"]).index(max(data["mjd"]))
    latest = data[latest_index]

    dt = format_date(Time(latest["mjd"], format="mjd"), atel=atel)

    logger.info(
        f"Most recent detection on {dt} UT at a magnitude of "
        f"{latest['filtercode'][1]}={latest['mag']:.2f}+/-{latest['magerr']:.2f}"
    )

    if check_obs:

        mro = get_most_recent_obs(
            ra=source_coords[0],
            dec=source_coords[1],
            lookback_weeks_max=check_obs_lookback_weeks,
            logger=logger,
        )

        if mro is not None:
            ot = format_date(Time(mro["obsjd"], format="jd"), atel=atel)
            logger.info(f"Most recent observation at {ot}")
        else:
            logger.info("No recent observation found.")

    for fc in ["zg", "zr", "zi"]:
        mask = data["filtercode"] == fc

        mags = data["mag"][mask] * u.ABmag

        magerrs = (data["magerr"][mask] + data["mag"][mask]) * u.ABmag

        if plot_mag:
            ax.errorbar(
                data["mjd"][mask],
                mags.value,
                yerr=data["magerr"][mask],
                marker="o",
                linestyle=" ",
                markersize=markersize,
                c=cmap[fc],
                label=f"{fc[-1]} ({wl[fc]:.0f} nm)",
            )

            if source_redshift is not None:
                ax1b.errorbar(
                    data["mjd"][mask],
                    mags.value - dist_mod,
                    yerr=data["magerr"][mask],
                    marker="o",
                    linestyle=" ",
                    markersize=markersize,
                    c=cmap[fc],
                    label=f"{fc[-1]} ({wl[fc]:.0f} nm)",
                )

        else:

            flux_j = mags.to(u.Jansky)

            f = (const.c / (wl[fc] * u.nm)).to("Hz")

            flux = (flux_j * f).to("erg cm-2 s-1")

            jerrs = magerrs.to(u.Jansky)
            ferrs = (jerrs * f).to("erg cm-2 s-1").value - flux.value

            ax.errorbar(
                data["mjd"][mask],
                flux.to("erg cm-2 s-1").value,
                yerr=ferrs,
                marker="o",
                linestyle=" ",
                markersize=markersize,
                c=cmap[fc],
                label=f"{fc[-1]} ({wl[fc]:.0f} nm)",
            )

            if source_redshift is not None:
                l = flux * conversion_factor

                ax1b.errorbar(
                    data["mjd"][mask],
                    l.to("erg s-1"),
                    marker="o",
                    linestyle=" ",
                    markersize=markersize,
                    c=cmap[fc],
                    label=f"{fc[-1]} ({wl[fc]:.0f} nm)",
                )

    if plot_mag:
        ax2.set_ylabel(r"Apparent magnitude [AB]", fontsize=big_fontsize)

        ax2.invert_yaxis()

        if source_redshift is not None:
            ax1b.set_ylabel(fr"Absolute magnitude [AB]", fontsize=big_fontsize)

            y_min, y_max = ax.get_ylim()

            ax1b.invert_yaxis()

            ax1b.set_ylim(y_min - dist_mod, y_max - dist_mod)

    else:
        ax2.set_ylabel(
            r"$\nu$F$_{\nu}$ [erg cm$^{-2}$ s$^{-1}$]", fontsize=big_fontsize
        )

        ax2.set_yscale("log")

        if source_redshift is not None:
            ax1b.set_ylabel(r"$\nu$L$_{\nu}$ [erg s$^{-1}$]", fontsize=big_fontsize)
            ax1b.set_yscale("log")

            y_min, y_max = ax.get_ylim()

            ax1b.set_ylim(
                y_min * conversion_factor.value, y_max * conversion_factor.value
            )

    ax.set_xlabel("Date (MJD)", fontsize=big_fontsize)

    # Add neutrino

    if nu_name is None:
        nu_name = []

    if not isinstance(nu_name, list):
        nu_name = [nu_name]

    for j, nu in enumerate(nu_name):
        gcn_no = find_gcn_no(nu)
        gcn_info = parse_gcn_circular(gcn_no)

        ax.axvline(gcn_info["time"].mjd, linestyle=":", label=nu, color=f"C{j}")

    # Set up ISO dates

    lmjd, umjd = ax.get_xlim()

    lt = Time(lmjd, format="mjd")
    ut = Time(umjd, format="mjd")

    nt = Time.now()
    nt.format = "fits"

    mjds = []
    labs = []

    for year in range(2016, int(nt.value[:4]) + 1):
        for k, month in enumerate([1, 7]):

            t = Time(f"{year}-{month}-01T00:00:00.01", format="isot", scale="utc")
            t.format = "fits"
            t.out_subfmt = "date"

            if np.logical_and(t > lt, t < ut):
                mjds.append(t.mjd)
                labs.append(t.value)

    ax2.set_xticks(mjds)
    ax2.set_xticklabels(labels=labs, rotation=80)

    ax2.set_xlim(lmjd, umjd)

    ax.set_title(f'ZTF Lightcurve of {plot_title.replace("J", " J")}', y=1.4)

    ax.tick_params(axis="both", which="major", labelsize=big_fontsize)
    ax2.tick_params(axis="both", which="major", labelsize=big_fontsize)

    if source_redshift is not None:
        ax1b.tick_params(axis="both", which="major", labelsize=big_fontsize)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.42),
        ncol=3 + len(nu_name),
        fancybox=True,
        fontsize=big_fontsize,
    )

    filename = f"{source_name.replace(' ', '')}_lightcurve{['_flux', ''][plot_mag]}.png"

    output_path = os.path.join(plot_folder, f"{filename}")

    logger.info(f"Saving to {output_path}")

    plt.savefig(output_path, bbox_inches="tight", pad_inches=0.05)

    if extra_folder is not None:
        extra_path = os.path.join(extra_folder, f"{filename}")
        logger.info(f"Saving to {extra_path}")
        plt.savefig(extra_path, bbox_inches="tight", pad_inches=0.5)
