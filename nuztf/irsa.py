#!/usr/bin/env python3
# coding: utf-8

import os
import logging

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from astropy.time import Time
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.exceptions import RemoteServiceError
from astroquery.ipac.ned import Ned

from ztfquery.io import LOCALSOURCE
from ztfquery.lightcurve import LCQuery

from nuztf.ampel_api import ampel_api_name
from nuztf.style import plot_dir, big_fontsize, base_width, base_height, dpi
from nuztf.utils import cosmo, is_ztf_name, is_tns_name, query_tns_by_name
from nuztf.observations import get_most_recent_obs
from nuztf.parse_nu_gcn import find_gcn_no, parse_gcn_circular

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


def load_irsa(ra_deg: float, dec_deg: float, radius_arcsec: float = 0.5, **kwargs):
    """
    Get lightcuve from IPAC
    """

    logger.debug("Querying IPAC")
    df = LCQuery.from_position(ra_deg, dec_deg, radius_arcsec, **kwargs).data

    logger.debug(f"Found {len(df)} datapoints")

    if len(df) == 0:
        logger.info("No data found.")
        return None

    else:
        mask = df.catflags > 0

        flags = list(set(df.catflags))

        logger.info(
            f"Found {len(df)} datapoints, masking {np.sum(mask)} datapoints with bad flags."
        )

        for flag in sorted(flags):
            logger.debug(f"{np.sum(df.catflags == flag)} datapoints with flag {flag}")

        df = df.drop(df[mask].index)
        return df


def plot_irsa_lightcurve(
    source_name: str,
    nu_name: list = None,
    source_coords: list = None,
    source_redshift: float = None,
    plot_mag: bool = False,
    atel: bool = True,
    plot_folder: str = plot_dir,
    extra_folder: str = None,
    check_obs: bool = True,
    check_obs_lookback_weeks: float = 4,
    from_cache: bool = False,
    cache_dir: str = os.path.join(LOCALSOURCE, "cache/"),
    expanded_labels: bool = True,
    ylim: tuple = None,
    radius_arcsec: float = 0.5,
    query_irsa_for_logs: bool = True,
) -> None:
    plot_title = source_name

    # If there are no coordinates, try name resolve to get coordinates!

    if source_coords is None:

        logger.info(f"Trying to resolve {source_name} and query for coordinates.")
        logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

        # Try ampel to find ZTF coordinates

        if is_ztf_name(name=source_name):
            logger.info("Source name is a ZTF name.")
            res = ampel_api_name(source_name, with_history=False)[0]
            source_coords = [res["candidate"]["ra"], res["candidate"]["dec"]]
            logger.info(f"Found ZTF coordinates for source {source_name}")

        # Try TNS

        elif is_tns_name(name=source_name):
            logger.info("Source name is a TNS name.")
            result_dict = query_tns_by_name(name=source_name)

            if not result_dict:
                logger.warning(f"{source_name} is not in TNS.")

            if result_dict:
                logger.info(f"Found {source_name} on TNS.")
                res = result_dict["data"]["reply"]
                ra = res["radeg"]
                dec = res["decdeg"]
                source_coords = [ra, dec]
                if np.logical_and("redshift" in res.keys(), source_redshift is None):
                    source_redshift = res["redshift"]

        # Otherwise try NED

        else:
            logger.info(
                "Source name is neither as ZTF, nor a TNS name. Querying NED instead."
            )
            result_table = Ned.query_object(source_name)

            if len(result_table) == 0:
                logger.warning(
                    f"Failed to resolve name {source_name} in NED. Trying to be clever instead."
                )

                querystring = "".join(
                    [
                        x
                        for x in source_name
                        if x in [str(i) for i in range(10)] + ["+", "-"]
                    ]
                )

                result_table = Ned.query_object(
                    "".join(
                        [
                            x
                            for x in source_name
                            if x in [str(i) for i in range(10)] + ["+", "-"]
                        ]
                    )
                )

            if len(result_table) == 1:
                source_coords = [result_table["RA"][0], result_table["DEC"][0]]

                if "ZTF" in plot_title:
                    plot_title += f' ({result_table["Object Name"][0]})'

                if np.logical_and(
                    str(result_table["Redshift"][0]) != "--", source_redshift is None
                ):
                    source_redshift = result_table["Redshift"]

                logger.info(
                    f"Using AStroquery NED query result for name {source_name} ({source_coords})"
                )

            if source_coords is None:
                sc = SkyCoord.from_name(source_name)
                logger.info(
                    f"Using Astroquery CDS query result for name {source_name} (RA={sc.ra}, Dec={sc.dec})"
                )
                source_coords = (sc.ra.value, sc.dec.value)

    # Try to find a catalogue source nearby using coordinates

    if ("ZTF" in source_name) and source_coords:

        c = SkyCoord(source_coords[0], source_coords[1], unit=u.deg, frame="icrs")

        r = 0.5 * u.arcsecond

        logger.info("Querying NED")
        result_table = Ned.query_region(c, radius=r)

        if len(result_table) == 1:
            if "ZTF" in plot_title:
                plot_title += f' ({result_table["Object Name"][0]})'

            source_coords = [result_table["RA"][0], result_table["DEC"][0]]

            if np.logical_and(
                str(result_table["Redshift"][0]) != "--", source_redshift is None
            ):
                source_redshift = result_table["Redshift"]

            logger.info(
                f"Found likely match to {source_name}"
                f"(type = '{result_table['Type'][0]}'. "
                f"distance = {result_table['Separation'][0]} arcsec')"
            )
        elif len(result_table) > 1:
            logger.warning(
                f"Found multiple possible cross-matches: {result_table['Object Name']}"
            )
        else:
            logger.info("No NED crossmatch found.")

    # Query IRSA, or load from cache

    try:
        os.makedirs(cache_dir)
    except OSError:
        pass

    cache_path = os.path.join(cache_dir, f'{source_name.replace(" ", "")}.csv')

    if from_cache:
        logger.debug(f"Loading from {cache_path}")
        df = pd.read_csv(cache_path)

    else:
        logger.info("Querying IPAC for a lightcurve")
        df = load_irsa(source_coords[0], source_coords[1], radius_arcsec=radius_arcsec)

        logger.debug(f"Saving to {cache_path}")
        df.to_csv(cache_path)

    data = Table.from_pandas(df)

    logger.info(f"There are a total of {len(data)} detections for {source_name}")

    # Start Figure

    plt.figure(figsize=(base_width * 1.1, base_height), dpi=dpi)

    if expanded_labels:

        ax2 = plt.subplot(111)

        ax = ax2.twiny()

    else:
        ax = plt.subplot(111)

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
        except (RemoteServiceError, IndexError) as e:
            logger.info("No redshift found")

    elif np.isnan(source_redshift):
        source_redshift = None

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

    # If you want, you can check the most recent observation

    if check_obs:

        logger.info(
            f"Getting most recent ZTF observation, looking back {check_obs_lookback_weeks} weeks."
        )
        mro = get_most_recent_obs(
            ra=source_coords[0],
            dec=source_coords[1],
            lookback_weeks_max=check_obs_lookback_weeks,
        )

        if mro is not None:
            ot = format_date(Time(mro["obsjd"], format="jd"), atel=atel)
            logger.info(f"Most recent observation at {ot}")
        else:
            logger.info("No recent observation found.")

    # Plot each band (g/r/i)

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

    # You can force the y limits if you want

    if ylim is not None:
        ax.set_ylim(ylim)

    if plot_mag:
        ax.set_ylabel(r"Apparent magnitude [AB]", fontsize=big_fontsize)

        ax.invert_yaxis()

        if source_redshift is not None:
            ax1b.set_ylabel(rf"Absolute magnitude [AB]", fontsize=big_fontsize)

            y_min, y_max = ax.get_ylim()

            ax1b.invert_yaxis()

            ax1b.set_ylim(y_min - dist_mod, y_max - dist_mod)

    else:
        ax.set_ylabel(r"$\nu$F$_{\nu}$ [erg cm$^{-2}$ s$^{-1}$]", fontsize=big_fontsize)

        ax.set_yscale("log")

        if source_redshift is not None:
            ax1b.set_ylabel(r"$\nu$L$_{\nu}$ [erg s$^{-1}$]", fontsize=big_fontsize)
            ax1b.set_yscale("log")

            y_min, y_max = ax.get_ylim()

            ax1b.set_ylim(
                y_min * conversion_factor.value, y_max * conversion_factor.value
            )

    ax.set_xlabel("Date [MJD]", fontsize=big_fontsize)

    # Add neutrino

    if nu_name is None:
        nu_name = []

    if not isinstance(nu_name, list):
        nu_name = [nu_name]

    for j, nu in enumerate(nu_name):
        gcn_no = find_gcn_no(nu)
        gcn_info = parse_gcn_circular(gcn_no)

        ax.axvline(gcn_info["time"].mjd, linestyle=":", label=nu, color=f"C{j}")

    if expanded_labels:

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

        ax2.tick_params(axis="both", which="major", labelsize=big_fontsize)

    ax.tick_params(axis="both", which="both", labelsize=big_fontsize, left=True)

    if source_redshift is not None:
        ax1b.tick_params(axis="both", which="major", labelsize=big_fontsize)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.22 + 0.2 * float(expanded_labels)),
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
        plt.savefig(extra_path, bbox_inches="tight", pad_inches=0.05)
