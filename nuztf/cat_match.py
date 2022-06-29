#!/usr/bin/env python
# coding: utf-8
import logging
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.ned import Ned
from astroquery.exceptions import RemoteServiceError
from astropy import units as u
from nuztf.ampel_api import ampel_api_catalog, ampel_api_name
from nuztf.utils import is_ztf_name, is_tns_name, query_tns_by_name

logger = logging.getLogger(__name__)


def query_ned_for_z(
    ra_deg: float, dec_deg: float, searchradius_arcsec: float = 20, logger=None
):
    """Function to obtain redshifts from NED (via the AMPEL API)"""

    z = None
    dist_arcsec = None

    query = ampel_api_catalog(
        catalog="NEDz_extcats",
        catalog_type="extcats",
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        searchradius_arcsec=searchradius_arcsec,
        searchtype="nearest",
        logger=logger,
    )

    if query:
        z = query["body"]["z"]
        dist_arcsec = query["dist_arcsec"]

    return z, dist_arcsec


def query_ned_astroquery(
    ra_deg: float, dec_deg: float, searchradius_arcsec: float = 0.5
):
    """
    Function to obtain NED crossmatches via astroquery
    """
    c = SkyCoord(ra_deg, dec_deg, unit=u.deg, frame="icrs")

    r = searchradius_arcsec * u.arcsecond

    try:
        return Ned.query_region(c, radius=r)
    except RemoteServiceError:
        return None


def ampel_api_tns(
    ra_deg: float, dec_deg: float, searchradius_arcsec: float = 3, logger=None
):
    """Function to query TNS via the AMPEL API"""

    full_name = None
    discovery_date = None
    source_group = None

    res = ampel_api_catalog(
        catalog="TNS",
        catalog_type="extcats",
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        searchradius_arcsec=searchradius_arcsec,
        searchtype="nearest",
        logger=logger,
    )

    if res:
        response_body = res["body"]
        name = response_body["objname"]
        prefix = response_body["name_prefix"]
        full_name = prefix + name
        discovery_date = response_body["discoverydate"]
        if "source_group" in response_body.keys():
            source_group = response_body["source_group"]["group_name"]

    return full_name, discovery_date, source_group


def get_cross_match_info(raw: dict, logger=None):
    """ """
    alert = raw["candidate"]

    label = ""

    # Check if known variable star (https://arxiv.org/pdf/1405.4290.pdf)

    res = ampel_api_catalog(
        catalog="CRTS_DR1",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        searchradius_arcsec=5.0,
        logger=logger,
    )
    if res is not None:
        if logger:
            logger.info(res)
        label = f"[CRTS variable star: {res[0]['body']['name']} ({res[0]['dist_arcsec']:.2f} arsec)]"

    # Check if known QSO/AGN

    res = ampel_api_catalog(
        catalog="milliquas",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        searchradius_arcsec=1.5,
        logger=logger,
    )
    if res is not None:
        if len(res) == 1:

            if "q" in res[0]["body"]["broad_type"]:
                label = f"[MILLIQUAS: {res[0]['body']['name']} - Likely QSO (prob = {res[0]['body']['qso_prob']}%) ({res[0]['dist_arcsec']:.2f} arsec)]"
            else:
                label = f"[MILLIQUAS: {res[0]['body']['name']} - '{res[0]['body']['broad_type']}'-type source ({res[0]['dist_arcsec']:.2f} arsec)]"
        else:
            label = "[MULTIPLE MILLIQUAS MATCHES]"

    # Check if measured parallax in Gaia (i.e galactic)

    if label == "":
        res = ampel_api_catalog(
            catalog="GAIADR2",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=5.0,
            logger=logger,
        )
        if res is not None:
            if res[0]["body"]["Plx"] is not None:
                plx_sig = res[0]["body"]["Plx"] / res[0]["body"]["ErrPlx"]
                if plx_sig > 3.0:
                    label = f"[GAIADR2: {plx_sig:.1f}-sigma parallax ({res[0]['dist_arcsec']:.2f} arsec)]"

    # Check if classified as probable star in SDSS

    if label == "":
        res = ampel_api_catalog(
            catalog="SDSSDR10",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=1.5,
            logger=logger,
        )
        if res is not None:
            if len(res) == 1:
                if float(res[0]["body"]["type"]) == 6.0:
                    label = f"[SDSS Morphology: 'Star'-type source ({res[0]['dist_arcsec']:.2f} arsec)]"
            else:
                label = "[MULTIPLE SDSS MATCHES]"

    # WISE colour cuts (https://iopscience.iop.org/article/10.3847/1538-4365/)

    if label == "":
        res = ampel_api_catalog(
            catalog="wise_color",
            catalog_type="extcats",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=6.0,
            logger=logger,
        )
        if res is not None:
            if len(res) == 1:
                w1mw2 = res[0]["body"]["W1mW2"]
                if w1mw2 > 0.8:
                    label = (
                        f"[Probable WISE-selected quasar:W1-W2={w1mw2:.1f}>0.8  "
                        f"({res[0]['dist_arcsec']:.2f} arsec)]"
                    )
                elif w1mw2 > 0.5:
                    label = (
                        f"[Possible WISE-selected quasar:W1-W2={w1mw2:.1f}>0.5  "
                        f"({res[0]['dist_arcsec']:.2f} arsec)]"
                    )
                else:
                    label = "WISE DETECTIOM"
            else:
                label = "[MULTIPLE WISE MATCHES]"

    # Just check NED

    if label == "":

        res = query_ned_astroquery(
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=3.0,
        )

        if res is not None:
            if len(res) == 1:
                label = (
                    f"{res['Object Name'][0]} ['{res['Type'][0]}'-type source "
                    f"({res['Separation'][0]:.2f} arsec)]"
                )
            elif len(res) > 1:
                label += "[MULTIPLE NED MATCHES]"
                logger.debug(f"{res}")

    # Extra check to TNS, append to other info

    full_name, _, _ = ampel_api_tns(
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
    )

    if full_name is not None:
        label += f" [TNS NAME={full_name}]"

    return label


def check_cross_match_info_by_name(name: str, logger=None):
    """ """
    return get_cross_match_info(
        raw=ampel_api_name(name, with_history=False, logger=logger)[0], logger=logger
    )


def resolve_name(
    source_name: str, source_coords: list = None, source_redshift: float = None
):
    plot_title = source_name

    # If there are no coordinates, try name resolve to get coordinates!

    if source_coords is None:

        # Try ampel to find ZTF coordinates

        if is_ztf_name(name=source_name):
            logger.info("Source name is a ZTF name.")
            res = ampel_api_name(source_name, with_history=False)[0]
            source_coords = [res["candidate"]["ra"], res["candidate"]["dec"]]
            logger.info(f"Found ZTF coordinates for source {source_name}")

        # Try TNS

        elif is_tns_name(name=source_name):
            logger.info("Source name is a TNS name.")
            result_dict = query_tns_by_name(name=source_name, logger=logger)

            if not result_dict:
                logger.warning(f"{source_name} is not in TNS.")

            if result_dict:
                logger.info(f"Found {source_name} on TNS.")
                res = result_dict["data"]["reply"]
                ra = res["radeg"]
                dec = res["decdeg"]
                source_coords = [ra, dec]
                if "redshift" in res.keys():
                    source_redshift = res["redshift"]

        # Otherwise try NED

        else:
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

                if str(result_table["Redshift"][0]) != "--":
                    source_redshift = result_table["Redshift"]

                logger.info(
                    f"Using Astropy NED query result for name {source_name} ({source_coords})"
                )

            if source_coords is None:
                sc = SkyCoord.from_name(source_name)
                logger.info(
                    f"Using Astropy CDS query result for name {source_name} (RA={sc.ra}, Dec={sc.dec})"
                )
                source_coords = (sc.ra.value, sc.dec.value)

    # Try to find a catalogue source nearby using coordinates

    if np.logical_and("ZTF" in source_name, source_coords is not None):

        c = SkyCoord(source_coords[0], source_coords[1], unit=u.deg, frame="icrs")

        r = 0.5 * u.arcsecond

        result_table = Ned.query_region(c, radius=r)

        if len(result_table) == 1:
            if "ZTF" in plot_title:
                plot_title += f' ({result_table["Object Name"][0]})'

            source_coords = [result_table["RA"][0], result_table["DEC"][0]]

            if str(result_table["Redshift"][0]) != "--":
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

    return source_coords, source_redshift, plot_title
