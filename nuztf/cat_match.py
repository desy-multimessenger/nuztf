#!/usr/bin/env python
# coding: utf-8

import json
import logging
import warnings

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyWarning
from astroquery.exceptions import RemoteServiceError
from astroquery.ipac.irsa import Irsa
from astroquery.ipac.ned import Ned

from nuztf.ampel_api import ampel_api_catalog, ampel_api_name
from nuztf.paths import CROSSMATCH_CACHE


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
        search_radius_arcsec=searchradius_arcsec,
        search_type="nearest",
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


def query_wise_astroquery(
    ra_deg: float, dec_deg: float, searchradius_arcsec: float = 3.0
):
    """
    Function to obtain WISE crossmatches via astroquery

    :param ra_deg: Right ascension (deg)
    :param dec_deg: Declination (deg)
    :param searchradius_arcsec: Search radius (arcsec)
    :return: result of query
    """
    c = SkyCoord(ra_deg, dec_deg, unit=u.deg, frame="icrs")

    r = searchradius_arcsec * u.arcsecond

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        allwise = Irsa.query_region(c, catalog="allwise_p3as_psd", radius=r)
    return allwise


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
        search_radius_arcsec=searchradius_arcsec,
        search_type="nearest",
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

    cache_file = CROSSMATCH_CACHE.joinpath(f"{raw['objectId']}.json")

    if cache_file.exists():
        with open(cache_file) as f:
            res = json.load(f)
            label = res["data"]
            return label

    alert = raw["candidate"]

    label = ""

    if logger is None:
        logger = logging.getLogger()

    # Check if known variable star (https://arxiv.org/pdf/1405.4290.pdf)

    res = ampel_api_catalog(
        catalog="CRTS_DR1",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        search_radius_arcsec=5.0,
        logger=logger,
    )
    if res is not None:
        if logger:
            logger.info(res)
        label = (
            f"[CRTS variable star: "
            f"{res[0]['body']['name']} ({res[0]['dist_arcsec']:.2f} arsec)]"
        )

    # Check if known QSO/AGN

    res = ampel_api_catalog(
        catalog="milliquas",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        search_radius_arcsec=1.5,
        logger=logger,
    )
    if res is not None:
        if len(res) == 1:
            if "q" in res[0]["body"]["broad_type"]:
                label = (
                    f"[MILLIQUAS: {res[0]['body']['name']} - "
                    f"Likely QSO (prob = {res[0]['body']['qso_prob']}%) "
                    f"({res[0]['dist_arcsec']:.2f} arsec)]"
                )
            else:
                label = (
                    f"[MILLIQUAS: {res[0]['body']['name']} - "
                    f"'{res[0]['body']['broad_type']}'-type source "
                    f"({res[0]['dist_arcsec']:.2f} arsec)]"
                )
        else:
            label = "[MULTIPLE MILLIQUAS MATCHES]"

    # Check if measured parallax in Gaia (i.e galactic)

    if label == "":
        res = ampel_api_catalog(
            catalog="GAIADR2",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            search_radius_arcsec=5.0,
            logger=logger,
        )
        if res is not None:
            if res[0]["body"]["Plx"] is not None:
                plx_sig = res[0]["body"]["Plx"] / res[0]["body"]["ErrPlx"]
                if plx_sig > 3.0:
                    label = (
                        f"[GAIADR2: {plx_sig:.1f}-sigma parallax "
                        f"({res[0]['dist_arcsec']:.2f} arsec)]"
                    )

    # Check if classified as probable star in SDSS

    if label == "":
        res = ampel_api_catalog(
            catalog="SDSSDR10",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            search_radius_arcsec=1.5,
            logger=logger,
        )
        if res is not None:
            if len(res) == 1:
                if float(res[0]["body"]["type"]) == 6.0:
                    label = (
                        f"[SDSS Morphology: 'Star'-type source "
                        f"({res[0]['dist_arcsec']:.2f} arsec)]"
                    )
            else:
                label = "[MULTIPLE SDSS MATCHES]"

    # WISE colour cuts (https://iopscience.iop.org/article/10.3847/1538-4365/)

    if label == "":
        res = query_wise_astroquery(
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=3.0,
        )
        if res is not None:
            if len(res) > 0:
                w1mw2 = res["w1mpro"][0] - res["w2mpro"][0]

                if w1mw2 > 0.8:
                    label = (
                        f"[Probable WISE-selected quasar:W1-W2={w1mw2:.2f}>0.8  "
                        f"({res[0]['dist']:.2f} arsec)]"
                    )
                elif w1mw2 > 0.5:
                    label = (
                        f"[Possible WISE-selected quasar:W1-W2={w1mw2:.2f}>0.5  "
                        f"({res[0]['dist']:.2f} arsec)]"
                    )
                else:
                    label = (
                        f"WISE DETECTION: W1-W2={w1mw2:.2f} "
                        f"({res[0]['dist']:.2f} arsec)"
                    )

                if len(res) > 1:
                    label += "[MULTIPLE WISE MATCHES]"

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

    with open(cache_file, "w") as f:
        json.dump({"data": label}, f)

    return label


def check_cross_match_info_by_name(name: str, logger=None):
    """
    Utility function to check cross-match info for a given name

    :param name: ZTF name
    :param logger:
    :return:
    """
    return get_cross_match_info(
        raw=ampel_api_name(name, with_history=False, logger=logger)[0], logger=logger
    )
