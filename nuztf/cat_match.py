#!/usr/bin/env python3

from astropy.coordinates import SkyCoord
from astropy import units as u
from nuztf.ampel_api import ampel_api_catalog, ampel_api_name


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
            searchradius_arcsec=1.5,
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
                elif w1mw2 > 0.8:
                    label = (
                        f"[Possible WISE-selected quasar:W1-W2={w1mw2:.1f}>0.5  "
                        f"({res[0]['dist_arcsec']:.2f} arsec)]"
                    )
            else:
                label = "[MULTIPLE WISE MATCHES]"

    return label


def check_cross_match_info_by_name(name: str, logger=None):
    """ """
    return get_cross_match_info(
        raw=ampel_api_name(name, with_history=False, logger=logger)[0], logger=logger
    )
