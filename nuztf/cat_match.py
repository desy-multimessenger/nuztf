from nuztf.ampel_api import ampel_api_catalog


def query_ned_for_z(ra_deg: float, dec_deg: float, searchradius_arcsec: float = 20):

    z = None
    dist_arcsec = None

    query = ampel_api_catalog(
        catalog="NEDz_extcats",
        catalog_type="extcats",
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        searchradius_arcsec=searchradius_arcsec,
        searchtype="nearest",
    )

    if query:
        z = query["body"]["z"]
        dist_arcsec = query["dist_arcsec"]
    return z, dist_arcsec


def ampel_api_tns(ra_deg: float, dec_deg: float, searchradius_arcsec: float = 3):
    """Function to query TNS via ampel api"""

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


search_rad = 1.5


def get_cross_match_info(raw):
    alert = raw["candidate"]

    label = ""

    # Check if known QSO/AGN

    res = ampel_api_catalog(
        catalog="milliquas",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        searchradius_arcsec=search_rad
    )
    if res is not None:
        if len(res) == 1:

            if "q" in res[0]['body']['broad_type']:
                label = f"[MILLIQUAS: Likely QSO (prob = {res[0]['body']['qso_prob']}%) ({res[0]['dist_arcsec']:.2f} arsec)]"
            else:

                label = f"[MILLIQUAS: '{res[0]['body']['broad_type']}'-type source ({res[0]['dist_arcsec']:.2f} arsec)]"
        else:
            label = "[MULTIPLE MILLIQUAS MATCHES]"

    # Check if measured parallax in Gaia (i.e galactic)

    if label == "":
        res = ampel_api_catalog(
            catalog="GAIADR2",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=5.
        )
        if res is not None:
            if res[0]['body']['Plx'] is not None:
                plx_sig = res[0]['body']['Plx']/res[0]['body']['ErrPlx']
                if plx_sig > 1.5:
                    label = f"[GAIADR2: {plx_sig:.1f}-sigma parallax ({res[0]['dist_arcsec']:.2f} arsec)]"

    # Check if classified as probable star in SDSS

    if label == "":
        res = ampel_api_catalog(
            catalog="SDSSDR10",
            catalog_type="catsHTM",
            ra_deg=alert["ra"],
            dec_deg=alert["dec"],
            searchradius_arcsec=search_rad
        )
        if res is not None:
            if len(res) == 1:
                if float(res[0]['body']['type']) == 6.0:
                    label = f"[SDSS Morphology: 'Star'-type source ({res[0]['dist_arcsec']:.2f} arsec)]"
            else:
                label = "[MULTIPLE SDSS MATCHES]"

    return label
