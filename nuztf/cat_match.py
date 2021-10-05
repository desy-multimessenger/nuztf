import backoff
import requests
from json import JSONDecodeError
from nuztf.ampel_api import ampel_api_catalog, API_CATALOGMATCH_URL


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


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_tns(ra_deg: float, dec_deg: float, searchradius_arcsec: float = 3):
    """Function to query TNS via ampel api"""
    queryurl_catalogmatch = API_CATALOGMATCH_URL + f"/cone_search/nearest"

    # First, we create a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    query = {
        "ra_deg": ra_deg,
        "dec_deg": dec_deg,
        "catalogs": [
            {"name": "TNS", "rs_arcsec": searchradius_arcsec, "use": "extcats"}
        ],
    }

    # Now we retrieve results from the API
    response = requests.post(url=queryurl_catalogmatch, json=query, headers=headers)

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    full_name = None
    discovery_date = None
    source_group = None

    try:
        res = response.json()
    except JSONDecodeError:
        print(response)
        raise Exception

    if res[0]:
        response_body = res[0]["body"]
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

    res = ampel_api_catalog(
        catalog="milliquas",
        catalog_type="extcats",
        ra_deg=alert["ra"],
        dec_deg=alert["dec"],
        searchradius_arcsec=search_rad
    )
    if res is not None:
        if len(res) == 1:
            label = f"MILLIQUAS: '{res[0]['body']['broad_type']}'-type source ({res[0]['dist_arcsec']:.2f} arsec)"
        else:
            label = "MULTIPLE MILLIQUAS MATCHES"

    return label
