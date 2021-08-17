import requests
import backoff
from base64 import b64decode
from astropy.time import Time
from nuztf.credentials import load_credentials
from requests.auth import HTTPBasicAuth

# AMPEL API URLs

API_BASEURL = "https://ampel.zeuthen.desy.de"
API_ZTF_ARCHIVE_URL = API_BASEURL + "/api/ztf/archive"
API_CATALOGMATCH_URL = API_BASEURL + "/api/catalogmatch"
API_CUTOUT_URL = API_BASEURL + "/api/ztf/archive/cutouts"

api_user, api_pass = load_credentials("ampel_api")


def merge_alerts(alert_list):
    merged_list = []
    keys = list(set([x["objectId"] for x in alert_list]))

    for objectid in keys:
        alerts = [x for x in alert_list if x["objectId"] == objectid]
        if len(alerts) == 1:
            merged_list.append(alerts[0])
        else:
            jds = [x["candidate"]["jd"] for x in alerts]
            order = [jds.index(x) for x in sorted(jds)[::-1]]
            latest = alerts[jds.index(max(jds))]
            latest["candidate"]["jdstarthist"] = min(
                [x["candidate"]["jdstarthist"] for x in alerts]
            )

            for index in order[1:]:

                x = alerts[index]

                # Merge previous detections

                for prv in x["prv_candidates"] + [x["candidate"]]:
                    if prv not in latest["prv_candidates"]:
                        latest["prv_candidates"] = [prv] + latest["prv_candidates"]

            merged_list.append(latest)
    return merged_list

@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)


def ampel_api_cone(ra, dec, radius,
                   t_min_jd=Time(
                       '2018-04-01T00:00:00.123456789',
                       format='isot',
                       scale='utc'
                   ).jd,
                   t_max_jd=Time.now().jd, logger=None):
    """Function to query ampel via a cone search"""

    queryurl_conesearch = (
            API_ZTF_ARCHIVE_URL
            + f"/alerts/cone_search?ra={ra}&dec={dec}&"
              f"radius={radius}&jd_start={t_min_jd}&"
              f"jd_end={t_max_jd}&with_history=false&"
              f"with_cutouts=false&chunk_size=500"
    )

    if logger is not None:
        logger.debug(queryurl_conesearch)

    response = requests.get(
        queryurl_conesearch,
        auth=HTTPBasicAuth(api_user, api_pass),
    )
    if response.status_code == 503:
        raise requests.exceptions.RequestException

    query_res = [i for i in response.json()["alerts"]]

    return query_res

@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)


def ampel_api_name(ztf_name, logger=None):
    """Function to query ampel via name"""

    queryurl_ztf_name = (
        API_ZTF_ARCHIVE_URL + f"/object/{ztf_name}/alerts?with_history=true"
    )

    if logger is not None:
        logger.debug(queryurl_ztf_name)

    response = requests.get(
        queryurl_ztf_name,
        auth=HTTPBasicAuth(api_user, api_pass),
    )
    if response.status_code == 503:
        raise requests.exceptions.RequestException
    query_res = [i for i in response.json()]
    query_res = merge_alerts(query_res)
    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)

def ampel_api_cutout(candid: int, logger=None):
    """Function to query ampel for cutouts by candidate ID"""
    queryurl_cutouts = API_CUTOUT_URL + f"/{candid}"
    response = requests.get(
        queryurl_cutouts,
        auth=HTTPBasicAuth(api_user, api_pass),
    )
    if logger is not None:
        self.logger.debug(queryurl_cutouts)

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    cutouts = response.json()
    return cutouts


def reassemble_alert(mock_alert):
    """Function to recreate ztf alerts"""
    cutouts = ampel_api_cutout(mock_alert["candid"])

    for k in cutouts:
        mock_alert[f"cutout{k.title()}"] = {
            "stampData": b64decode(cutouts[k]),
            "fileName": "dunno",
        }

    mock_alert["schemavsn"] = "dunno"
    mock_alert["publisher"] = "dunno"
    for pp in [mock_alert["candidate"]] + mock_alert["prv_candidates"]:
        pp["pdiffimfilename"] = "dunno"
        pp["programpi"] = "dunno"
        pp["ssnamenr"] = "dunno"

    return mock_alert


def ampel_api_tns(ra: float, dec: float, searchradius_arcsec: float = 3):
    """Function to query TNS via ampel api"""
    queryurl_catalogmatch = API_CATALOGMATCH_URL + f"/cone_search/nearest"

    # First, we create a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    query = {
        "ra_deg": ra,
        "dec_deg": dec,
        "catalogs": [
            {"name": "TNS", "rs_arcsec": searchradius_arcsec, "use": "extcats"}
        ],
    }

    # Now we retrieve results from the API
    response = requests.post(url=queryurl_catalogmatch, json=query, headers=headers)

    full_name = None
    discovery_date = None
    source_group = None

    if response.json()[0]:
        response_body = response.json()[0]["body"]
        print(response_body)
        name = response_body["objname"]
        prefix = response_body["name_prefix"]
        full_name = prefix + name
        discovery_date = response_body["discoverydate"]
        if "source_group" in response_body.keys():
            source_group = response_body["source_group"]["group_name"]

    return full_name, discovery_date, source_group


def ampel_api_catalog(
        catalog: str,
        catalog_type: str,
        ra: float,
        dec: float,
        searchradius_arcsec: float = 10,
        searchtype: str = "all",
    ):
    """
    Method for querying catalogs via the Ampel API
    'catalog' must be the name of a supported catalog, e.g.
    SDSS_spec, PS1, NEDz_extcats...
    For a full list of catalogs, confer
    https://ampel.zeuthen.desy.de/api/catalogmatch/catalogs

    """
    assert catalog_type in ["extcats", "catsHTM"]
    assert searchtype in ["all", "nearest"]

    queryurl_catalogmatch = API_CATALOGMATCH_URL + "/cone_search/" + searchtype

    # First, we create a json body to post
    headers = {"accept": "application/json", "Content-Type": "application/json"}
    query = {
        "ra_deg": ra,
        "dec_deg": dec,
        "catalogs": [
            {"name": catalog, "rs_arcsec": searchradius_arcsec, "use": catalog_type}
        ],
    }

    response = requests.post(
        url=queryurl_catalogmatch, json=query, headers=headers
    ).json()[0]

    return response
