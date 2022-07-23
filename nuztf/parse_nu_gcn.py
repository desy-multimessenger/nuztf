#!/usr/bin/env python3
# License: BSD-3-Clause

import re, json, logging

import requests
import numpy as np
from astropy.time import Time

BASE_GCN_URL = "https://gcn.gsfc.nasa.gov/gcn3"


def gcn_url(gcn_number):
    """ """
    return f"{BASE_GCN_URL}/{gcn_number}.gcn3"


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


def parse_gcn_for_no(base_nu_name: str, url: str = f"{BASE_GCN_URL}_archive.html"):
    """ """
    if logger:
        logger.info(f"Checking for GCN on {url}")

    nu_name = str(base_nu_name)

    page = requests.get(url)

    gcn_no = None
    name = None

    while not nu_name[0].isdigit():
        nu_name = nu_name[1:]

    latest_archive_no = None

    for line in page.text.splitlines():
        if np.logical_and(
            "IceCube observation of a high-energy neutrino" in line, nu_name in line
        ):
            res = line.split(">")
            if gcn_no is None:
                gcn_no = "".join([x for x in res[2] if x.isdigit()])
                name = res[3].split(" - ")[0]
                print(f"Found match to {base_nu_name}: {name}")
            else:
                raise Exception(f"Multiple matches found to {base_nu_name}")

        elif np.logical_and("gcn3_arch_old" in line, latest_archive_no is None):
            url = line.split('"')[1]
            latest_archive_no = int(url[13:].split(".")[0])

    return gcn_no, name, latest_archive_no


def find_gcn_no(base_nu_name: str):
    """
    Trick the webpage into giving us results
    """
    endpoint = (
        "https://heasarc.gsfc.nasa.gov/wsgi-scripts/tach/gcn_v2/tach.wsgi/graphql_fast"
    )

    querystr = (
        '{ allEventCard( name: "'
        + base_nu_name
        + '" ) {edges {node { id_ event } } } }'
    )
    r = requests.post(
        endpoint,
        data={
            "query": querystr,
            "Content-Type": "application/json",
        },
    )
    res = json.loads(r.text)

    if res["data"]["allEventCard"]["edges"]:

        event_id = res["data"]["allEventCard"]["edges"][0]["node"]["id_"]

        querystr = (
            "{ allCirculars ( evtid:"
            + event_id
            + " ) { totalCount edges { node { id id_ received subject evtidCircular{ event } cid evtid oidCircular{ telescope detector oidEvent{ wavelength messenger } } } } } }"
        )

        r = requests.post(
            endpoint,
            data={
                "query": querystr,
                "Content-Type": "application/json",
            },
        )
        result = json.loads(r.text)

        received_date = []
        circular_nr = []

        for entry in result["data"]["allCirculars"]["edges"]:
            received_date.append(entry["node"]["received"])
            circular_nr.append(entry["node"]["cid"])

        """
        I don't trust this webserver, let's go with the
        earliest GCN, not the last in the list
        """
        gcn_no = circular_nr[
            np.argmin([Time(i, format="isot").mjd for i in received_date])
        ]
        logging.info(f"GCN found ({gcn_no})")

        return gcn_no

    else:
        logging.warning(f"No GCN found for {base_nu_name}")

        return None


def get_latest_gcn():
    """
    Get the last circular
    """
    endpoint = (
        "https://heasarc.gsfc.nasa.gov/wsgi-scripts/tach/gcn_v2/tach.wsgi/graphql_fast"
    )
    querystr = '{ allCirculars ( first:50after:"" ) { totalCount pageInfo{ hasNextPage hasPreviousPage startCursor endCursor } edges { node { id id_ received subject evtidCircular{ event } cid evtid oidCircular{ telescope detector oidEvent{ wavelength messenger } } } } } }'

    r = requests.post(
        endpoint,
        data={
            "query": querystr,
            "Content-Type": "application/json",
        },
    )
    result = json.loads(r.text)

    received_date = []
    circular_nr = []

    for entry in result["data"]["allCirculars"]["edges"]:
        received_date.append(entry["node"]["received"])
        circular_nr.append(entry["node"]["cid"])

    latest_gcn_no = circular_nr[
        np.argmin([Time(i, format="isot").mjd for i in received_date])
    ]

    logging.info(f"Most recent GCN is {latest_gcn_no}")

    return latest_gcn_no


def parse_radec(string: str):
    """ """
    regex_findall = re.findall(r"[-+]?\d*\.\d+|\d+", string)

    if len(regex_findall) == 2:
        pos = float(regex_findall[0])
        pos_upper = None
        pos_lower = None
    elif len(regex_findall) == 4:
        pos = float(regex_findall[0])
        pos_upper = float(regex_findall[1])
        pos_lower = float(regex_findall[1])
    elif len(regex_findall) == 5:
        pos, pos_upper, pos_lower = regex_findall[0:3]
        pos = float(pos)
        pos_upper = float(pos_upper.replace("+", ""))
        pos_lower = float(pos_lower.replace("-", ""))
    else:
        raise ParsingError(f"Could not parse GCN ra and dec")

    return pos, pos_upper, pos_lower


def parse_gcn_circular(gcn_number: int):
    """ """
    url = f"https://gcn.gsfc.nasa.gov/gcn3/{gcn_number}.gcn3"
    response = requests.get(url)
    returndict = {}
    mainbody_starts_here = 999
    splittext = response.text.splitlines()
    splittext = list(filter(None, splittext))
    for i, line in enumerate(splittext):
        if "SUBJECT" in line:
            name = line.split(" - ")[0].split(": ")[1]
            returndict.update({"name": name})
        elif "FROM" in line:
            base = line.split("at")[0].split(": ")[1].split(" ")
            author = [x for x in base if x != ""][1]
            returndict.update({"author": author})
        elif (
            ("RA" in line or "Ra" in line)
            and ("DEC" in splittext[i + 1] or "Dec" in splittext[i + 1])
            and i < mainbody_starts_here
        ):
            ra, ra_upper, ra_lower = parse_radec(line)
            dec, dec_upper, dec_lower = parse_radec(splittext[i + 1])
            if ra_upper and ra_lower:
                ra_err = [ra_upper, -ra_lower]
            else:
                ra_err = [None, None]

            if dec_upper and dec_lower:
                dec_err = [dec_upper, -dec_lower]
            else:
                dec_err = [None, None]
            returndict.update(
                {"ra": ra, "ra_err": ra_err, "dec": dec, "dec_err": dec_err}
            )
            mainbody_starts_here = i + 2
        elif ("Time" in line or "TIME" in line) and i < mainbody_starts_here:
            raw_time = [
                x for x in line.split(" ") if x not in ["Time", "", "UT", "UTC"]
            ][1]
            raw_time = "".join(
                [x for x in raw_time if np.logical_or(x.isdigit(), x in [":", "."])]
            )
            raw_date = name.split("-")[1][:6]
            ut_time = f"20{raw_date[0:2]}-{raw_date[2:4]}-{raw_date[4:6]}T{raw_time}"
            time = Time(ut_time, format="isot", scale="utc")
            returndict.update({"time": time})

    return returndict
