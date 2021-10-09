import re, logging

import requests
import numpy as np
from astropy.time import Time

base_gcn_url = "https://gcn.gsfc.nasa.gov/gcn3"


def gcn_url(gcn_number):
    """ """
    return f"{base_gcn_url}/{gcn_number}.gcn3"


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


def parse_gcn_archive():
    """ """
    page = requests.get(f"{base_gcn_url}_archive.html")

    nu_circulars = []

    for line in page.text.splitlines():
        if "IceCube observation of a high-energy neutrino" in line:
            res = line.split(">")
            gcn_no = "".join([x for x in res[2] if x.isdigit()])
            name = res[3].split(" - ")[0]
            nu_circulars.append((name, gcn_no))

    return nu_circulars


def parse_gcn_for_no(
    base_nu_name: str, url: str = f"{base_gcn_url}_archive.html", logger=None
):
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


def find_gcn_no(base_nu_name: str, logger=None):
    """ """
    gcn_no, name, latest_archive_no = parse_gcn_for_no(base_nu_name, logger=logger)

    if gcn_no is None:
        logging.info(
            f"No GCN found for {base_nu_name} on GCN page, checking archive instead. "
            f"The latest page is {latest_archive_no}"
        )

        while np.logical_and(latest_archive_no > 0, gcn_no is None):
            gcn_no, name, _ = parse_gcn_for_no(
                base_nu_name,
                url=f"{base_gcn_url}_arch_old{latest_archive_no}.html",
            )
            latest_archive_no -= 1

    # while

    if name is None:
        raise ParsingError("No GCN match found for {0}".format(base_nu_name))

    logging.info(f"Match is {name} (GCN #{gcn_no})")

    return gcn_no


def get_latest_gcn(logger=None):
    """ """
    latest = parse_gcn_archive()[0]
    if logger:
        logger.info(f"Latest GCN is {latest[0]} (GCN #{latest[1]})")
    return latest[1]


def parse_radec(str: str):
    """ """
    regex_findall = re.findall(r"[-+]?\d*\.\d+|\d+", str)
    if len(regex_findall) == 4:
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
            ra_err = [ra_upper, -ra_lower]
            dec_err = [dec_upper, -dec_lower]
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
