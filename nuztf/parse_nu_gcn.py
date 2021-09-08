import requests
import numpy as np
from ztf_plan_obs.gcn_parser import parse_gcn_circular


base_gcn_url = "https://gcn.gsfc.nasa.gov/gcn3"


def gcn_url(gcn_number):
    return f"{base_gcn_url}/{gcn_number}.gcn3"


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


def parse_gcn_archive():
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
        base_nu_name, url=f"{base_gcn_url}_archive.html"
):
    print(f"Checking for GCN on {url}")

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


def find_gcn_no(base_nu_name):
    gcn_no, name, latest_archive_no = parse_gcn_for_no(base_nu_name)

    if gcn_no is None:

        print(
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

    print(f"Match is {name} (GCN #{gcn_no})")

    return gcn_no


def get_latest_gcn():
    latest = parse_gcn_archive()[0]
    print(f"Latest GCN is {latest[0]} (GCN #{latest[1]})")
    return latest[1]
