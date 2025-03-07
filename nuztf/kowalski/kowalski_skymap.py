"""
This module contains functions to query Kowalski for objects in a skymap
"""

import healpy as hp
from astropy.time import Time
from penquins import Kowalski
from pydantic import BaseModel

from nuztf.kowalski.config import get_kowalski


class Cone(BaseModel):
    """
    Class to represent a cone in the sky
    """

    ra: float
    dec: float
    radius: float

    def __repr__(self):
        return f"Cone(ra={self.ra}, dec={self.dec}, radius={self.radius})"

    def __str__(self):
        return self.__repr__()


def get_cones_for_map(nside: int, cone_ids: list[id]) -> list[Cone]:
    """
    Function to get cones from a skymap file

    :param path: path to skymap file
    :param prob_threshold: cumulative probability threshold
    :return: list of cones
    """
    cones = []

    for value in cone_ids:
        ra, dec = hp.pix2ang(nside, value, lonlat=True, nest=True)
        r = hp.max_pixrad(nside, degrees=True)
        cones.append(Cone(ra=ra, dec=dec, radius=r))
    return cones


def kowalski_api_skymap(
    cone_nside: int,
    cone_ids: list[int],
    t_min: Time,
    t_max: Time,
    kowalski: Kowalski | None = None,
    max_n_threads: int = 10,
) -> list[dict]:
    """
    Query Kowalski for objects in a skymap

    :param cone_nside: nside of the skymap
    :param cone_ids: list of cone ids
    :param t_min: minimum time
    :param t_max: maximum time
    :param kowalski: Kowalski object
    :param max_n_threads: maximum number of threads
    """

    if kowalski is None:
        kowalski = get_kowalski()

    cones = get_cones_for_map(nside=cone_nside, cone_ids=cone_ids)

    filter_dict = {
        "candidate.jd": {"$gt": t_min.jd, "$lt": t_max.jd},
        # "candidate.jdstarthist": {"$gt": t_min.jd, "$lt": t_max.jd},
        # "candidate.jdendhist": {"$gt": t_min.jd, "$lt": t_max.jd},
        # "candidate.isdiffpos": {"$in": ["1", "t", "true", "True", "T", 1]},
        # "candidate.drb": {"$gt": 0.3},
        # "candidate.magpsf": {"$gt": 15},
        # "candidate.ndethist": {"$gt": 0, "$lte": max_n_detections},
    }

    queries = []
    for cone in cones:
        query = {
            "query_type": "cone_search",
            "filter": filter_dict,
            "query": {
                "object_coordinates": {
                    "cone_search_radius": cone.radius,
                    "cone_search_unit": "deg",
                    "radec": {"object": [cone.ra, cone.dec]},
                },
                "catalogs": {
                    "ZTF_alerts": {
                        "filter": filter_dict,
                        "projection": {
                            "_id": 0,
                            "cutoutScience": 0,
                            "cutoutTemplate": 0,
                            "cutoutDifference": 0,
                            "coordinates": 0,
                            "prv_candidates": 0,
                        },
                    },
                },
            },
        }

        queries.append(query)

    response = kowalski.query(
        queries=queries, use_batch_query=True, max_n_threads=max_n_threads
    )

    results = []
    candids = []

    # first we have one response per query. Each response contains a dict with one key per instance
    for name in list(response.keys()):
        for response in response[name]:
            data = response.get("data", None)
            if data is None:
                continue

            matches = data["ZTF_alerts"]["object"]
            if len(matches) > 0:
                for match in matches:
                    if match["candid"] in candids:
                        continue
                    match["prv_candidates"] = []
                    results.append(match)
                    candids.append(match["candid"])
    return results
