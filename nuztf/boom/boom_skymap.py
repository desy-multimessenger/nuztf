"""
Query a skymap using boom
"""

import logging

from astropy.time import Time
from blastwave import LSSTClient, ZTFClient
from blastwave.models.query import BOOMQuery
from tqdm import tqdm

from nuztf.kowalski.kowalski_skymap import get_cones_for_map

logger = logging.getLogger(__name__)


def boom_api_skymap(
    cone_nside: int,
    cone_ids: list[int],
    t_min: Time,
    t_max: Time,
    client: ZTFClient | LSSTClient | None = None,
) -> list[dict]:
    """
    Query BOOM for objects in a skymap

    :param cone_nside: nside of the skymap
    :param cone_ids: list of cone ids
    :param t_min: minimum time
    :param t_max: maximum time
    :param client: BOOM client
    """

    if client is None:
        client = LSSTClient()

    cones = get_cones_for_map(nside=cone_nside, cone_ids=cone_ids)

    time_cut = {
        "candidate.jd": {"$gt": t_min.jd, "$lt": t_max.jd},
    }

    rb_cut = {
        "candidate.reliability": {"$gt": 0.4},
        "candidate.isDipole": {"$eq": False},
        "candidate.pixelFlags": {"$eq": False},
        "candidate.glint_trail": {"$eq": False},
        "candidate.isNegative": {"$eq": False},
        "candidate.centroid_flag": {"$eq": False},
        "candidate.psfFlux_flag": {"$eq": False},
    }

    new_cuts = {
        "properties.rock": {"$eq": False},
        "properties.star": {"$eq": False},
        "properties.near_brightstar": {"$eq": False},
        "properties.stationary": {"$eq": True},
    }

    filter_dict: dict = {**time_cut, **rb_cut, **new_cuts}

    queries = []
    for cone in cones:
        query = BOOMQuery(
            catalog_name=client.catalog,
            filter={
                "coordinates.radec_geojson": {
                    "$nearSphere": {
                        "$geometry": {
                            "type": "Point",
                            "coordinates": [cone.ra - 180.0, cone.dec],
                        },
                        "$maxDistance": client.get_near_sphere_dist(
                            3600.0 * cone.radius
                        ),
                    }
                },
                **filter_dict,
            },
            projection={
                "objectId": 1,
                "candidate.ra": 1,
                "candidate.dec": 1,
                "candidate.jd": 1,
            },
        )
        queries.append(query)

    logger.info(f"Found {len(queries)} cones")

    res = []

    for query in tqdm(queries):
        res += client.query(query)
    return res
