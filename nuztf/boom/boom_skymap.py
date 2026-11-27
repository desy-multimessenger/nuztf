"""
Query a skymap using boom
"""

import logging

import numpy as np
from astropy.time import Time
from blastwave import LSSTClient, ZTFClient
from mocpy import MOC
from tqdm import tqdm

from nuztf.boom.boom_filters import get_quality_filter

logger = logging.getLogger(__name__)

# BOOM rejects a filter test spanning more than this many days
MAX_JD_WINDOW = 7.0

# Program IDs to ask for. BOOM intersects these with the grants the user
# actually has, so requesting all of them means "everything I may see".
PROGRAM_IDS = [1, 2, 3]

# Fields to return. BOOM only exposes the subset of the alert available at
# filter time, and requires the last pipeline stage to be a $project
# including objectId.
PROJECTION = {
    "objectId": 1,
    "candidate": 1,
    "properties": 1,
    "classifications": 1,
}


def ids_to_moc_ascii(nside: int, ids: list[int]) -> str:
    """
    Convert nested healpix indices to a MOC in IVOA ASCII form

    :param nside: nside of the healpix indices
    :param ids: nested healpix indices
    :return: MOC string, e.g. "4/278 5/1116 1118 6/4427-4428"
    """

    depth = int(np.log2(nside))

    if 2**depth != nside:
        raise ValueError(f"nside must be a power of two, not {nside}")

    moc = MOC.from_healpix_cells(
        ipix=np.array(ids),
        depth=np.full(len(ids), depth),
        max_depth=depth,
    )

    return moc.to_string(format="ascii")


def split_jd_window(t_min: Time, t_max: Time) -> list[tuple[float, float]]:
    """
    Split a time range into equal chunks that BOOM will accept

    :param t_min: minimum time
    :param t_max: maximum time
    :return: list of (start_jd, end_jd) pairs, each no longer than the cap
    """

    n_chunks = max(1, int(np.ceil((t_max.jd - t_min.jd) / MAX_JD_WINDOW)))

    # A window that is an exact multiple of the cap can land a hair over it
    # once the edges are rounded to float64, which BOOM would reject
    while True:
        edges = np.linspace(t_min.jd, t_max.jd, n_chunks + 1)
        if n_chunks == 1 or np.all(np.diff(edges) <= MAX_JD_WINDOW):
            break
        n_chunks += 1

    return [(float(x), float(y)) for x, y in zip(edges[:-1], edges[1:])]


def boom_api_skymap(
    cone_nside: int,
    cone_ids: list[int],
    t_min: Time,
    t_max: Time,
    client: ZTFClient | LSSTClient | None = None,
) -> list[dict]:
    """
    Query BOOM for objects in a skymap

    The healpix indices are sent as a single MOC, so the entire region is
    matched by healpix range in one query against the
    {coordinates.hpx, candidate.jd} index, rather than one cone per pixel.
    BOOM caps a filter test at 7 days, so longer windows are split.

    Unlike a cone search, the region is matched exactly, so alerts just
    outside the requested pixels are no longer returned.

    :param cone_nside: nside of the skymap
    :param cone_ids: list of cone ids
    :param t_min: minimum time
    :param t_max: maximum time
    :param client: BOOM client
    :return: list of alerts
    """

    if client is None:
        client = ZTFClient()

    if len(cone_ids) == 0:
        logger.warning("No healpix indices to query, skipping")
        return []

    moc_ascii = ids_to_moc_ascii(nside=cone_nside, ids=cone_ids)

    quality_cuts = get_quality_filter(catalog=client.catalog, jd=t_max.jd)

    new_cuts = {
        "properties.rock": {"$eq": False},
        "properties.star": {"$eq": False},
        "properties.near_brightstar": {"$eq": False},
        "properties.stationary": {"$eq": True},
    }

    # BOOM prepends the region and the time cut, and appends nothing, so the
    # pipeline is just our own cuts plus the projection it requires
    pipeline = [
        {"$match": {**quality_cuts, **new_cuts}},
        {"$project": PROJECTION},
    ]

    chunks = split_jd_window(t_min=t_min, t_max=t_max)

    logger.info(f"Querying {len(cone_ids)} pixels in {len(chunks)} time chunk(s)")

    res: dict[int, dict] = {}

    for start_jd, end_jd in tqdm(chunks):
        response = client.api(
            "post",
            "filters/test",
            data={
                "survey": client.base_catalog,
                "permissions": {client.base_catalog: PROGRAM_IDS},
                "start_jd": start_jd,
                "end_jd": end_jd,
                "moc_ascii": moc_ascii,
                "pipeline": pipeline,
            },
        )

        if not response.ok:
            logger.error(f"BOOM rejected the skymap query: {response.text}")

        response.raise_for_status()

        # chunk edges are inclusive, so an alert can be returned twice
        for alert in response.json()["data"]["results"]:
            res[alert["_id"]] = alert

    alerts = list(res.values())

    # propagate candid
    for alert in alerts:
        alert["candid"] = alert["_id"]
        alert["prv_candidates"] = []

    return alerts
