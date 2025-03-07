#!/usr/bin/env python3


import gzip
import logging
import time
from json import JSONDecodeError

import backoff
import requests
from ampel.util.json import load
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from astropy.io import fits  # type: ignore
from astropy.time import Time  # type: ignore

from nuztf import utils
from nuztf.ampel.urls import API_ZTF_ARCHIVE_URL
from nuztf.ampel.utils import get_ampel_token
from nuztf.credentials import load_credentials
from nuztf.paths import PREPROCESSED_CACHE_DIR

MAX_N_PIX = 1000

logger = logging.getLogger(__name__)


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=600,
)
def ampel_api_healpix(
    ipix: int,
    nside: int = 64,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    logger=None,
) -> list:
    """Function to query ampel based on a healpix pixel-index (nside is the pixelization degree)"""

    if logger is None:
        logger = logging.getLogger(__name__)

    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    queryurl_healpix = (
        API_ZTF_ARCHIVE_URL
        + f"/alerts/healpix?nside={nside}&ipix={ipix}&jd_start={t_min_jd}&jd_end={t_max_jd}&with_history={hist}"
        f"&with_cutouts={cutouts}&chunk_size={chunk_size}"
    )

    logger.debug(queryurl_healpix)

    headers = {"Authorization": f"Bearer {get_ampel_token()}"}

    response = requests.get(
        queryurl_healpix,
        headers=headers,
    )

    if response.status_code == 503:
        raise requests.exceptions.RequestException

    try:
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=1200,
)
def ampel_api_acknowledge_chunk(resume_token: str, chunk_id: int, logger=None):
    """
    After receiving a chunk, acknowledge that we got it
    (otherwise large alert queries will start looping)
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    endpoint = (
        API_ZTF_ARCHIVE_URL + f"/stream/{resume_token}/chunk/{chunk_id}/acknowledge"
    )

    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {get_ampel_token()}",
    }

    payload = {"resume_token": resume_token, "chunk_id": chunk_id}

    logger.debug(f"Acknowledging:\n{payload}")

    response = requests.post(url=endpoint, json=payload, headers=headers)


@backoff.on_exception(
    backoff.expo,
    requests.exceptions.RequestException,
    max_time=1200,
)
def ampel_api_skymap_single(
    pixels: list,
    nside: int = 64,
    t_min_jd=Time("2018-04-01T00:00:00.123456789", format="isot", scale="utc").jd,
    t_max_jd=Time.now().jd,
    with_history: bool = False,
    with_cutouts: bool = False,
    chunk_size: int = 500,
    resume_token: str = None,
    warn_exceeding_chunk: bool = True,
    program_id: int = None,
    logger=None,
) -> tuple:
    """
    Function to query ampel based on lists of healpix pixels of different resolution (nside is the respective resolution)
    """

    if logger is None:
        logger = logging.getLogger(__name__)

    queryurl_skymap = API_ZTF_ARCHIVE_URL + f"/alerts/healpix/skymap"

    # First, we create a json body to post
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {get_ampel_token()}",
    }
    if with_history:
        hist = "true"
    else:
        hist = "false"

    if with_cutouts:
        cutouts = "true"
    else:
        cutouts = "false"

    # if we have a resume_token to proceed to the next chunk, that's all we need
    if resume_token is not None:
        queryurl_stream = API_ZTF_ARCHIVE_URL + f"/stream/{resume_token}/chunk"
        response = requests.get(
            queryurl_stream, params={"with_history": hist}, headers=headers
        )

    # if we don't have a resume_token, we first need to create the full query
    else:
        # Now we reduce the query size
        regions = utils.deres(nside=nside, ipix=pixels)

        n_pix = 0
        for reg in regions:
            n_pix += len(reg["pixels"])

        logger.debug(f"This comprises {n_pix} individual pixels")

        if n_pix > MAX_N_PIX:
            logger.warning(
                f"Total number of pixels exceeds threshold ({MAX_N_PIX} pixels). Issuing a query for the full sky instead."
            )
            regions = [
                {"nside": 1, "pixels": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]}
            ]

        query = {
            "regions": regions,
            "jd": {
                "$lt": t_max_jd,
                "$gt": t_min_jd,
            },
            "latest": "false",
            "with_history": hist,
            "with_cutouts": cutouts,
            "chunk_size": chunk_size,
        }

        if resume_token:
            query["resume_token"] = resume_token

        if program_id is not None:
            query["programid"] = program_id

        logger.debug(f"Query url:\n{queryurl_skymap}")
        logger.debug(f"Query:\n{query}")

        response = requests.post(url=queryurl_skymap, json=query, headers=headers)

    logger.debug(response)
    logger.debug(response.status_code)

    if response.status_code in [424, 503]:
        raise requests.exceptions.RequestException

    try:
        res_json = response.json()
        remaining_chunks = res_json["remaining"]["chunks"]
        logger.debug(f"Remaining chunks: {remaining_chunks}")
        chunk_id = res_json.get("chunk", None)
        resume_token = response.json().get("resume_token", None)
        query_res = [i for i in response.json()["alerts"]]
    except JSONDecodeError:
        if response.headers:
            logger.debug(response.headers)
        raise requests.exceptions.RequestException

    nr_results = len(query_res)

    logger.debug(f"Found {nr_results} alerts.")

    if nr_results == chunk_size and warn_exceeding_chunk:
        logger.warning(
            f"Query result limited by chunk size! You will most likely be missing alerts!"
        )

    return query_res, resume_token, chunk_id, remaining_chunks


def ampel_api_skymap(
    t_min: Time,
    t_max: Time,
    cone_nside: int,
    cone_ids: list[int],
):
    logger.info("Commencing skymap scan")

    query_res = []

    resume = True
    chunk_size = 2000
    resume_token = None

    i = 0
    total_chunks = 0
    t0 = time.time()

    while resume:
        res, resume_token, chunk_id, remaining_chunks = ampel_api_skymap_single(
            pixels=cone_ids,
            nside=cone_nside,
            t_min_jd=t_min.jd,
            t_max_jd=t_max.jd,
            logger=logger,
            chunk_size=chunk_size,
            resume_token=resume_token,
            warn_exceeding_chunk=False,
        )
        query_res.extend(res)

        ampel_api_acknowledge_chunk(resume_token=resume_token, chunk_id=chunk_id)

        if i == 0:
            total_chunks = remaining_chunks + 1
            logger.info(f"Total chunks: {total_chunks}")

        if remaining_chunks % 50 == 0 and remaining_chunks != 0:
            t1 = time.time()
            processed_chunks = total_chunks - remaining_chunks
            time_per_chunk = (t1 - t0) / processed_chunks
            remaining_time = time_per_chunk * remaining_chunks
            logger.info(
                f"Remaining chunks: {remaining_chunks}. Estimated time to finish: "
                f"{remaining_time / 60:.0f} min"
            )

        if len(res) < chunk_size:
            resume = False
            logger.info("Done.")
        else:
            logger.debug(f"Chunk size reached ({chunk_size}), commencing next query.")
        i += 1

    return query_res


def get_preprocessed_results(file_basename: str, logger=None) -> None | list:
    """
    Access the DESY Cloud to look if there are precomputed results from an AMPEL run
    there
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    desy_cloud_token = load_credentials("desy_cloud_token", token_based=True)

    filename = PREPROCESSED_CACHE_DIR.joinpath(f"{file_basename}.json.gz")

    res = requests.get(
        f"https://syncandshare.desy.de/public.php/webdav/{filename.name}",
        headers={"X-Requested-With": "XMLHttpRequest"},
        auth=(desy_cloud_token, "bla"),
    )

    if res.status_code != 200:
        logger.warning(
            "\n\n-------------------- !! -------------------\nSomething went wrong with your query.\nCheck your credentials and make sure Ampel\nhas run correctly at Desy.\n-------------------- !! -------------------\n\n"
        )
        return None

    with open(f"{filename}", "wb") as f:
        f.write(res.content)

    res = []
    with gzip.open(filename, "rb") as f_in:
        data = load(f_in)
        for t in data:
            ztf_id = ZTFIdMapper.to_ext_id(t.stock.get("stock"))
            pp = t.get_photopoints()
            pp_reformatted = utils.reformat_downloaded_results(
                photopoints=pp, ztf_id=ztf_id
            )
            redshifts = t.get_latest_t2_body(unit="T2DigestRedshifts")
            kilonova_eval = t.get_latest_t2_body(unit="T2KilonovaEval")
            pp_reformatted.update({"kilonova_eval": kilonova_eval})
            pp_reformatted.update({"redshifts": redshifts})
            res.append(pp_reformatted)

    return res
