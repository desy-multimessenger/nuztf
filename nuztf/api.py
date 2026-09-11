import logging
import os

import dotenv
from astropy.time import Time

from nuztf.ampel import (
    ampel_api_cutout,
    ampel_api_name,
    ampel_api_skymap,
    ensure_ampel_cutouts,
)
from nuztf.boom import (
    boom_api_cutout,
    boom_api_name,
    boom_api_skymap,
)
from nuztf.kowalski import (
    ensure_kowalski_cutouts,
    kowalski_api_cutout,
    kowalski_api_name,
    kowalski_api_skymap,
)

OVERWRITE = False

# Load environment variables from .env file
dotenv.load_dotenv()

ZTF_BACKEND = str(os.getenv("ZTF_BACKEND", "ampel")).lower()
assert ZTF_BACKEND in [
    "ampel",
    "kowalski",
    "boom",
], f"Invalid ZTF backend: {ZTF_BACKEND}"


def api_name(
    ztf_name: str,
    with_history: bool = True,
    with_cutouts: bool = False,
    limit: int = 999999,
    backend: str = ZTF_BACKEND,
) -> list:
    """
    Get alert data from the specified backend.

    :param ztf_name: Name of source
    :param with_history: Whether to include history
    :param with_cutouts: Whether to include cutouts
    :param limit: Limit for the number of alerts
    :param backend: Backend to use for fetching data ("ampel" or "kowalski" or "boom")

    :return: Alert data
    """
    if backend == "ampel":
        return ampel_api_name(
            ztf_name,
            with_history=with_history,
            with_cutouts=with_cutouts,
            limit=limit,
        )
    elif backend == "kowalski":
        return kowalski_api_name(
            ztf_name,
            with_cutouts=with_cutouts,
        )
    elif backend == "boom":
        return boom_api_name(
            ztf_name,
            with_cutouts=with_cutouts,
        )
    else:
        raise ValueError(f"Backend {backend} not supported")


def api_skymap(
    t_min: Time,
    t_max: Time,
    cone_nside: int,
    cone_ids: list[int],
    backend: str = ZTF_BACKEND,
) -> list:
    """
    Get skymap data from the specified backend.

    :param t_min: Start time
    :param t_max: End time
    :param cone_nside: Nside for the skymap
    :param cone_ids: List of cone IDs
    :param backend: Backend to use for fetching data ("ampel" or "kowalski" or "boom")

    :return: Skymap data
    """
    if backend == "ampel":
        return ampel_api_skymap(
            t_min=t_min,
            t_max=t_max,
            cone_nside=cone_nside,
            cone_ids=cone_ids,
        )
    elif backend == "kowalski":
        return kowalski_api_skymap(
            t_min=t_min,
            t_max=t_max,
            cone_nside=cone_nside,
            cone_ids=cone_ids,
        )
    elif backend == "boom":
        return boom_api_skymap(
            t_min=t_min,
            t_max=t_max,
            cone_nside=cone_nside,
            cone_ids=cone_ids,
        )
    else:
        raise ValueError(f"Backend {backend} not supported")


def api_cutout(
    candid: int,
    backend: str = ZTF_BACKEND,
) -> dict:
    """
    Get cutout data from the specified backend.
    :param candid: Candid of the alert
    :param backend: Backend to use for fetching data ("ampel" or "kowalski" or "boom")
    :return: Cutout data
    """

    if backend == "ampel":
        return ampel_api_cutout(candid=candid)
    elif backend == "kowalski":
        return kowalski_api_cutout(candid=candid)
    elif backend == "boom":
        return boom_api_cutout(candid=candid)
    else:
        raise ValueError(f"Backend {backend} not supported")


def ensure_cutouts(
    alert: list,
    backend: str = ZTF_BACKEND,
):
    """
    Ensure cutouts for the alert data.

    :param alert: Alert data
    :param backend: Backend to use for fetching data ("ampel" or "kowalski" or "boom")
    :return: Alert data with cutouts
    """

    if backend == "ampel":
        return ensure_ampel_cutouts(alert)
    elif backend == "kowalski":
        return ensure_kowalski_cutouts(alert)
    elif backend == "boom":
        return ensure_boom_cutouts(alert)
    else:
        raise ValueError(f"Backend {backend} not supported")
