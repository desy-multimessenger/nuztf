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
from nuztf.kowalski import (
    ensure_kowalski_cutouts,
    kowalski_api_cutout,
    kowalski_api_name,
    kowalski_api_skymap,
)

OVERWRITE = False

# Load environment variables from .env file
dotenv.load_dotenv()

ZTF_BACKEND = os.getenv("ZTF_BACKEND", "ampel")
assert ZTF_BACKEND in ["ampel", "kowalski"], f"Invalid ZTF backend: {ZTF_BACKEND}"


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
    :param backend: Backend to use for fetching data ("ampel" or "kowalski")

    :return: Alert data
    """
    if backend == "ampel":
        return ampel_api_name(
            ztf_name,
            with_history=with_history,
            with_cutouts=with_cutouts,
            limit=limit,
        )
    else:
        return kowalski_api_name(
            ztf_name,
            with_cutouts=with_cutouts,
        )


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
    :param backend: Backend to use for fetching data ("ampel" or "kowalski")

    :return: Skymap data
    """
    if backend == "ampel":
        return ampel_api_skymap(
            t_min=t_min,
            t_max=t_max,
            cone_nside=cone_nside,
            cone_ids=cone_ids,
        )
    else:
        return kowalski_api_skymap(
            t_min=t_min,
            t_max=t_max,
            cone_nside=cone_nside,
            cone_ids=cone_ids,
        )


def api_cutout(
    candid: int,
    backend: str = ZTF_BACKEND,
) -> dict:
    """
    Get cutout data from the specified backend.
    :param candid: Candid of the alert
    :param backend: Backend to use for fetching data ("ampel" or "kowalski")
    :return: Cutout data
    """

    if backend == "ampel":
        return ampel_api_cutout(candid=candid)
    else:
        return kowalski_api_cutout(candid=candid)


def ensure_cutouts(
    alert: list,
    backend: str = ZTF_BACKEND,
):
    """
    Ensure cutouts for the alert data.

    :param alert: Alert data
    :param backend: Backend to use for fetching data ("ampel" or "kowalski")
    :return: Alert data with cutouts
    """

    if backend == "ampel":
        return ensure_ampel_cutouts(alert)
    else:
        return ensure_kowalski_cutouts(alert)
