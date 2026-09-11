"""
This module contains functions to interact with the BOOM API.
"""

from blastwave import LSSTClient, ZTFClient

from nuztf.kowalski.config import fp_mapping, get_kowalski


def boom_api_name(
    source_name: str,
    with_cutouts: bool = False,
    client: ZTFClient | LSSTClient | None = None,
):
    """
    Function to interact with the BOOM API.

    :param source_name: Source name
    :with_cutouts: Boolean to indicate if cutouts should be used
    :client: BoomClient to interact with the BOOM API
    """

    if client is None:
        is_ztf = "ztf" in source_name.lower()
        client = ZTFClient() if is_ztf else LSSTClient()

    src = client.get_source(object_id)
    return src.convert_to_ztfstyle()
