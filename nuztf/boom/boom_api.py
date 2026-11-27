"""
This module contains functions to interact with the BOOM API.
"""

from blastwave import LSSTClient, ZTFClient
from blastwave.models.query import BOOMQuery

MIN_FP_SNR = 3.0

#
# def boom_api_name(
#     source_name: str,
#     client: ZTFClient | LSSTClient | None = None,
# ) -> int:
#     """
#     Function to get the candid from the candidate name.
#
#     :param source_name: Name of the source
#     :param client: ZTFClient or LSSTClient to interact with the BOOM API
#     :return: Candid of the source
#     """
#     query = BOOMQuery(
#         catalog_name=client.catalog,
#         filter={"objectId": {"$eq": source_name}},
#         projection={"_id": 1},
#         limit=1
#     )
#     alert = client.query(query)[0]
#     print(alert)
#     raise


def boom_api_candid(
    candidate: dict,
    with_cutouts: bool = False,
    client: ZTFClient | LSSTClient | None = None,
):
    """
    Function to interact with the BOOM API.

    :param candidate: candidate dictionary containing the candidate information
    :with_cutouts: Boolean to indicate if cutouts should be used
    :client: BoomClient to interact with the BOOM API
    """

    object_id = candidate["objectId"]

    if client is None:
        is_ztf = "ztf" in object_id.lower()
        client = ZTFClient() if is_ztf else LSSTClient()

    query = BOOMQuery(
        catalog_name=client.catalog,
        filter={"_id": {"$eq": candidate["candid"]}},
    )
    alert = client.query(query)[0]
    alert["candid"] = candidate["candid"]

    query = BOOMQuery(
        catalog_name=client.aux_table,
        filter={"_id": {"$eq": object_id}},
    )
    aux_data = client.query(query)[0]

    prv_candidates = [
        x for x in aux_data["prv_candidates"] if x["jd"] < alert["candidate"]["jd"]
    ]

    jds = [x["jd"] for x in prv_candidates]

    fp_hists = [
        x
        for x in aux_data["fp_hists"]
        if (
            (x["jd"] < alert["candidate"]["jd"])
            & (x["jd"] not in jds)
            & (x["psfFlux"] / x["psfFluxErr"] > MIN_FP_SNR)
        )
    ]

    prv_candidates += fp_hists

    alert["prv_candidates"] = prv_candidates

    return [alert]
