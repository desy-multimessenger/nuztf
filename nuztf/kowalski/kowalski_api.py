"""
This module contains functions to interact with the Kowalski API.
"""

from penquins import Kowalski

from nuztf.kowalski.config import fp_mapping, get_kowalski


def kowalski_api_name(
    ztf_name: str,
    with_cutouts: bool = False,
    kowalski: Kowalski | None = None,
):
    """
    Download alert data from Kowalski

    :param ztf_name: Name of source
    :param with_cutouts: Whether to include cutouts
    :param kowalski: Kowalski object
    :return: Alert data
    """
    if kowalski is None:
        kowalski = get_kowalski()

    query_config = {
        "query_type": "find",
        "query": {
            "catalog": "ZTF_alerts",
            "filter": {
                "objectId": {"$eq": ztf_name},
            },
            "projection": {
                "_id": 0,
                "cutoutScience": int(with_cutouts),
                "cutoutTemplate": int(with_cutouts),
                "cutoutDifference": int(with_cutouts),
                "coordinates": 0,
            },
        },
    }

    query_result = kowalski.query(query_config)

    if "data" in query_result:
        alerts = query_result["data"]
    else:
        alerts = query_result.get("default").get("data")

    jds = [x["candidate"]["jd"] for x in alerts]

    max_idx = jds.index(max(jds))
    latest_alert = alerts[max_idx]

    jds = [latest_alert["candidate"]["jd"]]

    prv_alerts = []

    for prv_cand in alerts:
        if (prv_cand["candidate"]["jd"] not in jds) & (
            "magpsf" in prv_cand["candidate"]
        ):
            jds.append(prv_cand["candidate"]["jd"])
            prv_alerts.append(prv_cand["candidate"])

    # Query for prv_candidates/forced photometry

    query_config = {
        "query_type": "find",
        "query": {
            "catalog": "ZTF_alerts_aux",
            "filter": {
                "_id": {"$eq": ztf_name},
            },
            "projection": {"cross_matches": 0},
        },
    }

    query_result = kowalski.query(query_config)
    if "data" in query_result:
        out = query_result["data"]
    else:
        out = query_result.get("default").get("data")

    if len(out) > 0:
        for prv_cand in out[0]["prv_candidates"]:
            if (prv_cand["jd"] not in jds) & ("magpsf" in prv_cand):
                jds.append(prv_cand["jd"])
                prv_alerts.append(prv_cand)

        # Add forced photometry

        fp_dets = []

        if "fp_hists" in out[0]:
            for fp_dict in out[0]["fp_hists"]:
                if (
                    (fp_dict["jd"] not in jds)
                    & ("mag" in fp_dict)
                    & ("magerr" in fp_dict)
                ):
                    if fp_dict["snr"] > 3.0:
                        for old_key, new_key in fp_mapping.items():
                            fp_dict[new_key] = fp_dict.pop(old_key)
                        fp_dict["isdiffpos"] = "t"
                        fp_dict["fp_bool"] = True
                        fp_dets.append(fp_dict)

        jds += [x["jd"] for x in fp_dets]
        prv_alerts += fp_dets

    latest_alert["prv_candidates"] = prv_alerts
    latest_alert["candidate"]["jdstarthist"] = min(jds)

    if not len(jds) == len(set(jds)):
        raise ValueError(f"Duplicate JDs for {ztf_name}")
    return [latest_alert]
