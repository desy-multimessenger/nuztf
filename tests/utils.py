"""
Util functions for tests
"""

import numpy as np

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.skymap_scanner import SkymapScanner


def clip_scanner_data(
    scanner: SkymapScanner | NeutrinoScanner,
    max_jd: float | None = None,
) -> SkymapScanner | NeutrinoScanner:
    """
    Clip scanner data to only include detections before a certain JD.
    This is useful for testing, as it allows us to use old data
    without new detections changing the confidence interval.

    :param scanner: Scanner object
    :param max_jd: float - maximum JD to clip
    :return:
    """

    if max_jd is None:
        max_jd = scanner.default_t_max.jd

    for name, res in sorted(scanner.cache_candidates.items()):
        # Only use old data, so new detections do not change CI
        dets = [
            x
            for x in res["prv_candidates"] + [res["candidate"]]
            if ("isdiffpos" in x.keys()) & (x["jd"] < max_jd)
        ]
        dets = [x for x in dets if str(x["isdiffpos"]).lower() in ["t", "true", "1"]]
        max_jd_idx = np.argmax([x["jd"] for x in dets])
        cand = dets[max_jd_idx]
        res["candidate"] = cand
        res["prv_candidates"] = [x for x in dets if x != cand]

    return scanner
