_LSST_FILTER = {
    "candidate.reliability": {"$gt": 0.4},
    "candidate.isDipole": {"$eq": False},
    "candidate.pixelFlags": {"$eq": False},
    "candidate.glint_trail": {"$eq": False},
    "candidate.isNegative": {"$eq": False},
    "candidate.centroid_flag": {"$eq": False},
    "candidate.psfFlux_flag": {"$eq": False},
}
_ZTF_BASE_FILTER = {"candidate.isdiffpos": {"$eq": True}}
_ZTF_DRB_FILTER = {"candidate.drb": {"$gt": 0.65}}
_ZTF_RB_FILTER = {"candidate.rb": {"$gt": 0.5}}

ZTF_DRB_JD = 2458653.5


def get_quality_filter(catalog: str, jd: float) -> dict:
    """
    Get quality filter based on survey properties

    :param catalog: Catalog
    :param jd: JD of search window
    :return:
    """
    # LSST - reliability
    if catalog == "LSST_alerts":
        return {**_LSST_FILTER}
    # ZTF - drb or rb depending on date
    elif catalog == "ZTF_alerts":
        if jd > ZTF_DRB_JD:
            return {**_ZTF_BASE_FILTER, **_ZTF_DRB_FILTER}
        else:
            return {**_ZTF_BASE_FILTER, **_ZTF_RB_FILTER}
    else:
        raise ValueError(f"Unrecognized survey: {catalog}")
