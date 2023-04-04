#!/usr/bin/env python
# coding: utf-8

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rc

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

rc("font", **{"family": "serif", "serif": ["Palatino"]})

REJECTION_DIR = {
    0: "Pre-ZTF",
    1: "Alert\nretraction",
    2: "Proximity\nto the sun",
    3: "Telescope\nmaintenance",
    4: "Alert\nquality",
    5: "Southern sky",
    6: "Low altitude",
    7: "Bad weather",
    8: "Galactic plane",
    9: "Other",
}


def get_rectangular_area(
    ra_err: list | np.ndarray, dec_err: list | np.ndarray, dec: float
) -> float:
    rectangular_area = (
        (np.max(ra_err) - np.min(ra_err))
        * (np.max(dec_err) - np.min(dec_err))
        * abs(np.cos(np.radians(dec)))
    )
    return rectangular_area


def print_statistics(df: pd.DataFrame):
    """
    Print population statistics
    """
    if "Observed area (corrected for chip gaps)" in df.keys():
        df_type = "fu"
    else:
        df_type = "nofu"

    stats: dict[str, float | int | str] = {}
    stats["type"] = df_type

    stats["summed_signalness"] = df.Signalness.sum()

    if df_type == "nofu":
        areas: list[float] = []
        for i, ra_err in enumerate(df["RA Unc (rectangle)"].values):
            if not isinstance(ra_err, float) and not isinstance(
                dec_err := df["Dec Unc (rectangle)"].values[i], float
            ):
                ra_err = np.fromstring(ra_err.strip("[").strip("]"), sep=",")
                dec_err = np.fromstring(dec_err.strip("[").strip("]"), sep=",")
                dec = float(np.asarray(df["Dec"].values[i]))

                areas.append(get_rectangular_area(ra_err, dec_err, dec))
            else:
                areas.append(np.nan)
        df["Area (rectangle)"] = areas

        stats["alerts_unretracted"] = len(df.query("Code != 1"))
        stats["alerts_unretracted_ztf"] = len(df.query("Code != 0 and Code != 1"))

    stats["summed_area"] = df["Area (rectangle)"].sum()

    if df_type == "fu":
        stats["summed_obs_area"] = df["Observed area (corrected for chip gaps)"].sum()

    stats["alerts"] = len(df)

    print("---------------")
    for k, val in stats.items():
        print(f"{k}: {val}")
    print("---------------")

    return stats


def plot_overview(
    df_fu: pd.DataFrame,
    df_nofu: pd.DataFrame,
    include_retracted: bool = False,
    include_preztf: bool = False,
):
    """
    Plot pie chart for observation rejection reasons
    """
    counts = [n := len(df_fu)]
    labels = [f"Total\nfollow-up ({n})"]

    _df = df_nofu.copy(deep=True)

    if not include_preztf:
        _df.query("Code != 0", inplace=True)
    if not include_retracted:
        _df.query("Code != 1", inplace=True)

    _df.replace({"Code": {5: 9, 6: 9, 7: 9}}, inplace=True)

    for code in _df.Code.unique():
        counts.append(n := len(_df.query("Code == @code")))
        labels.append(f"{REJECTION_DIR[code]} ({n})")

    indices = list(np.argsort(counts))
    indices.reverse()
    counts = np.asarray(counts)[indices]
    labels = np.asarray(labels)[indices]

    explode = np.empty(len(counts))
    explode.fill(0.02)
    explode[0] = 0.1

    palette_color = sns.color_palette("dark")

    fig, ax = plt.subplots(figsize=(width := 4.5, width / 1.62))

    patches, texts, pcts = ax.pie(
        counts,
        labels=labels,
        explode=explode,
        startangle=90,
        labeldistance=1.13,
        autopct="%.0f%%",
        pctdistance=0.7,
    )
    for i, patch in enumerate(patches):
        texts[i].set_color(patch.get_facecolor())
    plt.setp(pcts, color="white", fontweight="bold", fontsize=10)
    plt.setp(texts, fontweight="bold")
    plt.tight_layout()

    outfile_base = "follow_up_overview"
    if include_preztf:
        outfile_base += "_preztf"
    if include_retracted:
        outfile_base += "_retracted"

    plt.savefig(outfile_base + ".pdf")
    plt.close()


infile_fu = Path.cwd() / "overview" / "neutrino_too_followup - OVERVIEW_FU.csv"
infile_nofu = Path.cwd() / "overview" / "neutrino_too_followup - OVERVIEW_NOT_FU.csv"

df_fu = pd.read_csv(infile_fu)
df_nofu = pd.read_csv(infile_nofu)

print_statistics(df=df_nofu)
print_statistics(df=df_fu)

plot_overview(df_fu=df_fu, df_nofu=df_nofu)
# plot_overview(df_fu=df_fu, df_nofu=df_nofu, include_preztf=True, include_retracted=True)
