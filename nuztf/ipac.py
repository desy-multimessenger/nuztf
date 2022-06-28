import subprocess
import os
import logging
import numpy as np
import pandas as pd
from ztfquery.io import LOCALSOURCE
from nuztf.credentials import load_credentials
from nuztf.ampel_api import ampel_api_name, calculate_mean_position

logger = logging.getLogger(__name__)


global_fp_user, global_fp_password = load_credentials("ipac_fp_global")
ztf_fp_user, ztf_fp_password = load_credentials("ipac_fp_personal")

fp_cache_dir = os.path.join(LOCALSOURCE, "cache/")


def check_wget():
    try:
        subprocess.run("command -v wget", check=True, capture_output=True, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "'wget' is not installed, but is required for requesting forced photometry."
        )


def request_forced_photometry(
    ztf_name: str, jd_start: float = None, jd_end: float = None
):
    check_wget()

    res = ampel_api_name(ztf_name, with_history=True, with_cutouts=False)

    ra, dec = calculate_mean_position(res)

    cmd = f'wget --http-user={global_fp_user} --http-passwd={global_fp_password} -O log.txt "https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?'
    cmd += f"ra={ra}&dec={dec}&"
    if jd_start is not None:
        cmd += f"jdstart={jd_start}&"
    if jd_end is not None:
        cmd += f"jdend={jd_end}&"
    cmd += f'email={ztf_fp_user}&userpass={ztf_fp_password}"'

    subprocess.run(cmd, check=True, capture_output=True, shell=True)


def plot_forced_photometry(ztf_name: str):

    cache_file = os.path.join(fp_cache_dir, f"{ztf_name}.txt")

    if not os.path.exists(cache_file):
        raise Exception(
            f"Failed to file {cache_file}. "
            f"Please download forced photometry for the source, and save it to this path"
        )

    df = pd.read_csv(cache_file, comment="#", sep=" ")

    acceptable_proc = [0, 62, 63, 255]

    mask = np.array([x in acceptable_proc for x in df["procstatus"]])

    logger.info(
        f"Found {np.sum(mask)} entries with procstatus in {acceptable_proc}. "
        f"Dropping {np.sum(~mask)} other entries."
    )

    df = df[mask]

    return df
