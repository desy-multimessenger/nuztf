# ztf_ForcedPhotometryRequestExample.py Updated: 2022-06-02
# Written for Python 3.9.7, likely compatible with other python3 versions
# Adapted by Robert Stein, originally written by Yashvi Sharma

# requests is used to query the Forced Photometry Service
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
from nuztf.credentials import load_credentials
from astropy.time import Time
import os
import logging
from ztfquery.io import LOCALSOURCE
from astropy.table import Table

from nuztf.ampel_api import ampel_api_name, calculate_mean_position
from nuztf.cat_match import context_from_name
from nuztf.plot import lightcurve_from_science_image
from nuztf.irsa import get_irsa_path

logger = logging.getLogger(__name__)

# Example to place a single request
# Easily iterable with for/while loops, though the service limits users to 100 requests at a time
FORCED_PHOTOMETRY_BASE_URL = "https://ztfweb.ipac.caltech.edu"

# URL used to place a request
IPAC_URL = f"{FORCED_PHOTOMETRY_BASE_URL}/cgi-bin/requestForcedPhotometry.cgi"

# URL for checking the status of jobs sent by user
# Note: This webpage only updates once an hour, on the hour
STATUS_URL = f"{FORCED_PHOTOMETRY_BASE_URL}/cgi-bin/getForcedPhotometryRequests.cgi"

sleep_time_s = 60.0

ipac_global_user, ipac_global_password = load_credentials("ipac_fp_global")
ztf_fp_user, ztf_fp_password = load_credentials("ipac_fp_personal")

fp_cache_dir = os.path.join(LOCALSOURCE, "cache/")


def submit_request(
    ra_deg: float,
    dec_deg: float,
    start_jd: float = 2458194.5,
    end_jd: float = None,
    save_path: str = None,
):

    if end_jd is None:
        end_jd = Time.now().jd

    # Send request and return unformatted HTML output of get() method
    request = requests.get(
        IPAC_URL,
        auth=(ipac_global_user, ipac_global_password),
        params={
            "ra": ra_deg,
            "dec": dec_deg,
            "jdstart": start_jd,
            "jdend": end_jd,
            "email": ztf_fp_user,
            "userpass": ztf_fp_password,
        },
    )

    if request.status_code not in [200]:
        logger.error(f"Error with request:\n {request}")

    # Formatted (Parseable) HTML of request
    req_soup = BeautifulSoup(request.text, "html.parser")

    # For one pending job (ie: one request), check forced photometry job tables

    req_ra = req_dec = req_jds = req_jde = None

    # Grab the table section of the HTML
    req_table = req_soup.find("table")
    # Limit to rows in said table
    for row in req_table.find_all("tr"):
        # Find the items in the row (omit header cells)
        cols = row.find_all("td")
        if len(cols) > 0:
            # Find RA, Dec, and JD values of request as recorded by the service
            # Use a nested for loop here if submitting multiple requests, append values to list/arr/etc
            req_ra = float(cols[0].text.strip())
            req_dec = float(cols[1].text.strip())
            req_jds = float(cols[2].text.strip())
            req_jde = float(cols[3].text.strip())

    if req_ra is None:
        raise Exception(f"Error parsing output of table: \n {req_table}")

    # Check if request is fulfilled using a request status check and parse with beautifulsoup4
    # Note: Requests older than 30 days will not show up here
    # Iterable, as beautifulsoup can parse html columns and output lists

    request_completed = False

    output_end_time = output_lc_path = ""

    # Open loop to periodically check the Forced Photometry job status page
    while not request_completed:
        # Query service for jobs sent in past 30 days
        output_status = requests.get(
            STATUS_URL,
            auth=(ipac_global_user, ipac_global_password),
            params={
                "email": ztf_fp_user,
                "userpass": ztf_fp_password,
                "option": "All recent jobs",
                "action": "Query Database",
            },
        )

        # Check if job has ended and lightcurve file was created
        # Note: If an exotic error occurs and the ended field is not populated, this will go on forever
        # Table has 11 cols, reqid=0, ended=7, lc=10

        # Format HTML
        output_soup = BeautifulSoup(output_status.text, "html.parser")
        # Get Job information in HTML table
        output_table = output_soup.find("table")

        # Parse Table rows
        if output_table not in ["", None]:
            for row in output_table.find_all("tr"):
                # Parse Table entries
                cols = row.find_all("td")
                if len(cols) > 0:
                    # Check if values contained in a given row coorespond to the current request
                    # Use a nested for loop here to check all requests submitted if there are more than one pending
                    output_ra = float(cols[1].text.strip())
                    output_dec = float(cols[2].text.strip())
                    output_jds = float(cols[3].text.strip())
                    output_jde = float(cols[4].text.strip())
                    output_end_time = cols[7].text.strip()

                    # Check if job is finished (output_end_time)
                    # Check for equality between recorded request params and what is known
                    # to be in the Job Status table, accounting for rounding
                    if (
                        output_end_time != ""
                        and abs(output_ra - req_ra) <= 0.00001
                        and abs(output_dec - req_dec) <= 0.00001
                        and abs(output_jds - req_jds) <= 0.00001
                        and abs(output_jde - req_jde) <= 0.00001
                    ):
                        # Get end time of job and lightcurve path
                        output_lc_path = cols[10].text.strip()
                        output_end_time = cols[7].text.strip()
                        # Set open loop hooks
                        request_completed = True

        if not request_completed:
            logger.info(
                f"Lightcurve not found, sleeping for {sleep_time_s} seconds. "
                f"Warning: the lightcurve webpage is only updated hourly."
            )
            time.sleep(sleep_time_s)

    logger.info(f"Lightcurve found! Job ended at: {output_end_time}")

    # If Job Status table values are not null,
    # set the path for file to be downloaded to the path recorded in Job Status Table
    if output_lc_path in [""]:
        raise Exception("Job is done but no lightcurve was produced, quitting")

    # Set link to download lightcurve from
    forced_photometry_req_url = f"{FORCED_PHOTOMETRY_BASE_URL}{output_lc_path}"

    # Set local path and filename for lightcurve (currently same as what is assigned by the service)
    if save_path is None:
        local_filename = forced_photometry_req_url.split("/")[-1]
        save_path = os.path.join(fp_cache_dir, local_filename)

    # Download the file with a get request
    with requests.get(
        forced_photometry_req_url,
        stream=True,
        auth=(ipac_global_user, ipac_global_password),
    ) as r:

        logger.info(f"Saving to {save_path}")

        # Write to the local file in chunks in case file is large
        with open(save_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


# fp_keymap = {
#     ""
# }


def load_forced_photometry(
    save_path: str, use_difference_flux: bool = True, snr_cut: float = 5.0
) -> pd.DataFrame:

    df = pd.read_csv(save_path, comment="#", sep=" ")

    acceptable_proc = ["0", "62", "63", "255"]
    mask = np.array([str(x) in acceptable_proc for x in df["procstatus"]])

    logger.info(
        f"Found {np.sum(mask)} entries with procstatus in {acceptable_proc}. "
        f"Dropping {np.sum(~mask)} other entries."
    )

    df = df[mask]
    df["mjd"] = df["jd,"] - 2400000.5
    df["filtercode"] = [f"z{x[-1].lower()}" for x in df["filter,"]]

    if not use_difference_flux:
        # Frank's formulae
        nearestrefflux = 10.0 ** (0.4 * (df["zpdiff,"] - df["nearestrefmag,"]))
        nearestreffluxunc = df["nearestrefmagunc,"] * nearestrefflux / 1.0857

        Fluxtot = df["forcediffimflux,"] + nearestrefflux
        Fluxunctot = np.sqrt(df["forcediffimflux,"] ** 2.0 - nearestreffluxunc**2.0)
        SNRtot = Fluxtot / Fluxunctot

        df["mag"] = df["zpdiff,"] - 2.5 * np.log10(Fluxtot)
        df["magerr"] = 1.0857 / SNRtot

        snr_mask = SNRtot > snr_cut
        df = df[snr_mask]

        logger.info(df)

    else:
        raise

    return df


def get_fp_cache_name(
    source_name: str,
) -> str:
    return os.path.join(fp_cache_dir, f'{source_name.replace(" ", "")}_fp.csv')


def plot_forced_photometry_lightcurve(
    source_name: str,
    start_jd: float = 2458194.5,
    end_jd: float = None,
    use_difference_flux: bool = True,
    snr_cut: float = 5.0,
    overwrite_cached_request: bool = False,
    nu_name: list = None,
    plot_mag: bool = False,
    atel: bool = True,
    extra_folder: str = None,
    check_obs: bool = False,
    check_obs_lookback_weeks: float = 4.0,
    expanded_labels: bool = True,
    ylim: tuple = None,
):

    # Query IRSA, or load from cache

    try:
        os.makedirs(fp_cache_dir)
    except OSError:
        pass

    cache_path = get_fp_cache_name(source_name)

    if not np.logical_and(os.path.exists(cache_path), not overwrite_cached_request):

        if "ZTF" in source_name:
            res = ampel_api_name(source_name, with_history=True, with_cutouts=False)
            ra_deg, dec_deg = calculate_mean_position(res)

        else:
            (ra_deg, dec_deg), _, _ = context_from_name(source_name=source_name)

        logger.info(
            f"Requesting IRSA photometry (rerequest is {overwrite_cached_request}, "
            f"cached file does {['not', ''][os.path.exists(cache_path)]}. "
            f"This will take some time!"
        )

        submit_request(
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            start_jd=start_jd,
            end_jd=end_jd,
            save_path=cache_path,
        )

    if not os.path.exists(cache_path):
        raise Exception(f"Something went wrong, file {cache_path} does not exist.")

    logger.debug(f"Loading from {cache_path}")

    df = load_forced_photometry(
        cache_path, use_difference_flux=use_difference_flux, snr_cut=snr_cut
    )
    data = Table.from_pandas(df)

    lightcurve_from_science_image(
        source_name=source_name,
        data=data,
        nu_name=nu_name,
        plot_mag=plot_mag,
        atel=atel,
        extra_folder=extra_folder,
        check_obs=check_obs,
        check_obs_lookback_weeks=check_obs_lookback_weeks,
        expanded_labels=expanded_labels,
        ylim=ylim,
    )


def plot_composite_lightcurve(source_name: str, snr_cut: float = 5.0, **kwargs):
    fp_path = get_fp_cache_name(source_name)
    fp_df = load_forced_photometry(fp_path, use_difference_flux=False, snr_cut=snr_cut)
    irsa_path = get_irsa_path(source_name)
    irsa_df = pd.read_csv(irsa_path)

    mask = fp_df["mjd"] > max(irsa_df["mjd"])

    df = pd.concat([irsa_df, fp_df[mask]], ignore_index=True)

    data = Table.from_pandas(df)
    lightcurve_from_science_image(source_name=source_name, data=data, **kwargs)


if __name__ == "__main__":
    logging.getLogger("nuztf").setLevel("DEBUG")
    plot_forced_photometry_lightcurve(
        "WISEA J145820.77+412101.9",
        plot_mag=True,
        use_difference_flux=False,
        snr_cut=10,
    )
