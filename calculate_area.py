#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np
import pandas as pd
from ampel.log.AmpelLogger import AmpelLogger

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.parse_nu_gcn import find_gcn_no, get_latest_gcn, parse_gcn_circular

RECALC = True
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if RECALC:
    rectangular_areas = []
    corrected_areas = []
    uncorrected_areas = []
    good_neutrinos = []
    bad_neutrinos = []

    neutrinos = [
        "IC190503A",
        "IC190619A",
        "IC190730A",
        "IC190922B",
        "IC191001A",
        "IC200107A",
        "IC200109A",
        "IC200117A",
        "IC200512A",
        "IC200530A",
        "IC200620A",
        "IC200916A",
        "IC200926A",
        "IC200929A",
        "IC201007A",
        "IC201021A",
        "IC201130A",
        "IC201209A",
        "IC201222A",
        "IC210210A",
        "IC210510A",
        "IC210629A",
        "IC210811A",
        "IC210922A",
        "IC220405A",
        "IC220405B",
        "IC220501A",
        "IC220513A",
        "IC220624A",
        "IC220822A",
        "IC220907A",
        "IC221216A",
        "IC221223A",
        "IC230112A",
    ]

    for i, neutrino in enumerate(neutrinos):
        logger.info(f"Processing {neutrino} ({i+1} of {len(neutrinos)})")

        try:
            nu = NeutrinoScanner(neutrino)
            gcn_no = find_gcn_no(base_nu_name=neutrino)
            gcn_info = parse_gcn_circular(gcn_no)
            RA = gcn_info["ra"]
            Dec = gcn_info["dec"]
            RA_max = nu.ra_max
            RA_min = nu.ra_min
            Dec_max = nu.dec_max
            Dec_min = nu.dec_min

            Dec_0 = np.mean([Dec_max, Dec_min])

            RA_width = RA_max - RA_min
            Dec_width = Dec_max - Dec_min
            correction = np.cos(np.radians(Dec_0))

            rectangular_area = RA_width * Dec_width * correction

            nu.plot_overlap_with_observations(first_det_window_days=4)

            corrected_area = nu.healpix_area
            covered_prob = nu.overlap_prob
            uncorrected_area = corrected_area / (covered_prob / 100)

            logger.info(f"Rectangular area = {rectangular_area}")
            logger.info(f"Observed area (no chipgap corr.) = {uncorrected_area}")
            logger.info(f"Observed area (corrected for chipgaps) = {corrected_area}")
            logger.info(f"Covered percentage = {covered_prob}")

            good_neutrinos.append(neutrino)
            rectangular_areas.append(rectangular_area)
            corrected_areas.append(corrected_area)
            uncorrected_areas.append(uncorrected_area)

        except:
            bad_neutrinos.append(neutrino)
            rectangular_areas.append(None)
            corrected_areas.append(None)
            uncorrected_areas.append(None)

    if bad_neutrinos:
        logger.info(f"Parsing or other error, please recheck {bad_neutrinos}:")

    df = pd.DataFrame()
    df["neutrino"] = neutrinos
    df["reported_area"] = rectangular_areas
    df["obs_area_uncorr"] = uncorrected_areas
    df["obs_area_corr"] = corrected_areas

    print(df)
    df.to_csv("areas.csv")

else:
    df = pd.read_csv("areas.csv").drop(columns=["Unnamed: 0"])
    print(df)
