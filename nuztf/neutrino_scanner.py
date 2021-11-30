#!/usr/bin/env python3
# coding: utf-8

import os, logging

from astropy.time import Time
import healpy as hp
import numpy as np
from tqdm import tqdm
from os import environ
from pathlib import Path
import requests

from ztfquery.io import LOCALSOURCE

from nuztf.base_scanner import BaseScanner
from nuztf.parse_nu_gcn import find_gcn_no, parse_gcn_circular, get_latest_gcn


nu_candidate_output_dir = os.path.join(LOCALSOURCE, "neutrino_candidates")

if not os.path.exists(nu_candidate_output_dir):
    os.makedirs(nu_candidate_output_dir)

NU_RUN_CONFIG = {
    "min_ndet": 1,  # Default:2
    "min_tspan": -1,  # Default 0, but that rejects everything!
    "max_tspan": 365,
    "min_rb": 0.0,
    "max_fwhm": 5.5,
    "max_elong": 1.4,
    "max_magdiff": 1.0,
    "max_nbad": 2,
    "min_sso_dist": 20,
    "min_gal_lat": -1.0,  # Default: 14
    "gaia_rs": 20,
    "gaia_pm_signif": 3,
    "gaia_plx_signif": 3,
    "gaia_veto_gmag_min": 9,
    "gaia_veto_gmag_max": 20,
    "gaia_excessnoise_sig_max": 999,
    "ps1_sgveto_rad": 1,
    "ps1_sgveto_th": 0.8,
    "ps1_confusion_rad": 3,
    "ps1_confusion_sg_tol": 0.1,
}


class NeutrinoScanner(BaseScanner):
    def __init__(
        self,
        nu_name: str = None,
        manual_args=None,
        gcn_no: int = None,
        cone_nside: int = 128,
        t_precursor: float = None,
        logger=None,
    ):

        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        self.prob_threshold = 0.9

        if manual_args is None:

            if nu_name is not None:
                gcn_no = find_gcn_no(nu_name, logger=self.logger)

            elif gcn_no is None:
                gcn_no = get_latest_gcn(logger=self.logger)

            gcn_info = parse_gcn_circular(gcn_no)

            self.logger.info(gcn_info)

            nu_name = gcn_info["name"]
            author = gcn_info["author"]
            ra = [gcn_info["ra"], gcn_info["ra_err"][0], gcn_info["ra_err"][1]]
            dec = [gcn_info["dec"], gcn_info["dec_err"][0], gcn_info["dec_err"][1]]
            nu_time = gcn_info["time"]

            if t_precursor is not None:
                nu_time = Time(nu_time.mjd - t_precursor, format="mjd")

        else:
            (nu_name, ra, dec, nu_time) = manual_args
            author = None

        self.nu_name = nu_name
        self.author = author
        self.gcn_no = gcn_no
        self.dist = None

        self.logger.info(f"Neutrino time: {nu_time}")

        self.ra_max = float(max(ra[1:]) + ra[0])
        self.ra_min = float(min(ra[1:]) + ra[0])
        self.dec_max = float(max(dec[1:]) + dec[0])
        self.dec_min = float(min(dec[1:]) + dec[0])

        self.logger.info(f"Coordinates: RA = {ra[0]} ({self.ra_min} - {self.ra_max})")
        self.logger.info(
            f"Coordinates: Dec = {dec[0]} ({self.dec_min} - {self.dec_max})"
        )

        self.summary_path = os.path.join(nu_candidate_output_dir, nu_name)

        BaseScanner.__init__(
            self,
            t_min=nu_time,
            run_config=NU_RUN_CONFIG,
            logger=logger,
            cone_nside=cone_nside,
        )
        self.prob_threshold = 0.9
        self.area = (
            (self.ra_max - self.ra_min)
            * (self.dec_max - self.dec_min)
            * abs(np.cos(np.radians(dec[0])))
        )
        self.logger.info(f"Projected Area: {self.area}")
        (
            self.map_coords,
            self.pixel_nos,
            self.nside,
            self.map_probs,
            self.data,
            self.key,
        ) = self.unpack_map()

    def get_name(self):
        """ """
        return self.nu_name

    def get_full_name(self):
        """ """
        return f"neutrino event {self.get_name()} ({self.author} et. al, GCN {self.gcn_no})"

    def get_overlap_line(self):
        """ """
        return (
            f"We covered {self.area:.1f} sq deg, corresponding to {self.overlap_prob:.1f}% of the reported localization region. "
            "This estimate accounts for chip gaps. "
        )

    def candidate_text(
        self, ztf_id: str, first_detection: float, lul_lim: float, lul_jd: float
    ):
        """ """
        fd = Time(first_detection, format="mjd")

        text = f"{ztf_id} was first detected on {fd.utc}. "

        return text

    @staticmethod
    def get_obs_line():
        """ """
        return "Each exposure was 300s with a typical depth of 21.0 mag."

    @staticmethod
    def remove_variability_line():
        """ """
        return ""

    def filter_f_no_prv(self, res: dict):
        """ """

        ztf_id = res["objectId"]

        # Positive detection
        if res["candidate"]["isdiffpos"] not in ["t", "1"]:
            self.logger.debug(f"{ztf_id}: Negative subtraction")
            return False

        try:
            if res["candidate"]["drb"] < 0.3:
                self.logger.debug(f"{ztf_id}: DRB too low")
                return False
        except (KeyError, TypeError) as e:
            pass

        # Check contour
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            self.logger.debug(f"{ztf_id}: Not in contour")
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            self.logger.debug(
                f"{ztf_id}: Does not have 2 detections separated  by >15 mins"
            )
            return False

        return True

    def filter_f_history(self, res: dict):
        """Filter based on 2 detection requirement and probability contour requirement"""

        ztf_id = res["objectId"]

        n_detections = len(
            [x for x in res["prv_candidates"] if "isdiffpos" in x.keys()]
        )

        if n_detections < 1:
            self.logger.debug(f"{ztf_id}: Has insufficient detection")
            return False

        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            self.logger.debug(f"{ztf_id}: Not in contour")
            return False

        return True

    def find_cone_coords(self):
        """ """
        cone_coords = []
        cone_ids = []

        scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))

        self.logger.info("Finding search pixels:")

        for i in tqdm(range(hp.nside2npix(self.cone_nside))):
            ra, dec = self.extract_ra_dec(self.cone_nside, i)
            ra_rad = np.radians(ra)
            dec_rad = np.radians(dec)
            if np.logical_and(
                ra > self.ra_min - scan_radius, ra < self.ra_max + scan_radius
            ):
                if np.logical_and(
                    dec > self.dec_min - scan_radius,
                    dec < self.dec_max + scan_radius,
                ):
                    cone_coords.append((ra_rad, dec_rad))
                    cone_ids.append(i)

        cone_coords = np.array(
            cone_coords, dtype=np.dtype([("ra", float), ("dec", float)])
        )

        return cone_ids, cone_coords

    def in_contour(self, ra_deg, dec_deg):

        in_ra = np.logical_and(ra_deg > self.ra_min, ra_deg < self.ra_max)
        in_dec = np.logical_and(dec_deg > self.dec_min, dec_deg < self.dec_max)

        return np.logical_and(in_ra, in_dec)

    def unpack_map(self):
        """ """
        # nside = self.cone_nside
        nside = 1024
        map_coords = []
        pixel_nos = []

        center_ra = np.radians(np.mean([self.ra_max, self.ra_min]))
        center_dec = np.radians(np.mean([self.dec_max, self.dec_min]))
        rad = (
            np.radians(max(self.ra_max - self.ra_min, self.dec_max - self.dec_min))
            / 2.0
        )

        nearish_pixels = list(
            hp.query_disc(
                nside=nside,
                vec=hp.ang2vec(np.pi / 2.0 - center_dec, center_ra),
                radius=rad,
                nest=True,
            )
        )

        for i in tqdm(nearish_pixels):
            ra, dec = self.extract_ra_dec(nside, i)
            if self.in_contour(ra, dec):
                map_coords.append((ra, dec))
                pixel_nos.append(i)

        map_probs = np.ones_like(pixel_nos, dtype=float)
        map_probs /= np.sum(map_probs)

        key = "PROB"

        data = np.zeros(hp.nside2npix(nside), dtype=np.dtype([(key, float)]))
        data[np.array(pixel_nos)] = map_probs

        return map_coords, pixel_nos, nside, map_probs, data, key


if __name__ == "__main__":

    import logging

    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    nu = NeutrinoScanner(logger=logger)
    nu.scan_cones()
