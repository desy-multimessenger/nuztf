#!/usr/bin/env python
# coding: utf-8

import json
import logging
import os
import time

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import yaml
from astropy.time import Time
from astropy_healpix import HEALPix
from tqdm import tqdm
from ztfquery.io import LOCALSOURCE

from nuztf.ampel_api import (
    ampel_api_acknowledge_chunk,
    ampel_api_lightcurve,
    ampel_api_skymap,
    get_preprocessed_results,
)
from nuztf.base_scanner import BaseScanner
from nuztf.skymap import Skymap


class RetractionError(Exception):
    """Base class for retracted event"""

    pass


class SkymapScanner(BaseScanner):
    def __init__(
        self,
        event: str = None,
        rev: int = None,
        prob_threshold: float = 0.9,
        cone_nside: int = 64,
        n_days: float = 3.0,  # By default, accept things detected within 72 hours of event time
        custom_prefix: str = "",
        config: dict = None,
    ):
        self.logger = logging.getLogger(__name__)
        self.prob_threshold = prob_threshold
        self.n_days = n_days

        if config:
            self.config = config
        else:
            config_path = os.path.join(
                os.path.dirname(__file__), "config", "gw_run_config.yaml"
            )
            with open(config_path) as f:
                self.config = yaml.safe_load(f)

        self.skymap = Skymap(
            event=event,
            rev=rev,
            prob_threshold=prob_threshold,
            custom_prefix=custom_prefix,
        )

        self.t_min = Time(self.skymap.t_obs, format="isot", scale="utc")
        self.cache_dir = os.path.join(
            self.skymap.candidate_cache, self.skymap.event_name
        )
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

        BaseScanner.__init__(
            self,
            run_config=self.config,
            t_min=self.t_min,
            cone_nside=cone_nside,
        )

        self.default_t_max = Time(self.t_min.jd + self.n_days, format="jd")
        self.summary_path = self.skymap.summary_path
        self.logger.info(f"Time-range is {self.t_min} -- {self.default_t_max.isot}")

    def get_full_name(self):
        if self.skymap.event_name is not None:
            return self.skymap.event_name
        else:
            return "?????"

    def download_results(self):
        """
        Retrieve computed results from the DESY cloud
        """
        self.logger.info("Retrieving results from the DESY cloud")
        file_basename = f"{self.skymap.event_name}_{self.skymap.rev}"

        res = get_preprocessed_results(file_basename=file_basename)

        if res is None:
            final_objects = []
        else:
            final_objects = [alert["objectId"] for alert in res]
            for alert in res:
                self.cache[alert["objectId"]] = alert

        final_objects = self.remove_duplicates(final_objects)

        self.logger.info(
            f"Retrieved {len(final_objects)} final objects for event {self.skymap.event_name} / map revision {self.skymap.rev} from DESY cloud."
        )

        self.final_candidates = final_objects

    def get_alerts(self):
        """Scan the skymap area and get ZTF transients"""
        self.logger.info("Commencing skymap scan")

        self.logger.debug(
            f"API skymap search: nside = {self.cone_nside} / # pixels = {len(self.cone_ids)} / timespan = {self.default_t_max.jd-self.t_min.jd:.1f} days."
        )

        time_healpix_start = time.time()

        self.queue = []

        resume = True
        chunk_size = 2000
        resume_token = None

        i = 0
        total_chunks = 0
        t0 = time.time()

        while resume:
            res, resume_token, chunk_id, remaining_chunks = ampel_api_skymap(
                pixels=self.cone_ids,
                nside=self.cone_nside,
                t_min_jd=self.t_min.jd,
                t_max_jd=self.default_t_max.jd,
                logger=self.logger,
                chunk_size=chunk_size,
                resume_token=resume_token,
                warn_exceeding_chunk=False,
            )
            self.queue.extend(res)

            ampel_api_acknowledge_chunk(resume_token=resume_token, chunk_id=chunk_id)

            if i == 0:
                total_chunks = remaining_chunks + 1
                self.logger.info(f"Total chunks: {total_chunks}")

            if remaining_chunks % 50 == 0 and remaining_chunks != 0:
                t1 = time.time()
                processed_chunks = total_chunks - remaining_chunks
                time_per_chunk = (t1 - t0) / processed_chunks
                remaining_time = time_per_chunk * remaining_chunks
                self.logger.info(
                    f"Remaining chunks: {remaining_chunks}. Estimated time to finish: {remaining_time/60:.0f} min"
                )

            if len(res) < chunk_size:
                resume = False
                self.logger.info("Done.")
            else:
                self.logger.debug(
                    f"Chunk size reached ({chunk_size}), commencing next query."
                )
            i += 1

        time_healpix_end = time.time()
        time_healpix = time_healpix_end - time_healpix_start

        cache_file = os.path.join(
            self.cache_dir, f"{self.skymap.event_name}_all_alerts.json"
        )

        outfile = open(cache_file, "w")
        json.dump(self.queue, outfile)
        outfile.close()

        self.n_alerts = len(self.queue)

        self.logger.info(
            f"Added {self.n_alerts} alerts found between {self.t_min} and {self.default_t_max.isot}"
        )
        self.logger.info(f"This took {time_healpix:.1f} s in total")

    def filter_alerts(self, load_cachefile=False):
        """ """
        self.logger.info(f"Commencing first stage filtering.")
        cache_file = os.path.join(
            self.cache_dir, f"{self.skymap.event_name}_all_alerts.json"
        )

        if load_cachefile:
            self.queue = json.load(open(cache_file, "r"))

        first_stage_objects = []
        filter_time_start = time.time()

        i_survived = []

        for i, res in enumerate(tqdm(self.queue)):
            ztf_id = res["objectId"]

            if self.filter_f_no_prv(
                res=res,
            ):
                self.logger.debug(
                    f"{ztf_id}: Passed first cut (does not have previous detections)."
                )
                if self.filter_ampel(res):
                    self.logger.debug(f"{ztf_id}: Passed AMPEL cut.")
                    i_survived.append(i)
                else:
                    self.logger.debug(f"{ztf_id}: Failed AMPEL cut.")
            else:
                self.logger.debug(
                    f"{ztf_id}: Failed first cut (has previous detections)."
                )

        first_stage_objects = [self.queue[i]["objectId"] for i in i_survived]
        first_stage_objects = self.remove_duplicates(first_stage_objects)

        filter_time_end = time.time()
        filter_time = filter_time_end - filter_time_start

        self.logger.info(
            f"First stage of filtering (based on predetections plus AMPEL cuts) took {filter_time:.1f} s in total. {len(first_stage_objects)} transients make the cut."
        )

        cache_file_first_stage = cache_file[:-15] + "_first_stage.json"

        outfile = open(cache_file_first_stage, "w")
        json.dump(first_stage_objects, outfile)
        outfile.close()

        # Second and final stage
        self.logger.info(
            f"Second stage commencing: Now we do additional filtering based on history."
        )

        start_secondfilter = time.time()

        final_objects = []

        for ztf_id in tqdm(first_stage_objects):
            # Get the full lightcurve from the API
            query_res = ampel_api_lightcurve(ztf_name=ztf_id, logger=self.logger)

            for res in query_res:
                _ztf_id = res["objectId"]

                if self.filter_f_history(res=res):
                    final_objects.append(_ztf_id)
                    self.cache[_ztf_id] = res
                    self.logger.debug(f"{_ztf_id}: Passed all filters.")
                else:
                    self.logger.debug(f"{_ztf_id}: Failed History.")

        end_secondfilter = time.time()
        filter_time = end_secondfilter - start_secondfilter

        final_objects = self.remove_duplicates(final_objects)

        cache_file_final_stage = cache_file[:-15] + "_final_stage.json"

        outfile = open(cache_file_final_stage, "w")
        json.dump(final_objects, outfile)
        outfile.close()

        self.logger.info(
            f"Final stage of filtering took {filter_time:.1f} s in total. {len(final_objects)} transients make the cut."
        )

        self.final_candidates = final_objects

    def remove_duplicates(self, ztf_ids: list):
        """ """
        return list(set(ztf_ids))

    def get_obs_line(self):
        """ """
        return "Each exposure was 30s with a typical depth of 20.5 mag."

    @staticmethod
    def remove_variability_line():
        """ """
        return (
            ", and removing candidates with history of "
            "variability prior to the merger time"
        )

    def candidate_text(
        self, name: str, first_detection: float, lul_lim: float, lul_jd: float
    ):
        """ """
        try:
            text = (
                "{0}, first detected {1:.1f} hours after merger, "
                "was not detected {2:.1f} days prior to a depth of {3:.2f}. ".format(
                    name,
                    24.0 * (first_detection - self.t_min.jd),
                    first_detection - lul_jd,
                    lul_lim,
                )
            )
        except TypeError:
            text = (
                f"{name} had upper limit problems. PLEASE FILL IN NUMBERS BY HAND!!! "
            )

        return text

    def filter_f_no_prv(self, res: dict, t_max_jd=None) -> bool:
        """First filtering stage"""

        if t_max_jd is None:
            t_max_jd = self.default_t_max.jd

        # Veto transients older than t_min_jd
        # as we don't expect detections before GRB or GW event time)
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            startdate_jd = res["candidate"]["jdstarthist"]
            startdate_date = Time(startdate_jd, format="jd").isot
            self.logger.debug(
                f"{res['objectId']}: Transient is too old "
                f"(jdstarthist predates event; first detection at {startdate_date})."
            )
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > t_max_jd:
            startdate_jd = res["candidate"]["jdstarthist"]
            startdate_date = Time(startdate_jd, format="jd").isot
            self.logger.debug(
                f"{res['objectId']}: Transient is too new "
                f"(jdstarthist too late after event; "
                f"first detection at {startdate_date}, "
                f"filter searches only up to {Time(t_max_jd, format='jd').isot})."
            )
            return False

        # Exclude negative detection
        if res["candidate"]["isdiffpos"] not in ["t", "1"]:
            self.logger.debug(f"{res['objectId']}: Negative subtraction")
            return False

        try:
            if res["candidate"]["drb"] < 0.3:
                self.logger.debug(f"{res['objectId']}: DRB too low")
                return False
        except (KeyError, TypeError):
            pass

        # Check contour
        if not self.skymap.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            self.logger.debug(f"{res['objectId']}: Outside of event contour.")
            return False

        self.logger.debug(f"{res['objectId']}: Passed first filter stage (no prv)")

        return True

    def filter_f_history(self, res: dict, t_max_jd=None):
        """Veto transients"""

        if t_max_jd is None:
            t_max_jd = self.default_t_max.jd

        # Veto old transients
        ztf_id = res["objectId"]
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            self.logger.debug(
                f"{ztf_id}: Transient is too old. (jdstarthist history predates event)"
            )
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > t_max_jd:
            self.logger.debug(
                f"{ztf_id}: Transient is too new. (jdstarthist too late after event)"
            )
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            self.logger.debug(f"{ztf_id}: Not passed mover cut")
            return False

        # Require 2 positive detections
        old_detections = [
            x
            for x in res["prv_candidates"]
            if np.logical_and("isdiffpos" in x.keys(), x["jd"] > self.t_min.jd)
        ]

        pos_detections = [x for x in old_detections if "isdiffpos" in x.keys()]

        if len(pos_detections) < 1:
            self.logger.debug(f"{ztf_id}: Does not have two detections")
            return False

        self.logger.debug(f"{ztf_id}: Passed the history filtering stage")

        return True

    def unpack_skymap(self):
        """ """

        nside = hp.npix2nside(len(self.skymap.data[self.skymap.key]))

        mask = self.skymap.data[self.skymap.key] > self.skymap.pixel_threshold

        map_coords = []

        pixel_nos = []

        self.logger.info("Checking which pixels are within the contour:")

        for i in tqdm(range(hp.nside2npix(nside))):
            if mask[i]:
                map_coords.append(self.extract_ra_dec(nside, i))
                pixel_nos.append(i)

        pixel_area = hp.nside2pixarea(nside, degrees=True) * float(len(map_coords))

        self.logger.info(f"Total pixel area: {pixel_area} degrees")

        map_coords = np.array(
            map_coords, dtype=np.dtype([("ra", float), ("dec", float)])
        )

        return (
            map_coords,
            pixel_nos,
            nside,
            self.skymap.data[self.skymap.key][mask],
            self.skymap.data,
            pixel_area,
            self.skymap.key,
        )

    def find_cone_coords(self):
        """ """

        cone_ids = []

        for ra, dec in self.map_coords:
            cone_ids.append(self.extract_npix(self.cone_nside, ra, dec))

        cone_ids = list(set(cone_ids))

        cone_coords = []

        for i in tqdm(cone_ids):
            cone_coords.append(self.extract_ra_dec(self.cone_nside, i))

        cone_coords = np.array(
            cone_coords, dtype=np.dtype([("ra", float), ("dec", float)])
        )

        return cone_ids, cone_coords

    def plot_skymap(self):
        """ """
        fig = plt.figure()
        plt.subplot(111, projection="aitoff")

        mask = self.data[self.key] > self.skymap.pixel_threshold

        size = hp.max_pixrad(self.nside, degrees=True) ** 2

        ra_map_rad = np.deg2rad(self.wrap_around_180(self.map_coords["ra"]))
        dec_map_rad = np.deg2rad(self.map_coords["dec"])

        plt.scatter(
            ra_map_rad,
            dec_map_rad,
            c=self.data[self.key][mask],
            vmin=0.0,
            vmax=max(self.data[self.key]),
            s=size,
        )

        plt.title("SKYMAP")

        outpath = os.path.join(
            self.skymap.base_skymap_dir, f"{self.skymap.event_name}.png"
        )
        plt.tight_layout()

        plt.savefig(outpath, dpi=300)

        return fig

    def plot_coverage(self, plot_candidates: bool = True):
        """Plot ZTF coverage of skymap region"""
        fig, message = self.plot_overlap_with_observations(
            first_det_window_days=self.n_days
        )

        if plot_candidates:
            for candidate, res in self.cache.items():
                ra = np.deg2rad(
                    self.wrap_around_180(np.array([res["candidate"]["ra"]]))
                )
                dec = np.deg2rad(res["candidate"]["dec"])

                plt.scatter(
                    ra, dec, color="white", marker="*", s=50.0, edgecolor="black"
                )

        plt.tight_layout()

        outpath = os.path.join(
            self.skymap.candidate_output_dir, f"{self.skymap.event_name}_coverage.png"
        )
        plt.savefig(outpath, dpi=300)

        self.logger.info(message)

        return fig, message


if __name__ == "__main__":
    import logging

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    scanner = SkymapScanner()
    scanner.plot_skymap()
    scanner.get_alerts()
    scanner.filter_alerts()
