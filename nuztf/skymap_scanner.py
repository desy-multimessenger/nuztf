#!/usr/bin/env python
# coding: utf-8

import os, datetime, wget, logging, time, json, math

from collections import Counter, defaultdict
from tqdm import tqdm
import requests
from numpy.lib.recfunctions import append_fields
import lxml.etree
import fitsio
import numpy as np
import matplotlib.pyplot as plt

from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

import healpy as hp
from ztfquery.io import LOCALSOURCE

from nuztf.base_scanner import BaseScanner
from nuztf.ampel_api import (
    ampel_api_healpix,
    ampel_api_name,
    ampel_api_lightcurve,
    ampel_api_skymap,
)

BASE_GW_DIR = os.path.join(LOCALSOURCE, "GW_skymaps")
BASE_GRB_DIR = os.path.join(LOCALSOURCE, "GRB_skymaps")
GW_CANDIDATE_OUTPUT_DIR = os.path.join(LOCALSOURCE, "GW_candidates")
GRB_CANDIDATE_OUTPUT_DIR = os.path.join(LOCALSOURCE, "GRB_candidates")
GW_CANDIDATE_CACHE = os.path.join(LOCALSOURCE, "GW_cache")
GRB_CANDIDATE_CACHE = os.path.join(LOCALSOURCE, "GRB_cache")

for entry in [
    BASE_GW_DIR,
    BASE_GRB_DIR,
    GW_CANDIDATE_OUTPUT_DIR,
    GRB_CANDIDATE_OUTPUT_DIR,
    GW_CANDIDATE_CACHE,
    GRB_CANDIDATE_CACHE,
]:
    if not os.path.exists(entry):
        os.makedirs(entry)

GW_RUN_CONFIG = {
    "min_ndet": 1,  # Default:2
    "min_tspan": -1,  # Default 0, but that rejects everything!
    "max_tspan": 365,
    "min_rb": 0.3,
    "max_fwhm": 5.5,
    "max_elong": 1.4,
    "max_magdiff": 1.0,
    "max_nbad": 2,
    "min_sso_dist": 20,
    "min_gal_lat": 0.0,  # Default: 14
    "gaia_rs": 10.0,
    "gaia_pm_signif": 3,
    "gaia_plx_signif": 3,
    "gaia_veto_gmag_min": 9,
    "gaia_veto_gmag_max": 20,
    "gaia_excessnoise_sig_max": 999,
    "ps1_sgveto_rad": 1.0,
    "ps1_sgveto_th": 0.8,
    "ps1_confusion_rad": 3.0,
    "ps1_confusion_sg_tol": 0.1,
}


class RetractionError(Exception):
    """Base class for retracted event"""

    pass


class SkymapScanner(BaseScanner):
    def __init__(
        self,
        event_name: str = None,
        skymap_file: str = None,
        scan_mode: str = "grb",
        rev: int = None,
        prob_threshold: float = 0.9,
        cone_nside: int = 64,
        n_days: int = 3,
        logger=None,
    ):
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)

        if not scan_mode in ["gw", "grb"]:
            raise ValueError(f"Scan mode must be either 'gw' or 'grb'.")

        self.prob_threshold = prob_threshold
        self.n_days = n_days
        self.scan_mode = scan_mode

        if skymap_file is None and self.scan_mode == "gw":
            self.skymap_path, self.summary_path, self.event_name = self.get_gw_skymap(
                event_name=event_name, rev=rev
            )

        elif skymap_file is None and self.scan_mode == "grb":
            self.get_grb_skymap(event_name=event_name)

        else:
            basename = os.path.basename(skymap_file)

            if self.scan_mode == "gw":
                self.skymap_path = os.path.join(BASE_GW_DIR, basename)
            else:
                self.skymap_path = os.path.join(BASE_GRB_DIR, basename)

            if skymap_file[:8] == "https://" and scan_mode == "gw":
                self.logger.info(f"Downloading from: {skymap_file}")
                self.skymap_path = os.path.join(
                    BASE_GW_DIR, os.path.basename(skymap_file[7:])
                )
                wget.download(skymap_file, self.skymap_path)

            if self.scan_mode == "gw":
                self.summary_path = os.path.join(
                    GW_CANDIDATE_OUTPUT_DIR,
                    f"{os.path.basename(skymap_file)}_{self.prob_threshold}",
                )

            else:
                self.summary_path = os.path.join(
                    GRB_CANDIDATE_OUTPUT_DIR,
                    f"{os.path.basename(skymap_file)}_{self.prob_threshold}",
                )

            self.event_name = os.path.basename(skymap_file[7:])

        self.data, t_obs, self.hpm, self.key, self.dist, self.dist_unc = self.read_map()

        t_min = Time(t_obs, format="isot", scale="utc")

        self.logger.info(f"Event time: {t_min}")
        self.logger.info("Reading map")

        self.pixel_threshold = self.find_pixel_threshold(self.data[self.key])
        (
            self.map_coords,
            self.pixel_nos,
            self.map_probs,
            self.nside,
            self.pixel_area,
        ) = self.unpack_skymap()

        BaseScanner.__init__(
            self,
            run_config=GW_RUN_CONFIG,
            t_min=t_min,
            logger=logger,
            cone_nside=cone_nside,
        )

        # By default, accept things detected within 72 hours of event time
        self.default_t_max = Time(self.t_min.jd + self.n_days, format="jd")

        self.logger.info(f"Time-range is {self.t_min} -- {self.default_t_max.isot}")

        if self.scan_mode == "gw":
            self.cache_dir = os.path.join(GW_CANDIDATE_CACHE, self.event_name)
        else:
            self.cache_dir = os.path.join(GRB_CANDIDATE_CACHE, self.event_name)

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

    def get_alerts(self):
        """Scan the skymap area and get ZTF transients"""
        self.logger.info("Commencing skymap scan")

        self.logger.debug(
            f"API skymap search: nside = {self.cone_nside} / # pixels = {len(self.cone_ids)} / timespan = {self.default_t_max.jd-self.t_min.jd:.1f} days."
        )

        time_healpix_start = time.time()

        self.queue = []

        resume = True
        chunk_size = 8000
        resume_token = None

        while resume:
            query_res, resume_token = ampel_api_skymap(
                pixels=self.cone_ids,
                nside=self.cone_nside,
                t_min_jd=self.t_min.jd,
                t_max_jd=self.default_t_max.jd,
                logger=self.logger,
                chunk_size=chunk_size,
                resume_token=resume_token,
                warn_exceeding_chunk=False,
            )
            self.queue.extend(query_res)

            if len(query_res) < chunk_size:
                resume = False
                self.logger.info("Done.")
            else:
                self.logger.info(
                    f"Chunk size reached ({chunk_size}), commencing next query."
                )

        time_healpix_end = time.time()
        time_healpix = time_healpix_end - time_healpix_start

        cache_file = os.path.join(self.cache_dir, f"{self.event_name}_all_alerts.json")

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
        cache_file = os.path.join(self.cache_dir, f"{self.event_name}_all_alerts.json")

        if load_cachefile:
            self.queue = json.load(open(cache_file, "r"))

        first_stage_objects = []
        filter_time_start = time.time()

        i_survived = []

        for i, res in enumerate(tqdm(self.queue)):

            ztf_id = res["objectId"]

            if self.filter_f_no_prv(
                res=res,
                t_min_jd=self.t_min.jd,
                t_max_jd=self.default_t_max.jd,
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

                if self.filter_f_history(
                    res=res, t_min_jd=self.t_min.jd, t_max_jd=self.default_t_max.jd
                ):
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

    def get_overlap_line(self):
        """ """
        return (
            "We covered {0:.1f}% of the enclosed probability "
            "based on the map in {1:.1f} sq deg. "
            "This estimate accounts for chip gaps. ".format(
                self.overlap_prob, self.area
            )
        )

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

    def filter_f_no_prv(self, res: dict, t_min_jd: float, t_max_jd: float) -> bool:
        """First filtering stage"""

        # Veto transients older than t_min_jd
        # as we don't expect detections before GRB or GW event time)
        if res["candidate"]["jdstarthist"] < t_min_jd:
            startdate_jd = res["candidate"]["jdstarthist"]
            startdate_date = Time(startdate_jd, format="jd").isot
            self.logger.debug(
                f"{res['objectId']}: Transient is too old (jdstarthist predates event; first detection at {startdate_date})."
            )
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > t_max_jd:
            startdate_jd = res["candidate"]["jdstarthist"]
            startdate_date = Time(startdate_jd, format="jd").isot
            self.logger.debug(
                f"{res['objectId']}: Transient is too new (jdstarthist too late after event; first detection at {startdate_date})"
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
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            self.logger.debug(f"{res['objectId']}: Outside of event contour.")
            return False

        self.logger.debug(f"{res['objectId']}: Passed first filter stage (no prv)")

        return True

    def filter_f_history(self, res: dict, t_min_jd, t_max_jd):
        """Veto transients"""

        # Veto old transients
        ztf_id = res["objectId"]
        if res["candidate"]["jdstarthist"] < t_min_jd:
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
            if np.logical_and("isdiffpos" in x.keys(), x["jd"] > t_min_jd)
        ]

        pos_detections = [x for x in old_detections if "isdiffpos" in x.keys()]

        if len(pos_detections) < 1:
            self.logger.debug(f"{ztf_id}: Does not have two detections")
            return False

        self.logger.debug(f"{ztf_id}: Passed the history filtering stage")

        return True

    def get_gw_skymap(self, event_name: str, rev: int):
        """ """
        from ligo.gracedb.rest import GraceDb

        ligo_client = GraceDb()

        self.logger.info("Obtaining skymap from GraceDB")

        if event_name is None:
            superevent_iterator = ligo_client.superevents("category: Production")
            superevent_ids = [
                superevent["superevent_id"] for superevent in superevent_iterator
            ]
            event_name = superevent_ids[0]

        voevents = ligo_client.voevents(event_name).json()["voevents"]

        if rev is None:
            rev = len(voevents)

        elif rev > len(voevents):
            raise Exception("Revision {0} not found".format(rev))

        latest_voevent = voevents[rev - 1]
        self.logger.info(f"Found voevent {latest_voevent['filename']}")

        if "Retraction" in latest_voevent["filename"]:
            raise RetractionError(
                f"The specified LIGO event, {latest_voevent['filename']}, was retracted."
            )

        response = requests.get(latest_voevent["links"]["file"])

        root = lxml.etree.fromstring(response.content)
        params = {
            elem.attrib["name"]: elem.attrib["value"]
            for elem in root.iterfind(".//Param")
        }

        latest_skymap = params["skymap_fits"]

        self.logger.info(f"Latest skymap URL: {latest_skymap}")

        base_file_name = os.path.basename(latest_skymap)
        savepath = os.path.join(
            BASE_GW_DIR,
            f"{event_name}_{latest_voevent['N']}_{base_file_name}",
        )

        self.logger.info(f"Saving to: {savepath}")
        response = requests.get(latest_skymap)

        with open(savepath, "wb") as f:
            f.write(response.content)

        summary_path = f"{GW_CANDIDATE_OUTPUT_DIR}/{event_name}_{latest_voevent['N']}_{self.prob_threshold}"

        return savepath, summary_path, event_name

    def get_grb_skymap(self, event_name: str):
        """ """
        if event_name is None:
            raise ValueError(
                "event_name must be provided for GRBs. They must have the form 'GRB210729A"
            )

        event_year_short = event_name[3:5]
        event_year = "20" + event_year_short
        event_month = event_name[5:7]
        event_day = event_name[7:9]
        event_letter = event_name[9]
        event_number = ord(event_letter) - 65

        event_date = datetime.datetime(
            int(event_year), int(event_month), int(event_day)
        )

        # get possible skymap URLs

        url = f"https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/{event_date.year}"

        from lxml import html

        page_overview = requests.get(url)
        webpage_overview = html.fromstring(page_overview.content)

        links_overview = webpage_overview.xpath("//a/@href")

        links_for_date = []

        for link in links_overview:
            if link[2:8] == f"{event_year_short}{event_month}{event_day}":
                links_for_date.append(url + "/" + link + "current/")

        if len(links_for_date) > 1:
            self.logger.info(
                f"Found multiple events. Will choose the one corresponding the GRB letter {event_letter}"
            )

        event_url = links_for_date[event_number]

        page_event = requests.get(event_url)
        webpage_event = html.fromstring(page_event.content)
        links_event = webpage_event.xpath("//a/@href")

        for link in links_event:
            if link[0:11] == "glg_healpix":
                final_link = event_url + link
                break

        self.skymap_path = os.path.join(BASE_GRB_DIR, link)

        if os.path.isfile(self.skymap_path):
            self.logger.info(
                f"Continuing with saved skymap. Located at {self.skymap_path}"
            )
        else:
            self.logger.info(f"Downloading skymap and saving to {self.skymap_path}")
            wget.download(final_link, self.skymap_path)

        self.summary_path = (
            f"{GRB_CANDIDATE_OUTPUT_DIR}/{event_name}_{self.prob_threshold}"
        )

        self.event_name = event_name

    def read_map(
        self,
    ):
        """Read the skymap"""

        self.logger.info(f"Reading file: {self.skymap_path}")

        data, h = fitsio.read(self.skymap_path, header=True)

        if "DISTMEAN" not in h:
            dist = None
        else:
            dist = h["DISTMEAN"]
        if "DISTSTD" not in h:
            dist_unc = None
        else:
            dist_unc = h["DISTSTD"]
        if "DATE-OBS" not in h:
            t_obs = fitsio.read_header(self.skymap_path)["DATE-OBS"]
        else:
            t_obs = h["DATE-OBS"]

        if "PROB" in data.dtype.names:
            key = "PROB"
        elif "PROBABILITY" in data.dtype.names:
            key = "PROB"
            prob = np.array(data["PROBABILITY"]).flatten()
            data = append_fields(data, "PROB", prob)
        else:
            raise Exception(
                "No recognised probability key in map. This is probably a weird one, right?"
            )

        if not isinstance(data[0], float):
            probs = np.array(data["PROB"]).flatten()
            data = np.array(probs, dtype=np.dtype([("PROB", float)]))

        self.logger.info(f"Summed probability is {100. * np.sum(data['PROB']):.1f}%")

        if h["ORDERING"] == "RING":
            data["PROB"] = hp.pixelfunc.reorder(data["PROB"], inp="RING", out="NESTED")
            h["ORDERING"] = "NESTED"

        hpm = HEALPix(nside=h["NSIDE"], order=h["ORDERING"], frame="icrs")

        return data, t_obs, hpm, key, dist, dist_unc

    def find_pixel_threshold(self, data):
        """ """

        ranked_pixels = np.sort(data)[::-1]
        int_sum = 0.0
        pixel_threshold = 0.0

        for i, prob in enumerate(ranked_pixels):
            int_sum += prob
            if int_sum > self.prob_threshold:
                self.logger.info(
                    f"Threshold found! \n To reach {int_sum * 100.0}% of probability, pixels with probability greater than {prob} are included."
                )
                pixel_threshold = prob
                break

        return pixel_threshold

    def unpack_skymap(self):
        """ """

        nside = hp.npix2nside(len(self.data[self.key]))

        threshold = self.find_pixel_threshold(self.data[self.key])

        mask = self.data[self.key] > threshold

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

        return map_coords, pixel_nos, self.data[self.key][mask], nside, pixel_area

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
        plt.subplot(211, projection="aitoff")

        mask = self.data[self.key] > self.pixel_threshold

        size = hp.max_pixrad(self.nside, degrees=True) ** 2

        ra_map_rad = np.deg2rad(self.wrap_around_180(self.map_coords["ra"]))
        dec_map_rad = np.deg2rad(self.map_coords["dec"])

        ra_cone_rad = np.deg2rad(self.wrap_around_180(self.cone_coords["ra"]))
        dec_cone_rad = np.deg2rad(self.cone_coords["dec"])

        plt.scatter(
            ra_map_rad,
            dec_map_rad,
            c=self.data[self.key][mask],
            vmin=0.0,
            vmax=max(self.data[self.key]),
            s=size,
        )

        plt.title("SKYMAP")

        plt.subplot(212, projection="aitoff")

        plt.scatter(ra_cone_rad, dec_cone_rad)
        plt.title("CONE REGION")

        outpath = os.path.join(BASE_GRB_DIR, f"{self.event_name}.png")
        plt.tight_layout()

        plt.savefig(outpath, dpi=300)

        return fig

    def plot_coverage(self):
        """Plot ZTF coverage of skymap region"""
        fig, message = self.plot_overlap_with_observations(
            first_det_window_days=self.n_days
        )
        plt.tight_layout()

        outpath = os.path.join(
            GRB_CANDIDATE_OUTPUT_DIR, f"{self.event_name}_coverage.png"
        )
        plt.savefig(outpath, dpi=300)

        self.logger.info(message)

        return fig, message

    def interpolate_map(self, ra_deg, dec_deg):
        """ """
        interpol_map = self.hpm.interpolate_bilinear_skycoord(
            SkyCoord(ra_deg * u.deg, dec_deg * u.deg), self.data[self.key]
        )
        return interpol_map

    def in_contour(self, ra_deg, dec_deg):
        """ """
        return self.interpolate_map(ra_deg, dec_deg) > self.pixel_threshold


if __name__ == "__main__":

    import logging

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    scanner = SkymapScanner(logger=logger)
    scanner.plot_skymap()
    scanner.get_alerts()
    scanner.filter_alerts()
