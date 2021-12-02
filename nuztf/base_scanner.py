#!/usr/bin/env python3
# coding: utf-8

import os, time, json, logging, datetime, logging, pickle

import backoff
import requests
import pandas
import healpy as hp
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.cosmology import FlatLambdaCDM

from ztfquery import alert
from ztfquery import fields as ztfquery_fields

from gwemopt.ztf_tiling import get_quadrant_ipix

from ampel.ztf.t0.DecentFilter import DecentFilter
from ampel.ztf.dev.DevAlertProcessor import DevAlertProcessor
from ampel.alert.PhotoAlert import PhotoAlert

from nuztf.ampel_api import (
    ampel_api_cone,
    ampel_api_timerange,
    ampel_api_name,
    ampel_api_lightcurve,
    ampel_api_skymap,
    ensure_cutouts,
)
from nuztf.cat_match import get_cross_match_info, ampel_api_tns, query_ned_for_z
from nuztf.observation_log import get_obs_summary
from nuztf.plot import lightcurve_from_alert

DEBUG = False
RATELIMIT_CALLS = 10
RATELIMIT_PERIOD = 1

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


class BaseScanner:
    def __init__(
        self,
        run_config,
        t_min,
        resource=None,
        filter_class=DecentFilter,
        cone_nside=64,
        cones_to_scan=None,
        logger=None,
    ):
        self.cone_nside = cone_nside
        self.t_min = t_min

        if not hasattr(self, "prob_threshold"):
            self.prob_threshold = None

        if resource is None:
            resource = {
                "ampel-ztf/catalogmatch": "https://ampel.zeuthen.desy.de/api/catalogmatch/",
            }

        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        self.logger.info("AMPEL run config:")
        self.logger.info(run_config)

        lvl = self.logger.level

        if lvl > 10:
            logger_ampel = logging.getLogger("AMPEL_filter")
            logger_ampel.setLevel(logging.WARNING)
        else:
            from ampel.log.AmpelLogger import AmpelLogger

            logger_ampel = AmpelLogger()

        self.ampel_filter_class = filter_class(
            logger=logger_ampel, resource=resource, **run_config
        )

        self.dap = DevAlertProcessor(self.ampel_filter_class)

        self.scanned_pixels = []

        if cones_to_scan is None:
            self.cone_ids, self.cone_coords = self.find_cone_coords()
        else:
            self.cone_ids, self.cone_coords = cones_to_scan

        self.cache = dict()
        self.default_t_max = t_min + 10.0

        self.overlap_prob = None
        self.overlap_fields = None
        self.first_obs = None
        self.last_obs = None
        self.n_fields = None
        self.area = None
        self.double_extragalactic_area = None

        if not hasattr(self, "dist"):
            self.dist = None

    def get_name(self):
        raise NotImplementedError

    def get_full_name(self):
        raise NotImplementedError

    @staticmethod
    def get_tiling_line():
        return ""

    @staticmethod
    def get_obs_line():
        raise NotImplementedError

    @staticmethod
    def remove_variability_line():
        raise NotImplementedError

    def get_overlap_line(self):
        raise NotImplementedError

    def filter_ampel(self, res):
        return (
            self.ampel_filter_class.apply(
                PhotoAlert(res["objectId"], res["objectId"], *self.dap._shape(res))
            )
            is not None
        )

    # @sleep_and_retry
    # @limits(calls=RATELIMIT_CALLS, period=RATELIMIT_PERIOD)
    @backoff.on_exception(
        backoff.expo,
        requests.exceptions.RequestException,
        max_time=600,
    )
    def add_res_to_cache(self, res):

        for res_alert in res:

            if res_alert["objectId"] not in self.cache.keys():
                self.cache[res_alert["objectId"]] = res_alert
            elif (
                res_alert["candidate"]["jd"]
                > self.cache[res_alert["objectId"]]["candidate"]["jd"]
            ):
                self.cache[res_alert["objectId"]] = res_alert

    def add_to_cache_by_names(self, *args):
        for ztf_name in args:
            query_res = ampel_api_name(ztf_name, logger=self.logger)
            self.add_res_to_cache(query_res)

    def check_ampel_filter(self, ztf_name):
        lvl = logging.getLogger().getEffectiveLevel()
        logging.getLogger().setLevel(logging.DEBUG)
        self.logger.info("Set logger level to DEBUG")
        all_query_res = ampel_api_name(ztf_name, logger=self.logger)
        pipeline_bool = False
        for query_res in all_query_res:
            self.logger.info("Checking filter f (no prv)")
            no_prv_bool = self.filter_f_no_prv(query_res)
            self.logger.info(f"Filter f (np prv): {no_prv_bool}")
            if no_prv_bool:
                self.logger.info("Checking ampel filter")
                bool_ampel = self.filter_ampel(query_res)
                self.logger.info(f"ampel filter: {bool_ampel}")
                if bool_ampel:
                    self.logger.info("Checking filter f (history)")
                    history_bool = self.filter_f_history(query_res)
                    self.logger.info(f"Filter f (history): {history_bool}")
                    if history_bool:
                        pipeline_bool = True
        self.logger.info(f"Setting logger back to {lvl}")
        logging.getLogger().setLevel(lvl)
        return pipeline_bool

    def plot_ztf_observations(self, **kwargs):
        self.get_multi_night_summary().show_gri_fields(**kwargs)

    def get_multi_night_summary(self, max_days=None):
        return get_obs_summary(self.t_min, max_days=max_days)

    def scan_area(
        self,
        t_min=None,
        t_max=None,
    ):
        """
        Retrieve alerts for the healpix map from AMPEL API,
        filter the candidates and create a summary
        """

        if t_max is None:
            t_max = self.default_t_max

        if t_min is None:
            t_min = self.t_min

        self.logger.info("Commencing skymap scan")

        self.logger.debug(
            f"API skymap search: nside = {self.cone_nside} / # pixels = {len(self.cone_ids)} / timespan = {t_max.jd-t_min.jd:.1f} days."
        )

        query_res = []

        resume = True
        chunk_size = 8000
        resume_token = None

        while resume:
            res, resume_token = ampel_api_skymap(
                pixels=self.cone_ids,
                nside=self.cone_nside,
                t_min_jd=t_min.jd,
                t_max_jd=t_max.jd,
                logger=self.logger,
                chunk_size=chunk_size,
                resume_token=resume_token,
                warn_exceeding_chunk=False,
            )
            query_res.extend(res)

            if len(res) < chunk_size:
                resume = False
                self.logger.info("Done.")
            else:
                self.logger.info(
                    f"Chunk size reached ({chunk_size}), commencing next query."
                )

        self.logger.info(f"Before filtering: Found {len(query_res)} candidates")

        ztf_ids_first_stage = []
        for res in query_res:
            if self.filter_f_no_prv(res):
                if self.filter_ampel(res):
                    ztf_ids_first_stage.append(res["objectId"])

        ztf_ids_first_stage = list(set(ztf_ids_first_stage))

        self.logger.info(f"{len(ztf_ids_first_stage)} candidates survive filtering")

        self.logger.info(f"Retrieving alert history from AMPEL")

        results = self.ampel_object_search(ztf_ids=ztf_ids_first_stage)

        for res in results:
            self.add_res_to_cache(res)

        self.logger.info(f"Found {len(self.cache)} candidates")

        self.create_candidate_summary()

    def filter_f_no_prv(self, res):
        raise NotImplementedError

    def fast_filter_f_no_prv(self, res):
        return self.filter_f_no_prv(res)

    def filter_f_history(self, res):
        raise NotImplementedError

    def fast_filter_f_history(self, res):
        return self.filter_f_history(res)

    def find_cone_coords(self):
        raise NotImplementedError

    @staticmethod
    def wrap_around_180(ra_deg: float):
        """ """
        ra_rad = np.deg2rad(ra_deg)
        ra_rad[ra_rad > np.pi] -= 2 * np.pi
        ra_deg = np.rad2deg(ra_rad)

        return ra_deg

    def ampel_object_search(self, ztf_ids: list, with_history: bool = True) -> list:
        """ """
        all_results = []

        for ztf_id in ztf_ids:

            # get the full lightcurve from the API
            query_res = ampel_api_lightcurve(ztf_name=ztf_id, logger=self.logger)

            final_res = []

            for res in query_res:
                if self.filter_f_history(res):
                    final_res.append(res)

            all_results.append(final_res)

        return all_results

    # @staticmethod
    def calculate_abs_mag(self, mag: float, redshift: float) -> float:
        """ """
        luminosity_distance = cosmo.luminosity_distance(redshift).value * 10 ** 6
        abs_mag = mag - 5 * (np.log10(luminosity_distance) - 1)

        return abs_mag

    def parse_candidates(self):

        table = (
            "+--------------------------------------------------------------------------------+\n"
            "| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n"
            "+--------------------------------------------------------------------------------+\n"
        )
        for name, res in sorted(self.cache.items()):

            jds = [x["jd"] for x in res["prv_candidates"]]

            if res["candidate"]["jd"] > max(jds):
                latest = res["candidate"]
            else:
                latest = res["prv_candidates"][jds.index(max(jds))]

            old_flag = ""

            second_det = [x for x in jds if x > min(jds) + 0.01]

            if len(second_det) > 0:
                if Time.now().jd - second_det[0] > 1.0:
                    old_flag = "(MORE THAN ONE DAY SINCE SECOND DETECTION)"

            tns_result = " ------- "
            tns_name, tns_date, tns_group = ampel_api_tns(
                latest["ra"], latest["dec"], searchradius_arcsec=3
            )
            if tns_name:
                tns_result = tns_name

            line = "| {0} | {1} | {2:011.7f} | {3:+011.7f} | {4}      | {5:.2f} | {6:.2f}   | {7} \n".format(
                name,
                tns_result,
                float(latest["ra"]),
                float(latest["dec"]),
                ["g", "r", "i"][latest["fid"] - 1],
                latest["magpsf"],
                latest["sigmapsf"],
                old_flag,
                sign="+",
                prec=7,
            )
            table += line

        table += "+--------------------------------------------------------------------------------+\n\n"
        return table

    def draft_gcn(self):

        first_obs_dt = self.first_obs.datetime
        pretty_date = first_obs_dt.strftime("%Y-%m-%d")
        pretty_time = first_obs_dt.strftime("%H:%M")

        text = (
            f"Astronomer Name (Institute of Somewhere), ............. report,\n"
            f"On behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n"
            f"We observed the localization region of the {self.get_full_name()} with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). {self.get_tiling_line()}"
            f"We started observations in the g-band and r-band beginning at {pretty_date} {pretty_time} UTC, "
            f"approximately {(self.first_obs.jd - self.t_min.jd) * 24.0:.1f} hours after event time. {self.get_overlap_line()}"
            f"{self.get_obs_line()} \n \n"
            "The images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). "
            "AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. "
            "We reject stellar sources (Tachibana and Miller 2018) and moving objects, and "
            f"apply machine learning algorithms (Mahabal et al. 2019) {self.remove_variability_line()}. We are left with the following high-significance transient "
            "candidates by our pipeline, all lying within the "
            f"{100 * self.prob_threshold}% localization of the skymap.\n\n{self.parse_candidates()} \n\n"
        )

        if self.dist:
            text += (
                "The GW distance estimate is {:.0f} [{:.0f} - {:.0f}] Mpc.\n\n".format(
                    self.dist, self.dist - self.dist_unc, self.dist + self.dist_unc
                )
            )
        else:
            pass

        text += f"Amongst our candidates, \n{self.text_summary()}\n\n"

        text += (
            "ZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\n"
            "GROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\n"
            "Alert distribution service provided by DIRAC@UW (Patterson et al. 2019).\n"
            "Alert database searches are done by AMPEL (Nordin et al. 2019).\n"
            "Alert filtering is performed with the AMPEL Follow-up Pipeline (Stein et al. 2021).\n"
        )
        if self.dist:
            text += "Alert filtering and follow-up coordination is being undertaken by the Fritz marshal system (FIXME CITATION NEEDED)."
        return text

    @staticmethod
    def extract_ra_dec(nside, index):
        """ """
        theta, phi = hp.pix2ang(nside, index, nest=True)
        ra_deg = np.rad2deg(phi)
        dec_deg = np.rad2deg(0.5 * np.pi - theta)

        return (ra_deg, dec_deg)

    @staticmethod
    def extract_npix(nside, ra, dec):
        """ " """
        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)

        return int(hp.ang2pix(nside, theta, phi, nest=True))

    def create_candidate_summary(self, outfile=None):
        """Create pdf with lightcurve plots of all candidates"""

        if outfile is None:
            pdf_path = self.summary_path + ".pdf"
        else:
            pdf_path = outfile

        self.logger.info(f"Saving to: {pdf_path}")

        with PdfPages(pdf_path) as pdf:
            for (name, alert) in tqdm(sorted(self.cache.items())):

                fig, _ = lightcurve_from_alert(
                    [alert], include_cutouts=True, logger=self.logger
                )
                pdf.savefig()
                plt.close()

    def create_overview_table(self, outfile=None):
        """Create csv table of all candidates"""
        if outfile is None:
            csv_path = self.summary_path + ".csv"
        else:
            csv_path = outfile

        self.logger.info(f"Saving to {csv_path}")

        ztf_ids = []
        ras = []
        decs = []
        mags = []
        crossmatches = []

        data = {"ztf_id": [], "RA": [], "Dec": [], "mag": [], "xmatch": []}

        for (ztf_id, alert) in tqdm(sorted(self.cache.items())):
            data["ztf_id"].append(ztf_id)
            data["RA"].append(alert["candidate"]["ra"])
            data["Dec"].append(alert["candidate"]["dec"])
            data["mag"].append(alert["candidate"]["magpsf"])
            data["xmatch"].append(get_cross_match_info(raw=alert, logger=self.logger))

        df = pandas.DataFrame.from_dict(data)

        df.to_csv(csv_path)

    @staticmethod
    def parse_ztf_filter(fid: int):
        """ """
        return ["g", "r", "i"][fid - 1]

    def tns_summary(self):
        """ """
        for name, res in sorted(self.cache.items()):
            detections = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ]
            detection_jds = [x["jd"] for x in detections]
            first_detection = detections[detection_jds.index(min(detection_jds))]
            latest = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ][-1]
            print(
                "Candidate:",
                name,
                res["candidate"]["ra"],
                res["candidate"]["dec"],
                first_detection["jd"],
            )
            try:
                last_upper_limit = [
                    x
                    for x in res["prv_candidates"]
                    if np.logical_and(
                        "isdiffpos" in x.keys(), x["jd"] < first_detection["jd"]
                    )
                ][-1]
                print(
                    "Last Upper Limit:",
                    last_upper_limit["jd"],
                    self.parse_ztf_filter(last_upper_limit["fid"]),
                    last_upper_limit["diffmaglim"],
                )
            except IndexError:
                last_upper_limit = None
                print("Last Upper Limit: None")
            print(
                "First Detection:",
                first_detection["jd"],
                self.parse_ztf_filter(first_detection["fid"]),
                first_detection["magpsf"],
                first_detection["sigmapsf"],
            )
            hours_after_merger = 24.0 * (first_detection["jd"] - self.t_min.jd)
            print(f"First observed {hours_after_merger} hours after merger")
            if last_upper_limit:
                print(
                    "It has risen",
                    -latest["magpsf"] + last_upper_limit["diffmaglim"],
                    self.parse_ztf_filter(latest["fid"]),
                    self.parse_ztf_filter(last_upper_limit["fid"]),
                )
            print(
                [
                    x["jd"]
                    for x in res["prv_candidates"] + [res["candidate"]]
                    if "isdiffpos" in x.keys()
                ]
            )
            print("\n")

    def peak_mag_summary(self):
        for name, res in sorted(self.cache.items()):

            detections = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ]
            detection_mags = [x["magpsf"] for x in detections]
            brightest = detections[detection_mags.index(min(detection_mags))]

            diff = 0.0
            df = None

            for fid in [1, 2, 3]:
                dets = [x["magpsf"] for x in detections if int(x["fid"]) == fid]
                if len(dets) > 1:
                    nd = max(dets) - min(dets)

                    if nd > diff:
                        diff = nd
                        df = self.parse_ztf_filter(fid)

            tns_result = ""
            tns_name, tns_date, tns_group = ampel_api_tns(
                brightest["ra"], brightest["dec"], searchradius_arcsec=3.0
            )
            if tns_name:
                tns_result = f"({tns_name})"

            xmatch_info = get_cross_match_info(raw=res, logger=self.logger)

            print(
                f"Candidate {name} peaked at {brightest['magpsf']:.1f} {tns_result}on "
                f"{brightest['jd']:.1f} with filter {self.parse_ztf_filter(brightest['fid'])}. "
                f"Max range of {diff:.1f} mag with filter {df}. {xmatch_info}"
            )

    def candidate_text(self, name, first_detection, lul_lim, lul_jd):
        raise NotImplementedError

    def text_summary(self):
        """ """
        text = ""
        for name, res in sorted(self.cache.items()):
            detections = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ]
            detection_jds = [x["jd"] for x in detections]
            first_detection = detections[detection_jds.index(min(detection_jds))]
            latest = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ][-1]
            try:
                last_upper_limit = [
                    x
                    for x in res["prv_candidates"]
                    if np.logical_and(
                        "isdiffpos" in x.keys(), x["jd"] < first_detection["jd"]
                    )
                ][-1]

                text += self.candidate_text(
                    name,
                    first_detection["jd"],
                    last_upper_limit["diffmaglim"],
                    last_upper_limit["jd"],
                )

            # No pre-detection upper limit
            except IndexError:
                text += self.candidate_text(name, first_detection["jd"], None, None)

            ned_z, ned_dist = query_ned_for_z(
                ra_deg=latest["ra"],
                dec_deg=latest["dec"],
                searchradius_arcsec=20,
                logger=self.logger,
            )

            if ned_z:
                ned_z = float(ned_z)
                absmag = self.calculate_abs_mag(latest["magpsf"], ned_z)
                if ned_z > 0:
                    z_dist = Distance(z=ned_z, cosmology=cosmo).value
                    text += f"It has a spec-z of {ned_z:.3f} [{z_dist:.0f} Mpc] and an abs. mag of {absmag:.1f}. Distance to SDSS galaxy is {ned_dist:.2f} arcsec. "
                    if self.dist:
                        gw_dist_interval = [
                            self.dist - self.dist_unc,
                            self.dist + self.dist_unc,
                        ]

            c = SkyCoord(res["candidate"]["ra"], res["candidate"]["dec"], unit="deg")
            g_lat = c.galactic.b.degree
            if abs(g_lat) < 15.0:
                text += f"It is located at a galactic latitude of {g_lat:.2f} degrees. "

            xmatch_info = get_cross_match_info(raw=res, logger=self.logger)
            text += xmatch_info
            text += "\n"

        return text

    def plot_overlap_with_observations(
        self, fields=None, pid=None, first_det_window_days=None, min_sep=0.01
    ):
        """ """

        fig = plt.figure()
        plt.subplot(projection="aitoff")

        double_in_plane_probs = []
        single_in_plane_prob = []

        if fields is None:
            mns = self.get_multi_night_summary(first_det_window_days)

        else:

            class MNS:
                def __init__(self, data):
                    self.data = pandas.DataFrame(
                        data, columns=["field", "ra", "dec", "datetime"]
                    )

            data = []

            for f in fields:
                ra, dec = ztfquery_fields.get_field_centroid(f)[0]
                for i in range(2):
                    t = Time(self.t_min.jd + 0.1 * i, format="jd").utc
                    t.format = "isot"
                    t = t.value
                    data.append([f, ra, dec, t])

            mns = MNS(data)

        # Skip all 64 simultaneous quadrant entries, we only need one per observation for qa log
        # data = mns.data.copy().iloc[::64]

        data = mns.data.copy()

        ras = np.ones_like(data["field"]) * np.nan
        decs = np.ones_like(data["field"]) * np.nan

        # Actually load up ra/dec

        veto_fields = []

        for field in list(set(data["field"])):

            mask = data["field"] == field

            res = ztfquery_fields.get_field_centroid(field)

            if len(res) > 0:

                ras[mask] = res[0][0]
                decs[mask] = res[0][1]

            else:
                veto_fields.append(field)

        if len(veto_fields) > 0:
            self.logger.info(
                f"No RA/Dec found by ztfquery for fields {veto_fields}. "
                f"These observation have to be ignored."
            )

        data["ra"] = ras
        data["dec"] = decs

        mask = np.array([~np.isnan(x) for x in data["ra"]])

        data = data[mask]

        if pid is not None:
            pid_mask = data["pid"] == str(pid)
            data = data[pid_mask]

        obs_times = np.array(
            [
                Time(
                    data["datetime"].iat[i].replace(" ", "T"),
                    format="isot",
                    scale="utc",
                )
                for i in range(len(data))
            ]
        )

        if first_det_window_days is not None:
            first_det_mask = [
                x < Time(self.t_min.jd + first_det_window_days, format="jd").utc
                for x in obs_times
            ]
            data = data[first_det_mask]
            obs_times = obs_times[first_det_mask]

        self.logger.info(f"Most recent observation found is {obs_times[-1]}")
        self.logger.info("Unpacking observations")

        pix_map = dict()
        pix_obs_times = dict()

        infile = os.path.join(
            "nuztf", "data", f"ztf_fields_ipix_nside={self.nside}.pickle"
        )

        # Generate a lookup table for field healpix
        # if none exists (because this is computationally costly)
        if not os.path.isfile(infile):
            self.generate_flatpix_file()

        with open(infile, "rb") as f:
            field_pix = pickle.load(f)
        f.close()

        for i, obs_time in enumerate(tqdm(obs_times)):

            field = data["field"].iat[i]
            ra = data["ra"].iat[i]
            dec = data["dec"].iat[i]

            flat_pix = field_pix[field]

            t = obs_time.jd

            for p in flat_pix:
                if p not in pix_obs_times.keys():
                    pix_obs_times[p] = [t]
                else:
                    pix_obs_times[p] += [t]

                if p not in pix_map.keys():
                    pix_map[p] = [field]
                else:
                    pix_map[p] += [field]

        npix = hp.nside2npix(self.nside)
        theta, phi = hp.pix2ang(self.nside, np.arange(npix), nest=False)
        radecs = SkyCoord(ra=phi * u.rad, dec=(0.5 * np.pi - theta) * u.rad)
        idx = np.where(np.abs(radecs.galactic.b.deg) <= 10.0)[0]

        double_in_plane_pixels = []
        double_in_plane_probs = []
        single_in_plane_pixels = []
        single_in_plane_prob = []
        veto_pixels = []
        plane_pixels = []
        plane_probs = []
        times = []
        double_no_plane_prob = []
        double_no_plane_pixels = []
        single_no_plane_prob = []
        single_no_plane_pixels = []

        overlapping_fields = []
        for i, p in enumerate(tqdm(hp.nest2ring(self.nside, self.pixel_nos))):

            if p in pix_obs_times.keys():

                if p in idx:
                    plane_pixels.append(p)
                    plane_probs.append(self.map_probs[i])

                obs = pix_obs_times[p]

                # check which healpix are observed twice
                if max(obs) - min(obs) > min_sep:
                    # is it in galactic plane or not?
                    if p not in idx:
                        double_no_plane_prob.append(self.map_probs[i])
                        double_no_plane_pixels.append(p)
                    else:
                        double_in_plane_probs.append(self.map_probs[i])
                        double_in_plane_pixels.append(p)

                else:
                    if p not in idx:
                        single_no_plane_pixels.append(p)
                        single_no_plane_prob.append(self.map_probs[i])
                    else:
                        single_in_plane_prob.append(self.map_probs[i])
                        single_in_plane_pixels.append(p)

                overlapping_fields += pix_map[p]

                times += list(obs)
            else:
                veto_pixels.append(p)

        overlapping_fields = sorted(list(set(overlapping_fields)))
        self.overlap_fields = list(set(overlapping_fields))

        self.overlap_prob = np.sum(double_in_plane_probs + double_no_plane_prob) * 100.0

        size = hp.max_pixrad(self.nside) ** 2 * 50.0

        veto_pos = np.array(
            [hp.pixelfunc.pix2ang(self.nside, i, lonlat=True) for i in veto_pixels]
        ).T

        if len(veto_pos) > 0:

            plt.scatter(
                self.wrap_around_180(np.radians(veto_pos[0])),
                np.radians(veto_pos[1]),
                color="red",
                s=size,
            )

        plane_pos = np.array(
            [hp.pixelfunc.pix2ang(self.nside, i, lonlat=True) for i in plane_pixels]
        ).T

        if len(plane_pos) > 0:

            plt.scatter(
                self.wrap_around_180(np.radians(plane_pos[0])),
                np.radians(plane_pos[1]),
                color="green",
                s=size,
            )

        single_pos = np.array(
            [
                hp.pixelfunc.pix2ang(self.nside, i, lonlat=True)
                for i in single_no_plane_pixels
            ]
        ).T

        if len(single_pos) > 0:
            plt.scatter(
                self.wrap_around_180(np.radians(single_pos[0])),
                np.radians(single_pos[1]),
                c=single_no_plane_prob,
                vmin=0.0,
                vmax=max(self.data[self.key]),
                s=size,
                cmap="gray",
            )

        plot_pos = np.array(
            [
                hp.pixelfunc.pix2ang(self.nside, i, lonlat=True)
                for i in double_no_plane_pixels
            ]
        ).T

        if len(plot_pos) > 0:
            plt.scatter(
                self.wrap_around_180(np.radians(plot_pos[0])),
                np.radians(plot_pos[1]),
                c=double_no_plane_prob,
                vmin=0.0,
                vmax=max(self.data[self.key]),
                s=size,
            )

        red_patch = mpatches.Patch(color="red", label="Not observed")
        gray_patch = mpatches.Patch(color="gray", label="Observed once")
        violet_patch = mpatches.Patch(
            color="green", label="Observed Galactic Plane (|b|<10)"
        )
        plt.legend(handles=[red_patch, gray_patch, violet_patch])

        message = (
            "In total, {0:.2f} % of the contour was observed at least once. \n "
            "This estimate includes {1:.2f} % of the contour "
            "at a galactic latitude <10 deg. \n "
            "In total, {2:.2f} % of the contour was observed at least twice. \n"
            "In total, {3:.2f} % of the contour was observed at least twice, "
            "and excluding low galactic latitudes. \n"
            "These estimates account for chip gaps.".format(
                100
                * (
                    np.sum(double_in_plane_probs)
                    + np.sum(single_in_plane_prob)
                    + np.sum(single_no_plane_prob)
                    + np.sum(double_no_plane_prob)
                ),
                100 * np.sum(plane_probs),
                100.0 * (np.sum(double_in_plane_probs) + np.sum(double_no_plane_prob)),
                100.0 * np.sum(double_no_plane_prob),
            )
        )

        all_pix = single_in_plane_pixels + double_in_plane_pixels

        n_pixels = len(
            single_in_plane_pixels
            + double_in_plane_pixels
            + double_no_plane_pixels
            + single_no_plane_pixels
        )
        n_double = len(double_no_plane_pixels + double_in_plane_pixels)
        n_plane = len(plane_pixels)

        self.area = hp.pixelfunc.nside2pixarea(self.nside, degrees=True) * n_pixels
        self.double_extragalactic_area = (
            hp.pixelfunc.nside2pixarea(self.nside, degrees=True) * n_double
        )
        plane_area = hp.pixelfunc.nside2pixarea(self.nside, degrees=True) * n_plane

        try:
            self.first_obs = Time(min(times), format="jd")
            self.first_obs.utc.format = "isot"
            self.last_obs = Time(max(times), format="jd")
            self.last_obs.utc.format = "isot"

        except ValueError:
            raise Exception(
                f"No observations of this field were found at any time between {self.t_min} and"
                f"{obs_times[-1]}. Coverage overlap is 0%, but recent observations might be missing!"
            )

        self.logger.info(f"Observations started at {self.first_obs.jd}")

        self.overlap_fields = overlapping_fields

        #     area = (2. * base_ztf_rad)**2 * float(len(overlapping_fields))
        #     n_fields = len(overlapping_fields)

        self.logger.info(
            f"{n_pixels} pixels were covered, covering approximately {self.area:.2g} sq deg."
        )
        self.logger.info(
            f"{n_double} pixels were covered at least twice (b>10), covering approximately {self.double_extragalactic_area:.2g} sq deg."
        )
        self.logger.info(
            f"{n_plane} pixels were covered at low galactic latitude, covering approximately {plane_area:.2g} sq deg."
        )
        return fig, message

    def generate_flatpix_file(self):
        """
        Generate and save the fields-healpix lookup table
        """
        import pickle
        from ztfquery.fields import FIELD_DATAFRAME

        self.logger.info(
            f"Generating field-healpix lookup table for nside={self.nside}"
        )

        FIELD_DATAFRAME = FIELD_DATAFRAME.reset_index()

        fields = FIELD_DATAFRAME["ID"].values
        ras = FIELD_DATAFRAME["RA"].values
        decs = FIELD_DATAFRAME["Dec"].values

        flat_pix_dict = dict()

        for i, field in tqdm(enumerate(fields), total=len(fields)):

            ra = ras[i]
            dec = decs[i]
            pix = get_quadrant_ipix(self.nside, ra, dec)

            flat_pix = []

            for sub_list in pix:
                for p in sub_list:
                    flat_pix.append(p)

            flat_pix = list(set(flat_pix))
            flat_pix_dict[field] = flat_pix

        outdir = os.path.join("nuztf", "data")
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outfile = os.path.join(outdir, f"ztf_fields_ipix_nside={self.nside}.pickle")
        with open(outfile, "wb") as f:
            pickle.dump(flat_pix_dict, f)
        f.close()

    def crosscheck_prob(self):

        try:
            nside = self.ligo_nside
        except AttributeError:
            nside = self.nside

        class MNS:
            def __init__(self, data):
                self.data = pandas.DataFrame(
                    data, columns=["field", "ra", "dec", "datetime"]
                )

        data = []

        for f in self.overlap_fields:
            ra, dec = ztfquery_fields.field_to_coords(float(f))[0]
            t = Time(self.t_min.jd, format="jd").utc
            t.format = "isot"
            t = t.value
            data.append([f, ra, dec, t])

            mns = MNS(data)

        data = mns.data.copy()

        self.logger.info("Unpacking observations")
        field_prob = 0.0

        ps = []

        for index, row in tqdm(data.iterrows()):
            pix = get_quadrant_ipix(nside, row["ra"], row["dec"])

            flat_pix = []

            for sub_list in pix:
                for p in sub_list:
                    flat_pix.append(p)

            flat_pix = list(set(flat_pix))
            ps += flat_pix

        ps = list(set(ps))

        for p in hp.ring2nest(nside, ps):
            field_prob += self.data[self.key][int(p)]

        self.logger.info(
            f"Intergrating all fields overlapping 90% contour gives {100*field_prob:.2g}%"
        )
