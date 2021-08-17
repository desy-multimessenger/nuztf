#!/usr/bin/env python3
# coding: utf-8

import os
import time
import backoff
import json
import requests
import pandas
import datetime
import logging
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.cosmology import Planck18 as cosmo
from ztfquery import alert, query, skyvision, io
from ztfquery import fields as ztfquery_fields
from matplotlib.backends.backend_pdf import PdfPages
from requests.auth import HTTPBasicAuth
from tqdm import tqdm
from ampel.ztf.t0.DecentFilter import DecentFilter
from ampel.ztf.dev.DevAlertProcessor import DevAlertProcessor
from ampel.alert.PhotoAlert import PhotoAlert
from ratelimit import limits, sleep_and_retry
from gwemopt.ztf_tiling import get_quadrant_ipix
from ampel.log.AmpelLogger import AmpelLogger
from nuztf.ampel_api import ampel_api_cone, ampel_api_name, reassemble_alert, ampel_api_catalog, ampel_api_tns

DEBUG = False
RATELIMIT_CALLS = 10
RATELIMIT_PERIOD = 1

class AmpelWizard:
    def __init__(
        self,
        run_config,
        t_min,
        logger=None,
        resource=None,
        filter_class=DecentFilter,
        cone_nside=64,
        cones_to_scan=None,
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
            self.logger = AmpelLogger()
        else:
            self.logger = logger

        self.logger.info("AMPEL run config:")
        self.logger.info(run_config)

        self.ampel_filter_class = filter_class(
            logger=self.logger, resource=resource, **run_config
        )

        self.dap = DevAlertProcessor(self.ampel_filter_class)

        if not hasattr(self, "output_path"):
            self.output_path = None

        self.scanned_pixels = []
        if cones_to_scan is None:
            self.cone_ids, self.cone_coords = self.find_cone_coords()
        else:
            self.cone_ids, self.cone_coords = cones_to_scan
        self.cache = dict()
        self.default_t_max = t_min + 10.

        self.mns_time = str(self.t_min).split("T")[0].replace("-", "")
        self.mns = None

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

    def get_avro_by_name(self, ztf_name):
        return ampel_api_name(ztf_name, logger=self.logger)

    def add_to_cache_by_names(self, *args):
        for ztf_name in args:
            self.cache[ztf_name] = self.get_avro_by_name(ztf_name)

    def check_ampel_filter(self, ztf_name):
        lvl = logging.getLogger().getEffectiveLevel()
        logging.getLogger().setLevel(logging.DEBUG)
        self.logger.info("Set logger level to DEBUG")
        query_res = self.get_avro_by_name(ztf_name)
        self.logger.info("Checking filter f (no prv)")
        no_prv_bool = self.filter_f_no_prv(query_res)
        self.logger.info(f"Filter f (np prv): {no_prv_bool}")
        self.logger.info("Checking ampel filter")
        bool_ampel = self.filter_ampel(query_res)
        self.logger.info(f"ampel filter: {bool_ampel}")
        self.logger.info("Checking filter f (history)")
        history_bool = self.filter_f_history(query_res)
        self.logger.info(f"Filter f (history): {history_bool}")
        self.logger.info(f"Setting logger back to {lvl}")
        logging.getLogger().setLevel(lvl)
        return bool_ampel

    def plot_ztf_observations(self, **kwargs):
        self.get_multi_night_summary().show_gri_fields(**kwargs)

    def get_multi_night_summary(self, max_days=None):
        now = datetime.datetime.now()

        if max_days is not None:
            date_1 = datetime.datetime.strptime(self.mns_time, "%Y%m%d")
            end_date = date_1 + datetime.timedelta(days=max_days)
            end_date = end_date.strftime("%Y-%m-%d")
        else:
            end_date = now.strftime("%Y-%m-%d")

        if self.mns is None:
            start_date_jd = self.t_min.jd
            start_date_jd = Time(start_date_jd, format="jd").jd
            end_date_jd = Time(now).jd

            self.mns = skyvision.CompletedLog.from_daterange(
                self.mns_time, end=end_date, verbose=False
            )

            self.mns.data["obsjd"] = Time(
                list(self.mns.data.datetime.values), format="isot"
            ).jd

            self.mns.data.query(f"obsjd > {start_date_jd}", inplace=True)

            self.mns.data.reset_index(inplace=True)
            self.mns.data.drop(columns=["index"], inplace=True)

        return self.mns

    # def simulate_observations(self, fields):

    def scan_cones(self, t_max=None, max_cones=None):
        """ """

        if max_cones is None:
            max_cones = len(self.cone_ids)

        scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))

        self.logger.info("Commencing Ampel queries!")
        self.logger.info(f"Scan radius is {scan_radius}")
        self.logger.info(
            f"So far, {len(self.scanned_pixels)} pixels out of {len(self.cone_ids)} have already been scanned."
        )

        all_ztf_names = []

        for i, cone_id in enumerate(tqdm(list(self.cone_ids)[:max_cones])):
            ra, dec = self.cone_coords[i]

            if cone_id not in self.scanned_pixels:
                ztf_names = self.ampel_cone_search(
                    ra=np.degrees(ra),
                    dec=np.degrees(dec),
                    radius=scan_radius,
                    t_max=t_max,
                )

                if ztf_names:
                    for ztf_name in ztf_names:
                        all_ztf_names.append(ztf_name)

                self.scanned_pixels.append(cone_id)

        self.logger.info(f"Scanned {len(self.scanned_pixels)} pixels")

        # remove duplicates
        all_ztf_names = list(set(all_ztf_names))

        self.logger.info(f"Before filtering: Found {len(all_ztf_names)} candidates")
        self.logger.info(f"Retrieving alert history from AMPEL")

        results = self.ampel_object_search(
            ztf_names=all_ztf_names
        )

        for res in results:

            for res_alert in res:

                if res_alert["objectId"] not in self.cache.keys():
                    self.cache[res_alert["objectId"]] = res_alert
                elif (
                    res_alert["candidate"]["jd"]
                    > self.cache[res_alert["objectId"]]["candidate"]["jd"]
                ):
                    self.cache[res_alert["objectId"]] = res_alert

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
    def wrap_around_180(ra):
        ra[ra > np.pi] -= 2 * np.pi
        return ra

    # The backoff-stuff is to handle serverside
    # rate-limits when querying the API

    # @sleep_and_retry
    # @limits(calls=RATELIMIT_CALLS, period=RATELIMIT_PERIOD)
    @backoff.on_exception(
        backoff.expo,
        requests.exceptions.RequestException,
        max_time=600,
    )
    def ampel_cone_search(
        self, ra: float, dec: float, radius: float, t_max=None
    ) -> list:
        """ """
        if t_max is None:
            t_max = self.default_t_max

        t_min = self.t_min

        query_res = ampel_api_cone(ra, dec, radius, t_min.jd, t_max.jd, logger=self.logger)

        ## LEGACY (KEPT FOR TIMING RESULTS FOR JVS)
        # result = ampel_client.get_alerts_in_cone(
        #     ra=ra, dec=dec, radius=radius, jd_min=self.t_min.jd, jd_max=t_max.jd, with_history=False, max_blocks=100)
        # query_res = [i for i in result]

        ztf_names = []
        for res in query_res:
            if self.filter_f_no_prv(res):
                if self.filter_ampel(res):
                    ztf_names.append(res["objectId"])

        ztf_names = list(set(ztf_names))

        return ztf_names

    # @sleep_and_retry
    # @limits(calls=RATELIMIT_CALLS, period=RATELIMIT_PERIOD)

    def ampel_object_search(self, ztf_names: list) -> list:
        """ """
        all_results = []

        ## LEGACY (KEPT FOR TIMING RESULTS FOR JVS)
        # ztf_object = ampel_client.get_alerts_for_object(objectids, with_history=True)
        # query_res = [i for i in ztf_object]
        # query_res = self.merge_alerts(query_res)

        for ztf_name in ztf_names:

            query_res = ampel_api_name(ztf_name)

            final_res = []

            for res in query_res:
                if self.filter_f_history(res):
                    final_res.append(res)

            all_results.append(final_res)

        return all_results

    @staticmethod
    def calculate_abs_mag(mag, redshift: float):
        luminosity_distance = cosmo.luminosity_distance(redshift).value * 10 ** 6
        abs_mag = mag - 5 * (np.log10(luminosity_distance) - 1)
        return abs_mag

    def query_catalog(
        self,
        catalog: str,
        catalog_type: str,
        ra: float,
        dec: float,
        searchradius_arcsec: float = 10,
        searchtype: str = "all",
    ):
        """
        Method for querying catalogs via the Ampel API
        'catalog' must be the name of a supported catalog, e.g.
        SDSS_spec, PS1, NEDz_extcats...
        For a full list of catalogs, confer
        https://ampel.zeuthen.desy.de/api/catalogmatch/catalogs

        """
        assert catalog_type in ["extcats", "catsHTM"]
        assert searchtype in ["all", "nearest"]

        queryurl_catalogmatch = API_CATALOGMATCH_URL + "/cone_search/" + searchtype

        # First, we create a json body to post
        headers = {"accept": "application/json", "Content-Type": "application/json"}
        query = {
            "ra_deg": ra,
            "dec_deg": dec,
            "catalogs": [
                {"name": catalog, "rs_arcsec": searchradius_arcsec, "use": catalog_type}
            ],
        }

        response = requests.post(
            url=queryurl_catalogmatch, json=query, headers=headers
        ).json()[0]

        return response

    def query_ned_for_z(self, ra: float, dec: float, searchradius_arcsec: float = 20):

        z = None
        dist_arcsec = None

        query = ampel_api_catalog(
            catalog="NEDz_extcats",
            catalog_type="extcats",
            ra=ra,
            dec=dec,
            searchradius_arcsec=searchradius_arcsec,
            searchtype="nearest",
        )

        if query:
            z = query["body"]["z"]
            dist_arcsec = query["dist_arcsec"]
        return z, dist_arcsec

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

        text = (
            "Astronomer Name (Institute of Somewhere), ............. report,\n"
            "On behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n"
            "We observed the localization region of the {0} with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). {1}"
            "We started observations in the g-band and r-band beginning at {2} UTC, "
            "approximately {3:.1f} hours after event time. {4}"
            "{5} \n \n"
            "The images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). "
            "AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. "
            "We reject stellar sources (Tachibana and Miller 2018) and moving objects, and "
            "apply machine learning algorithms (Mahabal et al. 2019) {6}. We are left with the following high-significance transient "
            "candidates by our pipeline, all lying within the "
            "{7}% localization of the skymap. \n\n{8} \n\n".format(
                self.get_full_name(),
                self.get_tiling_line(),
                self.first_obs.utc,
                (self.first_obs.jd - self.t_min.jd) * 24.0,
                self.get_overlap_line(),
                self.get_obs_line(),
                self.remove_variability_line(),
                100 * self.prob_threshold,
                self.parse_candidates(),
            )
        )

        if self.dist:
            text += (
                "The GW distance estimate is {:.0f} [{:.0f} - {:.0f}] Mpc.\n\n".format(
                    self.dist, self.dist - self.dist_unc, self.dist + self.dist_unc
                )
            )
        else:
            pass

        text += "Amongst our candidates, \n{0}. \n \n".format(self.text_summary())

        text += (
            "Based on observations obtained with the Samuel Oschin Telescope 48-inch and the 60-inch Telescope at the Palomar Observatory as part of the Zwicky Transient Facility project. ZTF is supported by the National Science Foundation under Grant No. AST-2034437 and a collaboration including Caltech, IPAC, the Weizmann Institute for Science, the Oskar Klein Center at Stockholm University, the University of Maryland, Deutsches Elektronen-Synchrotron and Humboldt University, the TANGO Consortium of Taiwan, the University of Wisconsin at Milwaukee, Trinity College Dublin, Lawrence Livermore National Laboratories, and IN2P3, France. Operations are conducted by COO, IPAC, and UW. \n"
            "GROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949. \n"
            "Alert distribution service provided by DIRAC@UW (Patterson et al. 2019). \n"
            "Alert database searches are done by AMPEL (Nordin et al. 2019). \n"
            "Alert filtering is performed with the AMPEL Follow-up Pipeline (Stein et al. 2021). \n"
        )
        if self.dist:
            text += "Alert filtering and follow-up coordination is being undertaken by the Fritz marshal system (FIXME CITATION NEEDED)."
        return text

    @staticmethod
    def extract_ra_dec(nside, index):
        (colat, ra) = hp.pix2ang(nside, index, nest=True)
        dec = np.pi / 2.0 - colat
        return (ra, dec)

    @staticmethod
    def extract_npix(nside, ra, dec):
        colat = np.pi / 2.0 - dec
        return hp.ang2pix(nside, colat, ra, nest=True)

    def create_candidate_summary(self):

        self.logger.info(f"Saving to: {self.output_path}")

        with PdfPages(self.output_path) as pdf:
            for (name, old_alert) in tqdm(sorted(self.cache.items())):
                mock_alert = reassemble_alert(old_alert)
                try:
                    fig = alert.display_alert(mock_alert, show_ps_stamp=True)
                    fig.text(0.01, 0.01, name)
                    pdf.savefig()
                    plt.close()
                except TypeError:
                    self.logger.info(
                        f"WARNING!!! {name} will be missing from the report pdf for some reason."
                    )
                    pass

    @staticmethod
    def parse_ztf_filter(fid):
        return ["g", "r", "i"][fid - 1]

    def tns_summary(self):
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

    def candidate_text(self, name, first_detection, lul_lim, lul_jd):
        raise NotImplementedError

    def text_summary(self):
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

            ned_z, ned_dist = self.query_ned_for_z(
                latest["ra"], latest["dec"], searchradius_arcsec=20
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
            text += "\n"
        return text

    def simple_plot_overlap_with_observations(
        self, fields=None, first_det_window_days=None
    ):

        try:
            nside = self.ligo_nside
        except AttributeError:
            nside = self.nside

        fig = plt.figure()
        plt.subplot(projection="aitoff")

        probs = []
        single_probs = []

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
                ra, dec = ztfquery_fields.field_to_coords(f)[0]
                t = Time(Time.now().jd + 1.0, format="jd").utc
                t.format = "isot"
                t = t.value
                for _ in range(2):
                    data.append([f, ra, dec, t])

            mns = MNS(data)

        ras = np.degrees(
            self.wrap_around_180(
                np.array([np.radians(float(x)) for x in mns.data["ra"]])
            )
        )
        fields = list(mns.data["field"])

        plot_ras = []
        plot_decs = []

        single_ras = []
        single_decs = []

        veto_ras = []
        veto_decs = []

        overlapping_fields = []

        base_ztf_rad = 3.5
        ztf_dec_deg = 30.0

        for j, (ra, dec) in enumerate(tqdm(self.map_coords)):
            ra_deg = np.degrees(self.wrap_around_180(np.array([ra])))
            # ra_deg = self.wrap_around_180(np.array(np.degrees(ra)))
            dec_deg = np.degrees(dec)
            # (np.cos(dec - np.radians(ztf_dec_deg))
            ztf_rad = base_ztf_rad
            ztf_height = 3.7

            n_obs = 0

            for i, x in enumerate(mns.data["dec"]):
                if np.logical_and(
                    not dec_deg < float(x) - ztf_height,
                    not dec_deg > float(x) + ztf_height,
                ):
                    if abs(dec_deg - ztf_dec_deg) < 70.0:
                        if np.logical_and(
                            not ra_deg < float(ras[i]) - ztf_rad / abs(np.cos(dec)),
                            not ra_deg > float(ras[i]) + ztf_rad / abs(np.cos(dec)),
                        ):
                            n_obs += 1
                            fid = fields[i]
                            if fid not in overlapping_fields:
                                overlapping_fields.append(fields[i])

            if n_obs > 1:
                probs.append(self.map_probs[j])
                plot_ras.append(ra)
                plot_decs.append(dec)

            elif n_obs > 0:
                single_probs.append(self.map_probs[j])
                single_ras.append(ra)
                single_decs.append(dec)

            else:
                veto_ras.append(ra)
                veto_decs.append(dec)

        overlapping_fields = list(set(overlapping_fields))

        obs_times = np.array(
            [
                Time(mns.data["datetime"].iat[i], format="isot", scale="utc")
                for i in range(len(mns.data))
                if mns.data["field"].iat[i] in overlapping_fields
            ]
        )

        self.first_obs = min(obs_times)
        self.last_obs = max(obs_times)

        size = hp.max_pixrad(nside, degrees=True) ** 2

        # print(hp.max_pixrad(self.ligo_nside, degrees=True)**2 * np.pi, size)

        plt.scatter(
            self.wrap_around_180(np.array([plot_ras])),
            plot_decs,
            c=probs,
            vmin=0.0,
            vmax=max(self.data[self.key]),
            s=size,
        )

        plt.scatter(
            self.wrap_around_180(np.array([single_ras])),
            single_decs,
            c=single_probs,
            vmin=0.0,
            vmax=max(self.data[self.key]),
            s=size,
            cmap="gray",
        )

        plt.scatter(
            self.wrap_around_180(np.array([veto_ras])), veto_decs, color="red", s=size
        )

        red_patch = mpatches.Patch(color="red", label="Not observed")
        gray_patch = mpatches.Patch(color="gray", label="Observed once")
        plt.legend(handles=[red_patch, gray_patch])

        self.overlap_prob = 100.0 * np.sum(probs)
        once_observed_prob = 100 * (np.sum(probs) + np.sum(single_probs))
        message = (
            f"In total, {once_observed_prob} % of the contour was observed at least once. \n "
            f"In total, {self.overlap_prob} % of the contour was observed at least twice. \n"
            "THIS DOES NOT INCLUDE CHIP GAPS!!!"
        )

        print(message)

        self.area = (2.0 * base_ztf_rad) ** 2 * float(len(overlapping_fields))
        self.n_fields = len(overlapping_fields)
        self.overlap_fields = overlapping_fields

        print(
            f"{self.n_fields} fields were covered, covering approximately {self.area} sq deg."
        )
        return fig, message

    # def plot_overlap_with_observations(self, fields=None, pid=None, first_det_window_days=None, min_sep=0.01):
    #
    #     try:
    #         nside = self.ligo_nside
    #     except AttributeError:
    #         nside = self.nside
    #
    #     fig = plt.figure()
    #     plt.subplot(projection="aitoff")
    #
    #     double_in_plane_probs = []
    #     single_in_plane_prob = []
    #
    #     if fields is None:
    #         mns = self.get_multi_night_summary(first_det_window_days)
    #
    #     else:
    #
    #         class MNS:
    #             def __init__(self, data):
    #                 self.data = pandas.DataFrame(data, columns=["field", "ra", "dec", "datetime"])
    #
    #         data = []
    #
    #         for f in fields:
    #             ra, dec = ztfquery_fields.field_to_coords(f)[0]
    #             for i in range(2):
    #                 t = Time(self.t_min.jd + 0.1*i, format="jd").utc
    #                 t.format = "isot"
    #                 t = t.value
    #                 data.append([f, ra, dec, t])
    #
    #         mns = MNS(data)
    #
    #     data = mns.data.copy()
    #
    #     if pid is not None:
    #         pid_mask = data["pid"] == str(pid)
    #         data = data[pid_mask]
    #
    #     obs_times = np.array([Time(data["datetime"].iat[i], format="isot", scale="utc")
    #                           for i in range(len(data))])
    #
    #
    #     if first_det_window_days is not None:
    #         first_det_mask = [x < Time(self.t_min.jd + first_det_window_days, format="jd").utc for x in obs_times]
    #         data = data[first_det_mask]
    #         obs_times = obs_times[first_det_mask]
    #
    #
    #     pix_obs_times = dict()
    #
    #     print("Unpacking observations")
    #     pix_map = dict()
    #     for i, obs_time in enumerate(tqdm(obs_times)):
    #         radec = [f"{data['ra'].iat[i]} {data['dec'].iat[i]}"]
    #         coords = SkyCoord(radec, unit=(u.hourangle, u.deg))
    #         # pix = get_quadrant_ipix(nside, data["ra"].iat[i], data["dec"].iat[i])
    #         pix = get_quadrant_ipix(nside, coords.ra.deg[0], coords.dec.deg[0])
    #
    #         field = data["field"].iat[i]
    #
    #         flat_pix = []
    #
    #         for sub_list in pix:
    #             for p in sub_list:
    #                 flat_pix.append(p)
    #
    #         flat_pix = list(set(flat_pix))
    #
    #         t = obs_time.jd
    #         for p in flat_pix:
    #             if p not in pix_obs_times.keys():
    #                 pix_obs_times[p] = [t]
    #             else:
    #                 pix_obs_times[p] += [t]
    #
    #             if p not in pix_map.keys():
    #                 pix_map[p] = [field]
    #             else:
    #                 pix_map[p] += [field]
    #
    #
    #     npix = hp.nside2npix(nside)
    #     theta, phi = hp.pix2ang(nside, np.arange(npix), nest=False)
    #     radecs = SkyCoord(ra=phi * u.rad, dec=(0.5 * np.pi - theta) * u.rad)
    #     idx = np.where(np.abs(radecs.galactic.b.deg) <= 10.0)[0]
    #
    #     double_in_plane_pixels = []
    #     double_in_plane_probs = []
    #     single_in_plane_pixels = []
    #     single_in_plane_prob = []
    #     veto_pixels = []
    #     plane_pixels = []
    #     plane_probs = []
    #     times = []
    #     double_no_plane_prob = []
    #     double_no_plane_pixels = []
    #     single_no_plane_prob = []
    #     single_no_plane_pixels = []
    #
    #     overlapping_fields = []
    #     for i, p in enumerate(tqdm(hp.nest2ring(nside, self.pixel_nos))):
    #         if p in pix_obs_times.keys():
    #             if p in idx:
    #                 plane_pixels.append(p)
    #                 plane_probs.append(self.map_probs[i])
    #
    #             obs = pix_obs_times[p]
    #
    #             # check which healpix are observed twice
    #             if max(obs) - min(obs) > min_sep:
    #                 # is it in galactic plane or not?
    #                 if p not in idx:
    #                     double_no_plane_prob.append(self.map_probs[i])
    #                     double_no_plane_pixels.append(p)
    #                 else:
    #                     double_in_plane_probs.append(self.map_probs[i])
    #                     double_in_plane_pixels.append(p)
    #
    #             else:
    #                 if p not in idx:
    #                     single_no_plane_pixels.append(p)
    #                     single_no_plane_prob.append(self.map_probs[i])
    #                 else:
    #                     single_in_plane_prob.append(self.map_probs[i])
    #                     single_in_plane_pixels.append(p)
    #
    #             overlapping_fields += pix_map[p]
    #
    #             times += obs
    #         else:
    #             veto_pixels.append(p)
    #
    #     # print(f"double no plane prob = {double_no_plane_prob}")
    #     # print(f"probs = {probs}")
    #     # print(f"single no plane prob = {single_no_plane_prob}")
    #     # print(f"single probs = {single_probs}")
    #
    #     overlapping_fields = sorted(list(set(overlapping_fields)))
    #     self.overlap_fields = list(set(overlapping_fields))
    #
    #     self.overlap_prob = np.sum(double_in_plane_probs + double_no_plane_prob) * 100.
    #
    #     size = hp.max_pixrad(nside) ** 2 * 50.
    #
    #     veto_pos = np.array([hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in veto_pixels]).T
    #
    #     if len(veto_pos) > 0:
    #
    #         plt.scatter(self.wrap_around_180(np.radians(veto_pos[0])), np.radians(veto_pos[1]),
    #                     color="red", s=size)
    #
    #     plane_pos = np.array([hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in plane_pixels]).T
    #
    #     if len(plane_pos) > 0:
    #
    #         plt.scatter(self.wrap_around_180(np.radians(plane_pos[0])), np.radians(plane_pos[1]),
    #                     color="green", s=size)
    #
    #     single_pos = np.array([hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in single_no_plane_pixels]).T
    #
    #     if len(single_pos) > 0:
    #         plt.scatter(self.wrap_around_180(np.radians(single_pos[0])), np.radians(single_pos[1]),
    #                     c=single_no_plane_prob, vmin=0., vmax=max(self.data[self.key]), s=size, cmap='gray')
    #
    #     plot_pos = np.array([hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in double_no_plane_pixels]).T
    #
    #     if len(plot_pos) > 0:
    #         plt.scatter(self.wrap_around_180(np.radians(plot_pos[0])), np.radians(plot_pos[1]),
    #                     c=double_no_plane_prob, vmin=0., vmax=max(self.data[self.key]), s=size)
    #
    #     red_patch = mpatches.Patch(color='red', label='Not observed')
    #     gray_patch = mpatches.Patch(color='gray', label='Observed once')
    #     violet_patch = mpatches.Patch(color='green', label='Observed Galactic Plane (|b|<10)')
    #     plt.legend(handles=[red_patch, gray_patch, violet_patch])
    #
    #     message = "In total, {0:.2f} % of the contour was observed at least once. \n " \
    #               "This estimate includes {1:.2f} % of the contour " \
    #               "at a galactic latitude <10 deg. \n " \
    #               "In total, {2:.2f} % of the contour was observed at least twice. \n" \
    #               "In total, {3:.2f} % of the contour was observed at least twice, " \
    #               "and excluding low galactic latitudes. \n" \
    #               "These estimates accounts for chip gaps.".format(
    #         100 * (np.sum(double_in_plane_probs) + np.sum(single_in_plane_prob) + np.sum(single_no_plane_prob) + np.sum(double_no_plane_prob)),
    #         100 * np.sum(plane_probs),
    #         100.*(np.sum(double_in_plane_probs) + np.sum(double_no_plane_prob)),
    #         100.*np.sum(double_no_plane_prob)
    #         )
    #
    #     all_pix = single_in_plane_pixels + double_in_plane_pixels
    #
    #     n_pixels = len(single_in_plane_pixels + double_in_plane_pixels + double_no_plane_pixels + single_no_plane_pixels)
    #     n_double = len(double_no_plane_pixels + double_in_plane_pixels)
    #     n_plane = len(plane_pixels)
    #
    #     self.area = hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_pixels
    #     self.double_extragalactic_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_double
    #     plane_area  = hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_plane
    #     try:
    #         self.first_obs = Time(min(times), format="jd")
    #         self.first_obs.utc.format = "isot"
    #         self.last_obs = Time(max(times), format="jd")
    #         self.last_obs.utc.format = "isot"
    #
    #     except ValueError:
    #         raise Exception(f"No observations of this field were found at any time after {self.t_min.jd:.2f} JD ({self.t_min}). "
    #                         "Coverage overlap is 0%!")
    #
    #     print("Observations started at {0}".format(self.first_obs.jd))
    #
    #     self.overlap_fields = overlapping_fields
    #
    #     #     area = (2. * base_ztf_rad)**2 * float(len(overlapping_fields))
    #     #     n_fields = len(overlapping_fields)
    #
    #     print("{0} pixels were covered, covering approximately {1:.2g} sq deg.".format(
    #         n_pixels, self.area))
    #     print("{0} pixels were covered at least twice (b>10), covering approximately {1:.2g} sq deg.".format(
    #         n_double, self.double_extragalactic_area))
    #     print("{0} pixels were covered at low galactic latitude, covering approximately {1:.2g} sq deg.".format(
    #         n_plane, plane_area))
    #     return fig, message

    def plot_overlap_with_observations(
        self, fields=None, pid=None, first_det_window_days=None, min_sep=0.01
    ):

        try:
            nside = self.ligo_nside
        except AttributeError:
            nside = self.nside

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

        self.logger.info(
            f"No RA/Dec found by ztfquery for fields {veto_fields}. These observation have to be ignored."
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

        pix_obs_times = dict()

        self.logger.info(f"Most recent observation found is {obs_times[-1]}")

        self.logger.info("Unpacking observations")
        pix_map = dict()

        for i, obs_time in enumerate(tqdm(obs_times)):

            # radec = [f"{data['ra'].iat[i]} {data['dec'].iat[i]}"]
            #
            # coords = SkyCoord(radec, unit=(u.deg, u.deg))
            pix = get_quadrant_ipix(nside, data["ra"].iat[i], data["dec"].iat[i])

            field = data["field"].iat[i]

            flat_pix = []

            for sub_list in pix:
                for p in sub_list:
                    flat_pix.append(p)

            flat_pix = list(set(flat_pix))

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

        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix), nest=False)
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
        for i, p in enumerate(tqdm(hp.nest2ring(nside, self.pixel_nos))):

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

        # print(f"double no plane prob = {double_no_plane_prob}")
        # print(f"probs = {probs}")
        # print(f"single no plane prob = {single_no_plane_prob}")
        # print(f"single probs = {single_probs}")

        overlapping_fields = sorted(list(set(overlapping_fields)))
        self.overlap_fields = list(set(overlapping_fields))

        self.overlap_prob = np.sum(double_in_plane_probs + double_no_plane_prob) * 100.0

        size = hp.max_pixrad(nside) ** 2 * 50.0

        veto_pos = np.array(
            [hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in veto_pixels]
        ).T

        if len(veto_pos) > 0:

            plt.scatter(
                self.wrap_around_180(np.radians(veto_pos[0])),
                np.radians(veto_pos[1]),
                color="red",
                s=size,
            )

        plane_pos = np.array(
            [hp.pixelfunc.pix2ang(nside, i, lonlat=True) for i in plane_pixels]
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
                hp.pixelfunc.pix2ang(nside, i, lonlat=True)
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
                hp.pixelfunc.pix2ang(nside, i, lonlat=True)
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
            "These estimates accounts for chip gaps.".format(
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

        self.area = hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_pixels
        self.double_extragalactic_area = (
            hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_double
        )
        plane_area = hp.pixelfunc.nside2pixarea(nside, degrees=True) * n_plane

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

        print(
            f"{n_pixels} pixels were covered, covering approximately {self.area:.2g} sq deg."
        )
        print(
            f"{n_double} pixels were covered at least twice (b>10), covering approximately {self.double_extragalactic_area:.2g} sq deg."
        )
        print(
            f"{n_plane} pixels were covered at low galactic latitude, covering approximately {plane_area:.2g} sq deg."
        )
        return fig, message

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

    def export_fields(self):
        mask = np.array([x in self.overlap_fields for x in self.mns.data["field"]])
        lim_mag = [20.5 for _ in range(np.sum(mask))]
        coincident_obs = self.mns.data[mask].assign(lim_mag=lim_mag)
        print(
            coincident_obs[["field", "pid", "datetime", "lim_mag", "exp"]].to_csv(
                index=False, sep=" "
            )
        )
