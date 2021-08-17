#!/usr/bin/env python3
# coding: utf-8

from nuztf.ampel_magic import AmpelWizard
from astropy.time import Time
import healpy as hp
import numpy as np
from tqdm import tqdm
import os
from pathlib import Path
import requests
import logging
from ztf_plan_obs import gcn_parser

nu_candidate_output_dir = os.path.join(Path(__file__).resolve().parents[1], "Neutrino_candidates")

nu_run_config = {
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


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


class NeutrinoScanner(AmpelWizard):
    def __init__(
        self,
        nu_name=None,
        manual_args=None,
        gcn_no=None,
        logger=None,
        cone_nside=128,
        t_precursor=None,
    ):

        self.prob_threshold = 0.9

        if manual_args is None:

            if nu_name is not None:
                gcn_no = self.find_gcn_no(nu_name)

            elif gcn_no is None:
                gcn_no = self.get_latest_gcn()

            # nu_name, author, ra, dec, nu_time = self.parse_gcn(gcn_no)

            gcn_info = gcn_parser.parse_gcn_circular(gcn_no)
            print(gcn_info)

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

        print(f"Neutrino time: {nu_time}")

        self.ra_max = float(max(ra[1:]) + ra[0])
        self.ra_min = float(min(ra[1:]) + ra[0])
        self.dec_max = float(max(dec[1:]) + dec[0])
        self.dec_min = float(min(dec[1:]) + dec[0])

        print(f"Coordinates: RA = {ra[0]} ({self.ra_min} - {self.ra_max})")
        print(f"Coordinates: DEC = {dec[0]} ({self.dec_min} - {self.dec_max})")

        self.output_path = f"{nu_candidate_output_dir}/{nu_name}.pdf"
        AmpelWizard.__init__(
            self,
            t_min=nu_time,
            run_config=nu_run_config,
            logger=logger,
            cone_nside=cone_nside,
        )
        self.prob_threshold = 0.9
        self.area = (
            (self.ra_max - self.ra_min)
            * (self.dec_max - self.dec_min)
            * abs(np.cos(np.radians(dec[0])))
        )
        print(f"Projected Area: {self.area}")
        (
            self.map_coords,
            self.pixel_nos,
            self.nside,
            self.map_probs,
            self.data,
            self.key,
        ) = self.unpack_map()

    @staticmethod
    def gcn_url(gcn_number):
        return f"https://gcn.gsfc.nasa.gov/gcn3/{gcn_number}.gcn3"

    def get_name(self):
        return self.nu_name

    def get_full_name(self):
        return f"neutrino event {self.get_name()} ({self.author} et. al, GCN {self.gcn_no})"

    def get_overlap_line(self):
        return (
            f"We covered {self.area:.1f} sq deg, corresponding to {self.overlap_prob:.1f}% of the reported localization region. "
            "This estimate accounts for chip gaps. "
        )

    def candidate_text(self, name, first_detection, lul_lim, lul_jd):
        fd = Time(first_detection, format="mjd")

        text = f"{name} was first detected on {fd.utc}. "

        return text

    # @staticmethod
    # def get_tiling_line():
    #     return ""

    @staticmethod
    def get_obs_line():
        return "Each exposure was 300s with a typical depth of 21.0 mag."

    @staticmethod
    def remove_variability_line():
        return ""

    @staticmethod
    def strip_numbers(line):
        vals = []
        line = line.replace("- ", "-")
        line = line.replace("-", " -")

        for x in line.replace(",", " ").replace("/", " ").split(" "):
            try:
                vals.append(
                    float(
                        "".join(
                            [y for y in x if y not in ["(", "+", ")", "[", "]", "/"]]
                        )
                    )
                )
            except ValueError:
                pass
        return vals

    # def parse_gcn(self, gcn_number):
    #     url = self.gcn_url(gcn_number)
    #     page = requests.get(url)
    #     print(f"Found GCN: {url}")
    #     name = author = ra = dec = time = None
    #     print(page.text)
    #     for line in page.text.splitlines():
    #         line = "".join([x for x in line if x not in ["Â"]])
    #         if "SUBJECT" in line:
    #             name = line.split(" - ")[0].split(": ")[1]
    #         elif "FROM" in line:
    #             base = line.split("at")[0].split(": ")[1].split(" ")
    #             author = [x for x in base if x != ""][1]
    #         elif np.logical_and(np.sum([x in line for x in ["Ra", "RA"]]) > 0, ra is None):
    #             ra = self.strip_numbers(line)
    #         elif np.logical_and(np.sum([x in line for x in ["Dec", "DEC"]]) > 0, dec is None):
    #             dec = self.strip_numbers(line)
    #         elif np.logical_and(np.sum([x in line for x in ["Time", "TIME"]]) > 0, dec is None):
    #             raw_time = [x for x in  line.split(" ") if x not in ["Time", "", "UT", "UTC"]][1]
    #             raw_time = "".join([x for x in raw_time if np.logical_or(
    #                 x.isdigit(), x in [":", "."]
    #             )])
    #             raw_date = name.split("-")[1][:6]
    #             ut_time = f"20{raw_date[0:2]}-{raw_date[2:4]}-{raw_date[4:6]}T{raw_time}"
    #             time = Time(ut_time, format='isot', scale='utc')

    #     try:

    #         if np.sum([x is not None for x in [name, author, ra, dec, time]]) == 5:
    #             ra2 = ra[1]
    #             ra2 = ra2*-1
    #             ra.append(ra2)
    #             if np.logical_and(len(ra) == 3, len(dec) == 3):
    #                 print(name, author, ra, dec, time)
    #                 return name, author, ra, dec, time

    #     except:
    #         pass

    #     raise ParsingError(f"Error parsing GCN {0}".format(url))

    def parse_gcn_archive(self):
        page = requests.get("https://gcn.gsfc.nasa.gov/gcn3_archive.html")

        nu_circulars = []

        for line in page.text.splitlines():
            if "IceCube observation of a high-energy neutrino" in line:
                res = line.split(">")
                gcn_no = "".join([x for x in res[2] if x.isdigit()])
                name = res[3].split(" - ")[0]
                nu_circulars.append((name, gcn_no))

        return nu_circulars

    @staticmethod
    def parse_gcn_for_no(
        base_nu_name, url="https://gcn.gsfc.nasa.gov/gcn3_archive.html"
    ):

        print(f"Checking for GCN on {url}")

        nu_name = str(base_nu_name)

        page = requests.get(url)

        gcn_no = None
        name = None

        while not nu_name[0].isdigit():
            nu_name = nu_name[1:]

        latest_archive_no = None

        for line in page.text.splitlines():
            if np.logical_and(
                "IceCube observation of a high-energy neutrino" in line, nu_name in line
            ):
                res = line.split(">")
                if gcn_no is None:
                    gcn_no = "".join([x for x in res[2] if x.isdigit()])
                    name = res[3].split(" - ")[0]
                    print(f"Found match to {base_nu_name}: {name}")
                else:
                    raise Exception(f"Multiple matches found to {base_nu_name}")

            elif np.logical_and("gcn3_arch_old" in line, latest_archive_no is None):
                url = line.split('"')[1]
                latest_archive_no = int(url[13:].split(".")[0])

        return gcn_no, name, latest_archive_no

    def find_gcn_no(self, base_nu_name):

        gcn_no, name, latest_archive_no = self.parse_gcn_for_no(base_nu_name)

        if gcn_no is None:

            print(
                f"No GCN found for {base_nu_name} on GCN page, checking archive instead. "
                f"The latest page is {latest_archive_no}"
            )

            while np.logical_and(latest_archive_no > 0, gcn_no is None):
                gcn_no, name, _ = self.parse_gcn_for_no(
                    base_nu_name,
                    url=f"https://gcn.gsfc.nasa.gov/gcn3_arch_old{latest_archive_no}.html",
                )
                latest_archive_no -= 1

        # while

        if name is None:
            raise ParsingError("No GCN match found for {0}".format(base_nu_name))

        print(f"Match is {name} (GCN #{gcn_no})")

        return gcn_no

    def get_latest_gcn(self):
        latest = self.parse_gcn_archive()[0]
        print(f"Latest GCN is {latest[0]} (GCN #{latest[1]})")
        return latest[1]

    def filter_f_no_prv(self, res):

        # Positive detection
        if res["candidate"]["isdiffpos"] not in ["t", "1"]:
            logging.debug("Negative subtraction")
            return False

        try:
            if res["candidate"]["drb"] < 0.3:
                logging.debug("DRB too low")
                return False
        except (KeyError, TypeError) as e:
            pass

        # Check contour
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            logging.debug("Not in contour")
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            logging.debug("Does not have 2 detections separated  by >15 mins")
            return False

        return True

    def filter_f_history(self, res):
        # Require 2 detections

        n_detections = len(
            [x for x in res["prv_candidates"] if "isdiffpos" in x.keys()]
        )

        if n_detections < 1:
            logging.debug("{0} has insufficient detection".format(res["objectId"]))
            return False

        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            logging.debug("{0} not in contour".format(res["objectId"]))
            return False

        return True

    def find_cone_coords(self):

        cone_coords = []
        cone_ids = []

        scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))

        print("Finding search pixels:")

        for i in tqdm(range(hp.nside2npix(self.cone_nside))):
            ra, dec = self.extract_ra_dec(self.cone_nside, i)
            ra_deg = np.degrees(ra)
            dec_deg = np.degrees(dec)
            if np.logical_and(
                ra_deg > self.ra_min - scan_radius, ra_deg < self.ra_max + scan_radius
            ):
                if np.logical_and(
                    dec_deg > self.dec_min - scan_radius,
                    dec_deg < self.dec_max + scan_radius,
                ):
                    cone_coords.append((ra, dec))
                    cone_ids.append(i)

        cone_coords = np.array(
            cone_coords, dtype=np.dtype([("ra", np.float), ("dec", np.float)])
        )

        return cone_ids, cone_coords

    def in_contour(self, ra_deg, dec_deg):

        in_ra = np.logical_and(ra_deg > self.ra_min, ra_deg < self.ra_max)
        in_dec = np.logical_and(dec_deg > self.dec_min, dec_deg < self.dec_max)

        return np.logical_and(in_ra, in_dec)

    def unpack_map(self):

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
            if self.in_contour(np.degrees(ra), np.degrees(dec)):
                map_coords.append((ra, dec))
                pixel_nos.append(i)

        map_probs = np.ones_like(pixel_nos, dtype=np.float)
        map_probs /= np.sum(map_probs)

        key = "PROB"

        data = np.zeros(hp.nside2npix(nside), dtype=np.dtype([(key, np.float)]))
        data[np.array(pixel_nos)] = map_probs

        return map_coords, pixel_nos, nside, map_probs, data, key


if __name__ == "__main__":

    import logging

    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    nu = NeutrinoScanner(logger=logger)
    nu.scan_cones()
