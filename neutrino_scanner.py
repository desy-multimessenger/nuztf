from ampel_magic import AmpelWizard
from astropy.time import Time
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from tqdm import tqdm
from ligo.gracedb.rest import GraceDb
import os
from pathlib import Path
import requests
import matplotlib.patches as mpatches
import lxml.etree
from astropy.io import fits

nu_candidate_output_dir = os.path.join(Path().absolute(), "Neutrino_candidates")

nu_run_config = {
    "MIN_NDET": 1,  # Default:2
    "MIN_TSPAN": -1,  # Default 0, but that rejects everything!
    "MAX_TSPAN": 365,
    "MIN_RB": 0.3,
    "MAX_FWHM": 5.5,
    "MAX_ELONG": 1.4,
    "MAX_MAGDIFF": 1.0,
    "MAX_NBAD": 2,
    "MIN_DIST_TO_SSO": 20,
    "MIN_GAL_LAT": -1.0,  # Default: 14
    "GAIA_RS": 20,
    "GAIA_PM_SIGNIF": 3,
    "GAIA_PLX_SIGNIF": 3,
    "GAIA_VETO_GMAG_MIN": 9,
    "GAIA_VETO_GMAG_MAX": 20,
    "GAIA_EXCESSNOISE_SIG_MAX": 999,
    "PS1_SGVETO_RAD": 1,
    "PS1_SGVETO_SGTH": 0.8,
    "PS1_CONFUSION_RAD": 3,
    "PS1_CONFUSION_SG_TOL": 0.1
}

class ParsingError(Exception):
   """Base class for parsing error"""
   pass


class NeutrinoScanner(AmpelWizard):

    def __init__(self, nu_name=None, manual_args=None, gcn_no=None, logger=None, cone_nside=128):

        if manual_args is None:

            if nu_name is not None:
                gcn_no = self.find_gcn_no(nu_name)

            elif gcn_no is None:
                gcn_no = self.get_latest_gcn()

            nu_name, author, ra, dec, nu_time = self.parse_gcn(gcn_no)

        else:
            (nu_name, ra, dec, nu_time) = manual_args
            author = None

        self.nu_name = nu_name
        self.author = author
        self.gcn_no = gcn_no

        print("Neutrino time: {0}".format(nu_time))

        self.ra_max = max(ra[1:]) + ra[0]
        self.ra_min = min(ra[1:]) + ra[0]
        self.dec_max = max(dec[1:]) + dec[0]
        self.dec_min = min(dec[1:]) + dec[0]

        print("Coordinates: RA = {0} ({1} - {2})".format(ra[0], self.ra_min, self.ra_max))
        print("Coordinates: Dec = {0} ({1} - {2})".format(dec[0], self.dec_min, self.dec_max))

        self.output_path = "{0}/{1}.pdf".format(nu_candidate_output_dir, nu_name)
        AmpelWizard.__init__(self, t_min=nu_time, run_config=nu_run_config, logger=logger, cone_nside=cone_nside)
        self.default_t_max = Time.now()
        self.prob_threshold = 0.9
        self.area = (self.ra_max - self.ra_min) * (self.dec_max - self.dec_min) * np.cos(np.radians(dec[0]))
        print("Projected Area: {0}".format(self.area))

    @staticmethod
    def gcn_url(gcn_number):
        return "https://gcn.gsfc.nasa.gov/gcn3/{0}.gcn3".format(gcn_number)

    def get_name(self):
        return self.nu_name

    def get_full_name(self):
        return "neutrino event {0} ({1} et. al, GCN {2})".format(self.get_name(), self.author, self.gcn_no)

    @staticmethod
    def get_tiling_line():
        return ""

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
        for x in line.replace(",", " ").split(" "):
            try:
                vals.append(float("".join([y for y in x if y not in ["(", "+", ")", "[", "]"]])))
            except ValueError:
                pass
        return vals

    def parse_gcn(self, gcn_number):
        url = self.gcn_url(gcn_number)
        page = requests.get(url)
        print("Found GCN: {0}".format(url))
        name = author = ra = dec = time = None
        for line in page.text.splitlines():
            if "SUBJECT" in line:
                name = line.split(" - ")[0].split(": ")[1]
            elif "FROM" in line:
                base = line.split("at")[0].split(": ")[1].split(" ")
                author = [x for x in base if x != ""][1]
            elif np.logical_and(np.sum([x in line for x in ["Ra", "RA"]]) > 0, ra is None):
                ra = self.strip_numbers(line)
            elif np.logical_and(np.sum([x in line for x in ["Dec", "DEC"]]) > 0, dec is None):
                dec = self.strip_numbers(line)
            elif np.logical_and(np.sum([x in line for x in ["Time", "TIME"]]) > 0, dec is None):
                raw_time = [x for x in  line.split(" ") if x not in ["Time", "", "UT", "UTC"]][1]
                raw_date = name.split("-")[1][:6]
                ut_time = "20{0}-{1}-{2}T{3}".format(raw_date[0:2], raw_date[2:4], raw_date[4:6], raw_time)
                time = Time(ut_time, format='isot', scale='utc')

        try:

            if np.sum([x is not None for x in [name, author, ra, dec, time]]) == 5:
                if np.logical_and(len(ra) == 3, len(dec) == 3):

                    return name, author, ra, dec, time

        except:
            pass

        raise ParsingError("Error parsing GCN {0}".format(url))


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

    def find_gcn_no(self, base_nu_name):
        page = requests.get("https://gcn.gsfc.nasa.gov/gcn3_archive.html")

        gcn_no = None
        name = None

        nu_name = str(base_nu_name)

        while not nu_name[0].isdigit():
            nu_name = nu_name[1:]

        for line in page.text.splitlines():
            if np.logical_and("IceCube observation of a high-energy neutrino" in line, nu_name in line):
                res = line.split(">")
                if gcn_no is None:
                    gcn_no = "".join([x for x in res[2] if x.isdigit()])
                    name = res[3].split(" - ")[0]
                    print("Found match to {0}: {1}".format(base_nu_name, name))
                else:
                    raise Exception("Multiple matches found to {0}".format(base_nu_name))

        if name is None:
            raise ParsingError("No GCN match found for {0}".format(base_nu_name))

        print("Match is {0} (GCN #{1})".format(name, gcn_no))

        return gcn_no

    def get_latest_gcn(self):
        latest = self.parse_gcn_archive()[0]
        print("Latest GCN is {0} (GCN #{1})".format(latest[0], latest[1]))
        return latest[1]


    def filter_f_no_prv(self, res):
        # Positive detection
        if res['candidate']['isdiffpos'] in ["t", "1"]:
            # if self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            return True

        # print("{0} has no positive detections".format(res["objectId"]))

        return False

    def filter_f_history(self, res):

        # print("Checking {0}".format(res["objectId"]))

        # Require 2 detections

        n_detections = len([x for x in res["prv_candidates"] if x["isdiffpos"] is not None])

        if n_detections < 1:
            # print("{0} has insufficient detection".format(res["objectId"]))
            return False

        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
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
                ra_deg > self.ra_min - scan_radius,
                ra_deg < self.ra_max + scan_radius
            ):
                if np.logical_and(
                    dec_deg > self.dec_min - scan_radius,
                    dec_deg < self.dec_max + scan_radius
                ):
                    cone_coords.append((ra, dec))
                    cone_ids.append(i)

        cone_coords = np.array(
            cone_coords, dtype=np.dtype([("ra", np.float), ("dec", np.float)])
        )

        return cone_ids, cone_coords

    def in_contour(self, ra_deg, dec_deg):

        in_ra = np.logical_and(
            ra_deg > self.ra_min,
            ra_deg < self.ra_max
        )
        in_dec = np.logical_and(
            dec_deg > self.dec_min,
            dec_deg < self.dec_max
        )

        return np.logical_and(in_ra, in_dec)

    def plot_overlap_with_observations(self):
        fig = plt.figure()
        plt.subplot(projection="aitoff")

        probs = []
        single_probs = []

        mns = self.get_multi_night_summary()

        ras = np.degrees(self.wrap_around_180(np.array([
            np.radians(float(x)) for x in mns.data["ra"]])))

        fields = list(mns.data["field"])

        plot_ras = []
        plot_decs = []

        single_ras = []
        single_decs = []

        veto_ras = []
        veto_decs = []

        overlapping_fields = []

        base_ztf_rad = 3.5
        ztf_dec_deg = 30.

        prob_cone = 1./float(len(self.cone_coords))

        for j, (ra, dec) in enumerate(tqdm(self.cone_coords)):
            ra_deg = np.degrees(self.wrap_around_180(np.array([ra])))
            # ra_deg = self.wrap_around_180(np.array(np.degrees(ra)))
            dec_deg = np.degrees(dec)
            ztf_rad = base_ztf_rad / (np.cos(dec - np.radians(ztf_dec_deg))*np.cos(dec))

            n_obs = 0

            for i, x in enumerate(self.get_multi_night_summary().data["dec"]):
                if np.logical_and(not dec_deg < float(x) - ztf_rad, not dec_deg > float(x) + ztf_rad):
                    if abs(dec_deg - ztf_dec_deg) < 70.:
                        if np.logical_and(not ra_deg < float(ras[i]) - ztf_rad, not ra_deg > float(ras[i]) + ztf_rad):
                            n_obs += 1
                            fid = fields[i]
                            if fid not in overlapping_fields:
                                overlapping_fields.append(fields[i])

            if n_obs > 1:
                probs.append(prob_cone)
                plot_ras.append(ra)
                plot_decs.append(dec)

            elif n_obs > 0:
                single_probs.append(prob_cone)
                single_ras.append(ra)
                single_decs.append(dec)

            else:
                veto_ras.append(ra)
                veto_decs.append(dec)

        overlapping_fields = list(set(overlapping_fields))

        obs_times = np.array([Time(mns.data["UT_START"].iat[i], format="isot", scale="utc")
                     for i in range(len(mns.data)) if mns.data["field"].iat[i] in overlapping_fields])

        self.first_obs = min(obs_times)
        self.last_obs = max(obs_times)

        size = hp.max_pixrad(self.cone_nside, degrees=True)**2

        # print(hp.max_pixrad(self.ligo_nside, degrees=True)**2 * np.pi, size)

        plt.scatter(self.wrap_around_180(np.array([plot_ras])), plot_decs,
                    c=probs, vmin=0., vmax=prob_cone, s=size)

        plt.scatter(self.wrap_around_180(np.array([single_ras])), single_decs,
                    c=single_probs, vmin=0., vmax=prob_cone, s=size, cmap='gray')

        plt.scatter(self.wrap_around_180(np.array([veto_ras])), veto_decs, color="red", s=size)

        red_patch = mpatches.Patch(color='red', label='Not observed')
        gray_patch = mpatches.Patch(color='gray', label='Observed once')
        plt.legend(handles=[red_patch, gray_patch])

        self.overlap_prob = 90.*np.sum(probs)

        message = "In total, {0} % of the contour was observed at least once. \n " \
                  "In total, {1} % of the contour was observed at least twice. \n" \
                  "THIS DOES NOT INCLUDE CHIP GAPS!!!".format(
            90. * (np.sum(probs) + np.sum(single_probs)), self.overlap_prob)

        print(message)

        self.area = (2. * base_ztf_rad)**2 * float(len(overlapping_fields))
        self.n_fields = len(overlapping_fields)

        print("{0} fields were covered, covering approximately {1} sq deg.".format(
            self.n_fields, self.area))
        return fig, message


if __name__=="__main__":

    import logging
    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    nu = NeutrinoScanner(logger=logger)
    nu.scan_cones()
