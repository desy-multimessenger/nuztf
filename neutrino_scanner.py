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
import lxml.etree
from astropy.io import fits

# Setup LIGO client

ligo_client = GraceDb()

try:
    r = ligo_client.ping()
except HTTPError as e:
    raise(e.message)

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


class NeutrinoScanner(AmpelWizard):

    def __init__(self, manual_args=None, gcn_no=None, logger=None, cone_nside=64):

        if manual_args is None:

            if gcn_no is None:
                gcn_no = self.get_latest_gcn()

            nu_name, ra, dec, nu_time = self.parse_gcn(gcn_no)

        else:
            (nu_name, ra, dec, nu_time) = manual_args

        self.ra_max = max(ra[1:]) + ra[0]
        self.ra_min = min(ra[1:]) + ra[0]
        self.dec_max = max(dec[1:]) + dec[0]
        self.dec_min = min(dec[1:]) + dec[0]

        self.output_path = "{0}/{1}.pdf".format(nu_candidate_output_dir, nu_name)
        self.t_max = Time.now()
        AmpelWizard.__init__(self, t_min=nu_time, run_config=nu_run_config, logger=logger, cone_nside=cone_nside)

    @staticmethod
    def gcn_url(gcn_number):
        return "https://gcn.gsfc.nasa.gov/gcn3/{0}.gcn3".format(gcn_number)

    @staticmethod
    def strip_numbers(line):
        vals = []
        for x in line.replace(",", " ").split(" "):
            try:
                vals.append(float("".join([y for y in x if y not in ["(", "+", ")", "[", "]"]])))
            except ValueError:
                pass

        return vals


    def parse_gcn(self, gcn_number, silent=True):
        url = self.gcn_url(gcn_number)
        page = requests.get(url)
        print("Found GCN: {0}".format(url))
        name = ra = dec = time = None
        for line in page.text.splitlines():
            if "SUBJECT" in line:
                name = line.split(" - ")[0].split(": ")[1]
            elif np.logical_and(np.sum([x in line for x in ["Ra", "RA"]]) > 0, ra is None):
                ra = self.strip_numbers(line)
            elif np.logical_and(np.sum([x in line for x in ["Dec", "DEC"]]) > 0, dec is None):
                dec = self.strip_numbers(line)
            elif np.logical_and(np.sum([x in line for x in ["Time", "TIME"]]) > 0, dec is None):
                raw_time = [x for x in  line.split(" ") if x not in ["Time", "", "UT", "UTC"]][1]
                raw_date = name.split("-")[1][:6]
                ut_time = "20{0}-{1}-{2}T{3}".format(raw_date[0:2], raw_date[2:4], raw_date[4:6], raw_time)
                time = Time(ut_time, format='isot', scale='utc')

        return name, ra, dec, time


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

    def get_latest_gcn(self):
        latest = self.parse_gcn_archive()[0]
        print("Latest GCN is {0} (GCN #{1})".format(latest[0], latest[1]))
        return latest[1]


    def filter_f_no_prv(self, res):
        # Positive detection
        if res['candidate']['isdiffpos'] in ["t", "1"]:
            if self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
                return True

        return False

    def filter_f_history(self, res):

        # Require 2 detections

        n_detections = len([x for x in res["prv_candidates"] if x["isdiffpos"] is not None])

        if n_detections < 1:
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


if __name__=="__main__":

    import logging
    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    nu = NeutrinoScanner(logger=logger)
    nu.scan_cones()
