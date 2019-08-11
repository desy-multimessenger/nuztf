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

base_ligo_dir = os.path.join(Path().absolute(), "LIGO_skymaps")
candidate_output_dir = os.path.join(Path().absolute(), "LIGO_candidates")

gw_run_config = {
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

class RetractionError(Exception):
   """Base class for retracted event"""
   pass


class GravWaveScanner(AmpelWizard):

    def __init__(self, gw_name=None, rev=None, logger=None, prob_threshold=0.9, cone_nside=64, t_max=None):
        self.gw_path, self.output_path = self.get_superevent(gw_name, rev)
        self.parsed_file = self.read_map()
        self.merger_time = Time(self.parsed_file[1].header["DATE-OBS"], format="isot", scale="utc")

        if t_max is not None:
            self.t_max = t_max
        else:
            self.t_max = Time.now()

        print("MERGER TIME: {0}".format(self.merger_time))

        self.data = self.parsed_file[1].data
        self.prob_map = hp.read_map(self.gw_path)
        self.prob_threshold = prob_threshold
        self.pixel_threshold = self.find_pixel_threshold(self.data["PROB"])
        self.map_coords = self.unpack_skymap()
        AmpelWizard.__init__(self, run_config=gw_run_config, logger=logger, cone_nside=cone_nside)


    def filter_f_no_prv(self, res):
        # Positive detection
        if res['candidate']['isdiffpos'] in ["t", "1"]:
            if self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
                return True

        return False

    def filter_f_history(self, res):
        # Veto past detections, but not past upper limits

        for prv_detection in res["prv_candidates"]:
            if np.logical_and(prv_detection["isdiffpos"] is not None, prv_detection["jd"] < self.merger_time.jd):
                return False

        # Require 2 detections

        n_detections = len([x for x in res["prv_candidates"] if np.logical_and(
            x["isdiffpos"] is not None, x["jd"] > self.merger_time.jd)])

        if n_detections < 1:
            return False

        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            return False

        return True

    def get_superevent(self, name, rev):
        if name is None:
            superevent_iterator = ligo_client.superevents('category: Production')
            superevent_ids = [superevent['superevent_id'] for superevent in superevent_iterator]
            name = superevent_ids[0]

        # latest_gw = ligo_client.superevent(name)

        voevents = ligo_client.voevents(name).json()["voevents"]

        if rev is None:
            rev = len(voevents)

        elif rev > len(voevents):
            raise ("Revision {0} not found".format(rev))

        latest_voevent = voevents[rev - 1]
        print("Found voevent {0}".format(latest_voevent["filename"]))

        if "Retraction" in latest_voevent["filename"]:
            raise RetractionError("The specified LIGO event, {0}, was retracted.".format(latest_voevent["filename"]))

        response = requests.get(latest_voevent["links"]["file"])

        root = lxml.etree.fromstring(response.content)
        params = {elem.attrib['name']:
                      elem.attrib['value']
                  for elem in root.iterfind('.//Param')}

        latest_skymap = params["skymap_fits"]

        print("Latest skymap URL: {0}".format(latest_skymap))

        base_file_name = os.path.basename(latest_skymap)
        savepath = os.path.join(base_ligo_dir, "{0}_{1}_{2}".format(
            name, latest_voevent["N"], base_file_name))

        print("Saving to: {0}".format(savepath))
        response = requests.get(latest_skymap)

        with open(savepath, "wb") as f:
            f.write(response.content)

        output_file = "{0}/{1}_{2}.pdf".format(candidate_output_dir, name, latest_voevent["N"])

        return savepath, output_file

    def read_map(self, ):
        print("Reading file: {0}".format(self.gw_path))
        f = fits.open(self.gw_path)
        return f

    def find_pixel_threshold(self, data):
        print("")
        ranked_pixels = np.sort(data)[::-1]
        int_sum = 0.0
        pixel_threshold = 0.0

        for i, prob in enumerate(ranked_pixels):
            int_sum += prob
            if int_sum > self.prob_threshold:
                print("Threshold found! \n To reach {0}% of probability, pixels with "
                      "probability greater than {1} are included".format(
                    int_sum * 100., prob))
                pixel_threshold = prob
                break

        return pixel_threshold

    @staticmethod
    def extract_ra_dec(nside, index):
        (colat, ra) = hp.pix2ang(nside, index, nest=True)
        dec = np.pi / 2. - colat
        return (ra, dec)

    @staticmethod
    def extract_npix(nside, ra, dec):
        colat = np.pi / 2. - dec
        return hp.ang2pix(nside, colat, ra, nest=True)

    def unpack_skymap(self):

        ligo_nside = hp.npix2nside(len(self.data["PROB"]))

        threshold = self.find_pixel_threshold(self.data["PROB"])

        mask = self.data["PROB"] > threshold

        map_coords = []

        print("Checking which pixels are within the contour:")

        for i in tqdm(range(hp.nside2npix(ligo_nside))):
            if mask[i]:
                map_coords.append(self.extract_ra_dec(ligo_nside, i))

        print("Total pixel area: {0} degrees".format(
            hp.nside2pixarea(ligo_nside, degrees=True) * float(len(map_coords))))

        map_coords = np.array(map_coords, dtype=np.dtype([("ra", np.float),
                                                          ("dec", np.float)]))

        return map_coords

    def find_cone_coords(self):
        cone_ids = []

        for ra, dec in self.map_coords:
            cone_ids.append(self.extract_npix(self.cone_nside, ra, dec))

        cone_ids = list(set(cone_ids))

        cone_coords = []

        for i in tqdm(cone_ids):
            cone_coords.append(self.extract_ra_dec(self.cone_nside, i))

        cone_coords = np.array(
            cone_coords, dtype=np.dtype([("ra", np.float), ("dec", np.float)])
        )

        return cone_ids, cone_coords

    @staticmethod
    def wrap_around_180(ra):
        ra[ra > np.pi] -= 2 * np.pi
        return ra

    def plot_skymap(self):

        plt.subplot(projection="aitoff")

        mask = self.data["PROB"] > self.pixel_threshold

        sc = plt.scatter(self.wrap_around_180(self.map_coords["ra"]), self.map_coords["dec"],
                         c=self.data["PROB"][mask], vmin=0., vmax=max(self.data["PROB"]), s=1e-4)
        plt.title("LIGO SKYMAP")
        plt.show()

        plt.subplot(projection="aitoff")

        sc = plt.scatter(self.wrap_around_180(self.cone_coords["ra"]), self.cone_coords["dec"])
        plt.title("CONE REGION")
        plt.show()

    def interpolate_map(self, ra_deg, dec_deg):
        colat = np.pi / 2. - np.radians(dec_deg)
        long = np.radians(ra_deg)
        return hp.pixelfunc.get_interp_val(self.prob_map, colat, long)

    def in_contour(self, ra_deg, dec_deg):
        return self.interpolate_map(ra_deg, dec_deg) > self.pixel_threshold


if __name__=="__main__":

    import logging
    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    gw = GravWaveScanner(logger=logger)
    gw.scan_cones()
