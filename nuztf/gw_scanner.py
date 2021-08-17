from nuztf.ampel_magic import AmpelWizard
from astropy.time import Time
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
from tqdm import tqdm
from ligo.gracedb.rest import GraceDb
import os
import requests
import lxml.etree
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord
import fitsio
from astropy import units as u
import wget
from pathlib import Path
from numpy.lib.recfunctions import append_fields

# Setup LIGO client

ligo_client = GraceDb()

try:
    r = ligo_client.ping()
except HTTPError as e:
    raise(e.message)

base_ligo_dir = os.path.join(Path(__file__).resolve().parents[1], "../LIGO_skymaps")
ligo_candidate_output_dir = os.path.join(Path(__file__).resolve().parents[1], "../LIGO_candidates")

gw_run_config = {
    "min_ndet": 1,  # Default:2
    "min_tspan": -1,  # Default 0, but that rejects everything!
    "max_tspan": 365,
    "min_rb": 0.3,
    "max_fwhm": 5.5,
    "max_elong": 1.4,
    "max_magdiff": 1.0,
    "max_nbad": 2,
    "min_sso_dist": 20,
    "min_gal_lat": 0.,  # Default: 14
    "gaia_rs": 10.,
    "gaia_pm_signif": 3,
    "gaia_plx_signif": 3,
    "gaia_veto_gmag_min": 9,
    "gaia_veto_gmag_max": 20,
    "gaia_excessnoise_sig_max": 999,
    "ps1_sgveto_rad": 1.,
    "ps1_sgveto_th": 0.8,
    "ps1_confusion_rad": 3.,
    "ps1_confusion_sg_tol": 0.1
}

class RetractionError(Exception):
   """Base class for retracted event"""
   pass


class GravWaveScanner(AmpelWizard):

    def __init__(self, gw_name=None, gw_file=None, rev=None, logger=None, prob_threshold=0.95, cone_nside=64,
                 fast_query=False, n_days=None):

        self.prob_threshold = prob_threshold

        if gw_file is None:
            self.gw_path, self.output_path, self.gw_name = self.get_superevent(gw_name, rev)

        else:
            basename = os.path.basename(gw_file)
            self.gw_path = "{0}/{1}".format(base_ligo_dir, basename)
            if gw_file[:8] == "https://":
                logger.info("Downloading from: {0}".format(gw_file))
                self.gw_path = "{0}/{1}".format(base_ligo_dir, os.path.basename(gw_file[7:]))
                wget.download(gw_file, self.gw_path)

            self.output_path = "{0}/{1}_{2}.pdf".format(
                ligo_candidate_output_dir, os.path.basename(gw_file), self.prob_threshold)
            self.gw_name = os.path.basename(gw_file[7:])
        self.data, t_obs, self.hpm, self.key, self.dist, self.dist_unc = self.read_map()

        t_min = Time(t_obs, format="isot", scale="utc")

        logging.info("MERGER TIME: {0}".format(t_min))
        logging.info("Reading map")

        self.pixel_threshold = self.find_pixel_threshold(self.data[self.key])
        self.map_coords, self.pixel_nos, self.map_probs, self.ligo_nside, self.pixel_area = self.unpack_skymap()
        AmpelWizard.__init__(self, run_config=gw_run_config, t_min=t_min, logger=logger, cone_nside=cone_nside,
                             fast_query=fast_query)

        # By default, accept things detected within 72 hours of merger
        if n_days is None:
            n_days = 3.

        self.default_t_max = Time(self.t_min.jd + n_days, format="jd")

    def get_name(self):
        return self.gw_name

    def get_full_name(self):
        return self.gw_name

    # def get_tiling_line(self):
    #
    #     too_bool = len(["ToO" in x for x in self.get_multi_night_summary().data.qcomment])> 0
    #
    #     if too_bool:
    #         return "The tiling was optimally determined and triggered using the GROWTH Target of Opportunity marshal (Coughlin et al. 2019a, Kasliwal et al. 2019b). "
    #     else:
    #         return ""

    def get_obs_line(self):
        return "Each exposure was 30s with a typical depth of 20.5 mag."

    def get_overlap_line(self):
        return "We covered {0:.1f}% of the enclosed probability " \
               "based on the map in {1:.1f} sq deg. " \
               "This estimate accounts for chip gaps. ".format(
            self.overlap_prob, self.area)

    @staticmethod
    def remove_variability_line():
        return ", and removing candidates with history of " \
               "variability prior to the merger time"

    def candidate_text(self, name, first_detection, lul_lim, lul_jd):

        try:
            text = "{0}, first detected {1:.1f} hours after merger, " \
            "was not detected {2:.1f} days prior to a depth of {3:.2f}. ".format(
                name, 24. * (first_detection - self.t_min.jd), first_detection - lul_jd, lul_lim
            )
        except TypeError:
            text = "{0} had upper limit problems. PLEASE FILL IN NUMBERS BY HAND!!! ".format(name)
        return text

    def filter_f_no_prv(self, res):

        # Positive detection
        if res['candidate']['isdiffpos'] not in ["t", "1"]:
            logging.debug("Negative subtraction.")
            return False

        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            logging.debug("Transient is too old. (jdstarthist history predates event)")
            return False

        # Check contour
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            logging.debug("Outside of event contour.")
            return False

        logging.debug("Passed filter f (no prv)")

        return True

    # def fast_filter_f_no_prv(self, res):
    #
    #     # Positive detection
    #     if res['candidate']['isdiffpos'] not in ["t", "1"]:
    #         return False
    #
    #     # Veto old transients
    #     if res["candidate"]["jdstarthist"] < self.t_min.jd:
    #         return False
    #
    #     # Require 2 detections separated by 15 mins
    #     if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
    #         return False
    #
    #     return True

    def filter_f_history(self, res):
        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            logging.debug("Transient is too old. (jdstarthist history predates event)")
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > self.default_t_max.jd:
            logging.debug("Transient is too new. (jdstarthist too late after event)")
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            logging.debug("Not passed mover cut")
            return False

        # Require 2 positive detections
        old_detections = [x for x in res["prv_candidates"] if np.logical_and(
            x["isdiffpos"] is not None,
            x["jd"] > self.t_min.jd
        )]

        pos_detections = [x for x in old_detections if x['isdiffpos'] in ["t", "1"]]

        if len(pos_detections) < 1:
            logging.debug("Does not have two detections")
            return False

        return True

    def get_superevent(self, name, rev):
        if name is None:
            superevent_iterator = ligo_client.superevents('category: Production')
            superevent_ids = [superevent['superevent_id'] for superevent in superevent_iterator]
            name = superevent_ids[0]

        voevents = ligo_client.voevents(name).json()["voevents"]

        if rev is None:
            rev = len(voevents)

        elif rev > len(voevents):
            raise Exception("Revision {0} not found".format(rev))

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

        output_file = "{0}/{1}_{2}_{3}.pdf".format(ligo_candidate_output_dir, name, latest_voevent["N"],
                                                   self.prob_threshold)

        return savepath, output_file, name

    def read_map(self, ):
        print("Reading file: {0}".format(self.gw_path))
        data, h = fitsio.read(self.gw_path, header=True)#columns=["PROB"],
        if "DISTMEAN" not in h:
            dist = None
        else:
            dist = h["DISTMEAN"]
        if "DISTSTD" not in h:
            dist_unc = None
        else:
            dist_unc = h["DISTSTD"]
        if "DATE-OBS" not in h:
            t_obs = fitsio.read_header(self.gw_path)["DATE-OBS"]
        else:
            t_obs = h["DATE-OBS"]

        if "PROB" in data.dtype.names:
            key = "PROB"
        elif 'PROBABILITY' in data.dtype.names:
            key = 'PROB'
            prob = np.array(data["PROBABILITY"]).flatten()
            data = append_fields(data, "PROB", prob)
        else:
            raise Exception("No recognised probability key in map. This is probably a weird one, right?")

        if not isinstance(data[0], float):
            probs = np.array(data["PROB"]).flatten()
            # # print(data)
            # # print(drop_fields(data, "PROB"))
            # print(data.dtype.names)
            # data = drop_fields(data, "PROB", asrecarray=True)
            # print(type(data))
            data = np.array(probs, dtype=np.dtype([("PROB", np.float)]))

        logging.info(f"Summed probability is {100. * np.sum(data['PROB']):.1f}%")

        if h["ORDERING"] == "RING":
            data["PROB"] = hp.pixelfunc.reorder(data["PROB"], inp="RING", out="NESTED")
            h["ORDERING"] = "NESTED"

        hpm = HEALPix(nside=h["NSIDE"], order=h["ORDERING"], frame='icrs')

        # with fits.open(self.gw_path) as hdul:
        #     print("Opened file")
        #     t_obs = hdul[0].header
        #     print(t_obs)
        #     print("read merger time")
        #     data = hdul[1].data
        #     print("Read data")
        return data, t_obs, hpm, key, dist, dist_unc

    def find_pixel_threshold(self, data):

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

    def unpack_skymap(self):

        ligo_nside = hp.npix2nside(len(self.data[self.key]))

        threshold = self.find_pixel_threshold(self.data[self.key])

        mask = self.data[self.key] > threshold

        map_coords = []

        pixel_nos = []

        print("Checking which pixels are within the contour:")

        for i in tqdm(range(hp.nside2npix(ligo_nside))):
            if mask[i]:
                map_coords.append(self.extract_ra_dec(ligo_nside, i))
                pixel_nos.append(i)

        pixel_area = hp.nside2pixarea(ligo_nside, degrees=True) * float(len(map_coords))

        print("Total pixel area: {0} degrees".format(pixel_area))

        map_coords = np.array(map_coords, dtype=np.dtype([("ra", np.float),
                                                          ("dec", np.float)]))

        return map_coords, pixel_nos, self.data[self.key][mask], ligo_nside, pixel_area

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

    def plot_skymap(self):
        fig = plt.figure()
        plt.subplot(211, projection="aitoff")

        mask = self.data[self.key] > self.pixel_threshold

        size = hp.max_pixrad(self.ligo_nside, degrees=True) ** 2

        plt.scatter(self.wrap_around_180(self.map_coords["ra"]), self.map_coords["dec"],
                         c=self.data[self.key][mask], vmin=0., vmax=max(self.data[self.key]), s=size)
        plt.title("LIGO SKYMAP")

        plt.subplot(212, projection="aitoff")

        plt.scatter(self.wrap_around_180(self.cone_coords["ra"]), self.cone_coords["dec"])
        plt.title("CONE REGION")
        return fig

    def interpolate_map(self, ra_deg, dec_deg):
        return self.hpm.interpolate_bilinear_skycoord(SkyCoord(ra_deg * u.deg, dec_deg * u.deg), self.data[self.key])

    def in_contour(self, ra_deg, dec_deg):
        return self.interpolate_map(ra_deg, dec_deg) > self.pixel_threshold

if __name__=="__main__":

    import logging
    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.INFO)

    gw = GravWaveScanner(logger=logger)
    gw.scan_cones()
