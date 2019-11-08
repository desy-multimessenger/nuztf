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
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord
import matplotlib.patches as mpatches
import fitsio
from astropy import units as u
import wget
from numpy.lib.recfunctions import append_fields
import pandas
from ztfquery import fields as ztfquery_fields
from gwemopt.ztf_tiling import get_quadrant_ipix
import logging

# Setup LIGO client

ligo_client = GraceDb()

try:
    r = ligo_client.ping()
except HTTPError as e:
    raise(e.message)

base_ligo_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "LIGO_skymaps")
ligo_candidate_output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "LIGO_candidates")

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
    "MIN_GAL_LAT": 0.,  # Default: 14
    "GAIA_RS": 10.,
    "GAIA_PM_SIGNIF": 3,
    "GAIA_PLX_SIGNIF": 3,
    "GAIA_VETO_GMAG_MIN": 9,
    "GAIA_VETO_GMAG_MAX": 20,
    "GAIA_EXCESSNOISE_SIG_MAX": 999,
    "PS1_SGVETO_RAD": 1.,
    "PS1_SGVETO_SGTH": 0.8,
    "PS1_CONFUSION_RAD": 3.,
    "PS1_CONFUSION_SG_TOL": 0.1
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
                print("Downloading from: {0}".format(gw_file))
                self.gw_path = "{0}/{1}".format(base_ligo_dir, os.path.basename(gw_file[7:]))
                wget.download(gw_file, self.gw_path)

            self.output_path = "{0}/{1}_{2}.pdf".format(
                ligo_candidate_output_dir, os.path.basename(gw_file), self.prob_threshold)
            self.gw_name = os.path.basename(gw_file[7:])
        self.data, t_obs, self.hpm, self.key = self.read_map()

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

    def get_tiling_line(self):

        too_bool = len(["ToO" in x for x in self.get_multi_night_summary().data.qcomment])> 0

        if too_bool:
            return "The tiling was optimally determined and triggered using the GROWTH Target of Opportunity marshal (Coughlin et al. 2019a, Kasliwal et al. 2019b). "
        else:
            return ""

    def get_obs_line(self):
        return "Each exposure was 30s with a typical depth of 20.5 mag."

    def get_overlap_line(self):
        return "We covered {0:.1f}% of the enclosed probability " \
               "based on the bayestar map in {1} sq deg. " \
               "This estimate accounts for chip gaps. ".format(
            self.overlap_prob, self.area)

    @staticmethod
    def remove_variability_line():
        return ", and removing candidates with history of " \
               "variability prior to the merger time"

    def candidate_text(self, name, first_detection, lul_lim, lul_jd):
        text = "{0}, first detected {1:.1f} hours after merger, " \
        "was not detected {2:.1f} days prior to a depth of {3:.2f}. ".format(
            name, 24. * (first_detection - self.t_min.jd), first_detection - lul_jd, lul_lim
        )
        return text

    def filter_f_no_prv(self, res):

        # Positive detection
        if res['candidate']['isdiffpos'] not in ["t", "1"]:
            logging.debug("Negative subtraction")
            return False

        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            logging.debug("Transient is too old")
            return False

        # Check contour
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            logging.debug("Outside of contour")
            return False

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
            logging.debug("Transient is too old")
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > self.default_t_max.jd:
            logging.debug("Transient is too new")
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

        output_file = "{0}/{1}_{2}_{3}.pdf".format(ligo_candidate_output_dir, name, latest_voevent["N"],
                                                   self.prob_threshold)

        return savepath, output_file, name

    def read_map(self, ):
        print("Reading file: {0}".format(self.gw_path))
        data, h = fitsio.read(self.gw_path, header=True)#columns=["PROB"],
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

        hpm = HEALPix(nside=h["NSIDE"], order=h["ORDERING"], frame='icrs')
        # with fits.open(self.gw_path) as hdul:
        #     print("Opened file")
        #     t_obs = hdul[0].header
        #     print(t_obs)
        #     print("read merger time")
        #     data = hdul[1].data
        #     print("Read data")
        return data, t_obs, hpm, key

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

    def plot_overlap_with_observations(self, fields=None, pid=None, first_det_window_days=None):
        fig = plt.figure()
        plt.subplot(projection="aitoff")

        probs = []
        single_probs = []

        if fields is None:
            mns = self.get_multi_night_summary(first_det_window_days)

        else:

            class MNS:
                def __init__(self, data):
                    self.data = pandas.DataFrame(data, columns=["field", "ra", "dec", "UT_START"])

            data = []

            for f in fields:
                ra, dec = ztfquery_fields.field_to_coords(f)[0]
                t = Time(Time.now().jd + 1., format="jd").utc
                t.format = "isot"
                t = t.value
                for _ in range(2):
                    data.append([f, ra, dec, t])

            mns = MNS(data)

        data = mns.data.copy()

        if pid is not None:
            pid_mask = data["pid"] == str(pid)
            data = data[pid_mask]

        obs_times = np.array([Time(data["UT_START"].iat[i], format="isot", scale="utc")
                              for i in range(len(data))])

        if first_det_window_days is not None:
            first_det_mask = [x < Time(self.t_min.jd + first_det_window_days, format="jd").utc for x in obs_times]
            data = data[first_det_mask]
            obs_times = obs_times[first_det_mask]

        pix_obs_times = dict()

        print("Unpacking observations")

        for i, obs_time in enumerate(tqdm(obs_times)):
            pix = get_quadrant_ipix(self.ligo_nside, data["ra"].iat[i], data["dec"].iat[i])

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

        plot_pixels = []
        probs = []
        single_pixels = []
        single_probs = []
        veto_pixels = []
        times = []

        print("Checking skymap")

        # hp.nest2ring(gw.ligo_nside,

        for i, p in enumerate(tqdm(hp.nest2ring(self.ligo_nside, self.pixel_nos))):

            if p in pix_obs_times.keys():

                obs = pix_obs_times[p]

                if max(obs) - min(obs) > 0.01:
                    plot_pixels.append(p)
                    probs.append(self.map_probs[i])
                else:
                    single_pixels.append(p)
                    single_probs.append(self.map_probs[i])

                times += obs
            else:
                veto_pixels.append(p)

        self.overlap_prob = np.sum(probs) * 100.

        size = hp.max_pixrad(self.ligo_nside) ** 2 * 50.

        print(size, self.ligo_nside)

        veto_pos = np.array([hp.pixelfunc.pix2ang(self.ligo_nside, i, lonlat=True) for i in veto_pixels]).T

        plt.scatter(self.wrap_around_180(np.radians(veto_pos[0])), np.radians(veto_pos[1]),
                    color="red", s=size)

        single_pos = np.array([hp.pixelfunc.pix2ang(self.ligo_nside, i, lonlat=True) for i in single_pixels]).T

        plt.scatter(self.wrap_around_180(np.radians(single_pos[0])), np.radians(single_pos[1]),
                    c=single_probs, vmin=0., vmax=max(self.data[self.key]), s=size, cmap='gray')

        plot_pos = np.array([hp.pixelfunc.pix2ang(self.ligo_nside, i, lonlat=True) for i in plot_pixels]).T

        if len(plot_pos) > 0:
            plt.scatter(self.wrap_around_180(np.radians(plot_pos[0])), np.radians(plot_pos[1]),
                        c=probs, vmin=0., vmax=max(self.data[self.key]), s=size)

        red_patch = mpatches.Patch(color='red', label='Not observed')
        gray_patch = mpatches.Patch(color='gray', label='Observed once')
        plt.legend(handles=[red_patch, gray_patch])

        message = "In total, {0} % of the LIGO contour was observed at least once. \n " \
                  "In total, {1} % of the LIGO contour was observed at least twice. \n" \
                  "This estimate accounts for chip gaps.".format(
            100 * (np.sum(probs) + np.sum(single_probs)), self.overlap_prob)

        all_pix = single_pixels + plot_pixels

        n_pixels = len(single_pixels + plot_pixels)

        self.area = hp.pixelfunc.nside2pixarea(self.ligo_nside, degrees=True) * n_pixels

        self.first_obs = Time(min(times), format="jd")
        self.last_obs = Time(max(times), format="jd")

        print("Observations started at {0}".format(self.first_obs.jd))

        #     area = (2. * base_ztf_rad)**2 * float(len(overlapping_fields))
        #     n_fields = len(overlapping_fields)

        print("{0} pixels were covered, covering approximately {1} sq deg.".format(
            n_pixels, self.area))
        return fig, message

    def simple_plot_overlap_with_observations(self, fields=None, first_det_window_days=None):
        fig = plt.figure()
        plt.subplot(projection="aitoff")

        probs = []
        single_probs = []

        if fields is None:
            mns = self.get_multi_night_summary(first_det_window_days)

        else:

            class MNS:
                def __init__(self, data):
                    self.data = pandas.DataFrame(data, columns=["field", "ra", "dec", "UT_START"])

            data = []

            for f in fields:
                ra, dec = ztfquery_fields.field_to_coords(f)[0]
                t = Time(Time.now().jd + 1., format="jd").utc
                t.format = "isot"
                t = t.value
                for _ in range(2):
                    data.append([f, ra, dec, t])

            mns = MNS(data)

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

        for j, (ra, dec) in enumerate(tqdm(self.map_coords)):
            ra_deg = np.degrees(self.wrap_around_180(np.array([ra])))
            # ra_deg = self.wrap_around_180(np.array(np.degrees(ra)))
            dec_deg = np.degrees(dec)
            # (np.cos(dec - np.radians(ztf_dec_deg))
            ztf_rad = base_ztf_rad
            ztf_height = 3.7

            n_obs = 0

            for i, x in enumerate(mns.data["dec"]):
                if np.logical_and(not dec_deg < float(x) - ztf_height, not dec_deg > float(x) + ztf_height):
                    if abs(dec_deg - ztf_dec_deg) < 70.:
                        if np.logical_and(not ra_deg < float(ras[i]) - ztf_rad/ abs(np.cos(dec)),
                                          not ra_deg > float(ras[i]) + ztf_rad/ abs(np.cos(dec))):
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

        obs_times = np.array([Time(mns.data["UT_START"].iat[i], format="isot", scale="utc")
                     for i in range(len(mns.data)) if mns.data["field"].iat[i] in overlapping_fields])

        self.first_obs = min(obs_times)
        self.last_obs = max(obs_times)

        size = hp.max_pixrad(self.ligo_nside, degrees=True)**2

        # print(hp.max_pixrad(self.ligo_nside, degrees=True)**2 * np.pi, size)

        plt.scatter(self.wrap_around_180(np.array([plot_ras])), plot_decs,
                    c=probs, vmin=0., vmax=max(self.data[self.key]), s=size)

        plt.scatter(self.wrap_around_180(np.array([single_ras])), single_decs,
                    c=single_probs, vmin=0., vmax=max(self.data[self.key]), s=size, cmap='gray')

        plt.scatter(self.wrap_around_180(np.array([veto_ras])), veto_decs, color="red", s=size)

        red_patch = mpatches.Patch(color='red', label='Not observed')
        gray_patch = mpatches.Patch(color='gray', label='Observed once')
        plt.legend(handles=[red_patch, gray_patch])

        self.overlap_prob = 100.*np.sum(probs)

        message = "In total, {0} % of the LIGO contour was observed at least once. \n " \
                  "In total, {1} % of the LIGO contour was observed at least twice. \n" \
                  "THIS DOES NOT INCLUDE CHIP GAPS!!!".format(
            100 * (np.sum(probs) + np.sum(single_probs)), self.overlap_prob)

        print(message)

        self.area = (2. * base_ztf_rad)**2 * float(len(overlapping_fields))
        self.n_fields = len(overlapping_fields)
        self.overlap_fields = overlapping_fields

        print("{0} fields were covered, covering approximately {1} sq deg.".format(
            self.n_fields, self.area))
        return fig, message

    def interpolate_map(self, ra_deg, dec_deg):
        return self.hpm.interpolate_bilinear_skycoord(SkyCoord(ra_deg * u.deg, dec_deg * u.deg), self.data[self.key])

    def in_contour(self, ra_deg, dec_deg):
        return self.interpolate_map(ra_deg, dec_deg) > self.pixel_threshold

if __name__=="__main__":

    import logging
    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    gw = GravWaveScanner(logger=logger)
    gw.scan_cones()
