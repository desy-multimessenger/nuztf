# coding: utf-8

from pathlib import Path
import requests

import healpy as hp
from astropy.time import Time
import numpy as np
from ligo.skymap.postprocess.util import find_greedy_credible_levels, smooth_ud_grade
from ligo.skymap.io.fits import read_sky_map

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.paths import SKYMAP_DIR


class NeutrinoKafkaScanner(NeutrinoScanner):
    def __init__(
        self,
        alert: dict,
        prob_threshold: float = 0.9,
    ):
        map_path = self.download_map(alert["healpix_url"])
        hpx_map, header = read_sky_map(str(map_path))

        # to be compatible with code relying on the 90% rectangle
        # we parse accordingly from the header
        ra = [alert["RA"], header["RA_ERR_MINUS"], header["RA_ERR_PLUS"]]
        dec = [alert["DEC"], header["DEC_ERR_MINUS"], header["DEC_ERR_PLUS"]]

        self.skymap = np.array(hpx_map, dtype=[("PROB", float)])
        self.skymap_header = header

        # the map contains the probability per pixel so we need to convert
        # to cumulative probability contained within a contour
        self.credible_levels = find_greedy_credible_levels(self.skymap)

        self.prob_threshold = prob_threshold

        self.alert = alert

        super().__init__(
            manual_args=(alert["event_name"][0], ra, dec, Time(alert["trigger_time"])),
            cone_nside=128,
            t_precursor=None,
            config=None,
            output_nside=hp.npix2nside(len(self.skymap)),
        )

    @staticmethod
    def download_map(url: str) -> Path:
        response = requests.get(url)
        response.raise_for_status()
        local_path = SKYMAP_DIR / Path(url).name
        with open(local_path, "wb") as f:
            f.write(response.content)
        return local_path

    def unpack_skymap(self, output_nside: None | int = None):
        map_nside = hp.npix2nside(len(self.skymap))

        # interpolate skymap to output nside if necessary
        if output_nside is None:
            output_nside = map_nside
        if output_nside != map_nside:
            self.logger.info(
                f"Interpolating input skymap from nside {map_nside} to {output_nside}"
            )
            skymap = smooth_ud_grade(self.skymap, output_nside)
            credible_levels = find_greedy_credible_levels(self.skymap)
        else:
            skymap = self.skymap
            credible_levels = self.credible_levels

        # find healpix indices inside credible region
        healpix_indices = np.where(credible_levels <= self.prob_threshold)[0]
        map_coords = hp.pix2ang(map_nside, healpix_indices, lonlat=True)

        # calculate pixel area
        total_pixel_area = hp.nside2pixarea(output_nside, degrees=True) * float(
            len(healpix_indices)
        )
        return (
            map_coords,
            healpix_indices,
            map_nside,
            skymap[healpix_indices],
            skymap,
            total_pixel_area,
            "PROB",
        )

    def in_contour(self, ra_deg, dec_deg):
        # check whether the position is inside the credible region defined by the threshold
        return (
            hp.get_interp_val(self.credible_levels, ra_deg, dec_deg, lonlat=True)
            <= self.prob_threshold
        )
