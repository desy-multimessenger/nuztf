# coding: utf-8
import json
import logging
from pathlib import Path
import requests
from hashlib import sha256

import healpy as hp
from astropy.time import Time
import numpy as np
from ligo.skymap.postprocess.util import find_greedy_credible_levels, smooth_ud_grade
from ligo.skymap.io.fits import read_sky_map
from gcn_kafka import Consumer

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.paths import SKYMAP_DIR, GCN_KAFKA_CACHE


logger = logging.getLogger(__name__)

TEST_GCN_KAFKA_TOPIC = "gcn.notices.icecube.test.gold_bronze_track_alerts"


class NeutrinoKafkaScanner(NeutrinoScanner):
    def __init__(
        self,
        alert: dict,
        prob_threshold: float = 0.9,
        map_path: Path = None,
    ):
        map_path = map_path or self.download_map(alert["healpix_url"])
        hpx_map, header = read_sky_map(str(map_path))

        # to be compatible with code relying on the 90% rectangle
        # we parse accordingly from the header
        # TODO:
        #  The position values are taken from the header for now
        #  because the position does not match for the example alert
        #  between alert json and header!
        ra = [header["RA"], header["RA_ERR_MINUS"], header["RA_ERR_PLUS"]]
        dec = [header["DEC"], header["DEC_ERR_MINUS"], header["DEC_ERR_PLUS"]]

        self.skymap = np.array(hpx_map, dtype=[("PROB", float)])
        self.skymap_header = header
        self.prob_threshold = prob_threshold
        self.alert = alert

        # set up an attribute to store credible levels
        # calculated in self.unpack_skymap
        self._credible_levels = None

        super().__init__(
            manual_args=(alert["event_name"][0], ra, dec, Time(alert["trigger_time"])),
            cone_nside=128,
            t_precursor=None,
            config=None,
            output_nside=hp.npix2nside(len(self.skymap)),
        )

    @property
    def credible_levels(self):
        # the private attribute is calculated in self.unpack_skymap
        return self._credible_levels

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
        else:
            skymap = self.skymap

        # the map contains the probability per pixel so we need to convert
        # to cumulative probability contained within a contour
        credible_levels = find_greedy_credible_levels(skymap["PROB"])
        self._credible_levels = credible_levels

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


def alert_filename(alert: dict) -> Path:
    nu_name = alert["event_name"][0]
    h = sha256(json.dumps(alert).encode()).hexdigest()
    return GCN_KAFKA_CACHE / f"{nu_name}_{h}.json"


def save_alert(alert: dict, overwrite: bool = False):
    fn = alert_filename(alert)
    if fn.exists() and not overwrite:
        new_fn = fn.parent / (fn.stem + "_new" + fn.suffix)
        logger.warning(f"{fn} already exists!")
        fn = new_fn
    with fn.open("w") as f:
        json.dump(alert, f)


def load_alert(nu_name: str) -> dict:
    files = list(GCN_KAFKA_CACHE.glob(f"{nu_name}*.json"))
    if len(files) == 0:
        raise FileNotFoundError(f"No file found in {GCN_KAFKA_CACHE} for {nu_name}!")
    if len(files) > 1:
        raise RuntimeError(
            f"More than one file found in {GCN_KAFKA_CACHE} for {nu_name}!"
        )
    with open(GCN_KAFKA_CACHE / files[0], "r") as f:
        return json.load(f)


def listen(
    client_id: str,
    client_secret: str,
    topics: list[str] = None,
    console=None,
    gcn_filename: str | None = None,
):
    consumer = Consumer(
        client_id=client_id,
        client_secret=client_secret,
    )

    # Subscribe to topics and receive alerts
    consumer.subscribe(topics or [TEST_GCN_KAFKA_TOPIC])
    logger.info("Connected ...")
    while True:
        for message in consumer.consume(timeout=1):
            if message.error():
                logger.error(message.error())
                continue
            # Print the topic and message ID
            logger.info("Found alert!")
            logger.info(f"topic={message.topic()}, offset={message.offset()}")
            value = message.value()
            logger.debug(value)

            # save alert for future reference
            save_alert(message)

            # instantiate scanner
            nu = NeutrinoKafkaScanner(alert=message)
            nu.scan(console=console, gcn_filename=gcn_filename)


def scan_saved(nu_name: str, console=None, gcn_filename: str | None = None):
    alert = load_alert(nu_name)
    nu = NeutrinoKafkaScanner(alert)
    nu.scan(console=console, gcn_filename=gcn_filename)
