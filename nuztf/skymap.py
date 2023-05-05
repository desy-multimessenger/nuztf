#!/usr/bin/env python
# coding: utf-8

import json
import logging
import os
import time

import healpy as hp
import lxml.etree
import numpy as np
import requests
import wget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy_healpix import HEALPix
from ligo.gracedb.exceptions import HTTPError
from ligo.skymap.moc import rasterize
from lxml import html
from numpy.lib.recfunctions import append_fields
from ztfquery.io import LOCALSOURCE


class EventNotFound(Exception):
    """Base class for non-existing event"""

    pass


class Skymap:
    def __init__(
        self,
        event: str = None,
        rev: int = None,
        prob_threshold: float = 0.9,
        custom_prefix: str = "",
    ):
        self.base_skymap_dir = os.path.join(LOCALSOURCE, f"{custom_prefix}skymaps")
        self.candidate_output_dir = os.path.join(
            LOCALSOURCE, f"{custom_prefix}candidates"
        )
        self.candidate_cache = os.path.join(LOCALSOURCE, f"{custom_prefix}cache")

        for entry in [
            self.base_skymap_dir,
            self.candidate_output_dir,
            self.candidate_cache,
        ]:
            if not os.path.exists(entry):
                os.makedirs(entry)

        self.logger = logging.getLogger(__name__)

        self.prob_threshold = prob_threshold

        if ".fit" in event:
            basename = os.path.basename(event)

            self.skymap_path = os.path.join(self.base_skymap_dir, basename)

            if event[:8] == "https://":
                self.logger.info(f"Downloading from: {event}")
                self.skymap_path = os.path.join(
                    self.base_skymap_dir, os.path.basename(event[7:])
                )
                wget.download(event, self.skymap_path)

                self.summary_path = os.path.join(
                    self.candidate_output_dir,
                    f"{os.path.basename(event)}_{self.prob_threshold}",
                )

            else:
                self.summary_path = os.path.join(
                    self.candidate_output_dir,
                    f"{os.path.basename(event)}_{self.prob_threshold}",
                )

            self.event_name = os.path.basename(event[7:])

        # elif np.sum([x in event for x in ["grb", "GRB"]]) > 0:
        elif "grb" in event or "GRB" in event:
            self.get_grb_skymap(event_name=event)

        elif np.sum([x in event for x in ["s", "S", "gw", "GW"]]) > 0:
            self.skymap_path, self.summary_path, self.event_name = self.get_gw_skymap(
                event_name=event, rev=rev
            )
        else:
            raise Exception(
                f"Event {event} not recognised as a fits file, a GRB or a GW event."
            )

        (
            self.data,
            self.t_obs,
            self.hpm,
            self.key,
            self.dist,
            self.dist_unc,
        ) = self.read_map()

        t_min = Time(self.t_obs, format="isot", scale="utc")

        self.logger.info(f"Event time: {t_min}")
        self.logger.info("Reading map")

        self.pixel_threshold = self.find_pixel_threshold(self.data[self.key])

        # self.cache_dir = os.path.join(skymap_candidate_cache, self.event_name)

        # if not os.path.exists(self.cache_dir):
        #     os.makedirs(self.cache_dir)

    def get_grb_skymap(self, event_name: str):
        """ """
        if event_name is None:
            raise ValueError(
                "event_name must be provided for GRBs. They must have the form 'GRB210729A"
            )

        event_year_short = event_name[3:5]
        event_year = "20" + event_year_short
        event_month = event_name[5:7]
        event_day = event_name[7:9]
        event_letter = event_name[9]
        event_number = ord(event_letter) - 65

        # get possible skymap URLs

        url = f"https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/{event_year}"

        page_overview = requests.get(url)
        webpage_overview = html.fromstring(page_overview.content)

        links_overview = webpage_overview.xpath("//a/@href")

        links_for_date = []

        for link in links_overview:
            if link[2:8] == f"{event_year_short}{event_month}{event_day}":
                links_for_date.append(url + "/" + link + "current/")

        if len(links_for_date) > 1:
            self.logger.info(
                f"Found multiple events. Will choose the one corresponding the GRB letter {event_letter}"
            )

        event_url = links_for_date[event_number]

        page_event = requests.get(event_url)
        webpage_event = html.fromstring(page_event.content)
        links_event = webpage_event.xpath("//a/@href")

        for link in links_event:
            if link[0:11] == "glg_healpix":
                final_link = event_url + link
                break

        self.skymap_path = os.path.join(self.base_skymap_dir, link)

        if os.path.isfile(self.skymap_path):
            self.logger.info(
                f"Continuing with saved skymap. Located at {self.skymap_path}"
            )
        else:
            self.logger.info(f"Downloading skymap and saving to {self.skymap_path}")
            wget.download(final_link, self.skymap_path)

        self.summary_path = f"{self.base_skymap_dir}/{event_name}_{self.prob_threshold}"

        self.event_name = event_name

    def get_gw_skymap(self, event_name: str, rev: int):
        """ """
        from ligo.gracedb.rest import GraceDb

        ligo_client = GraceDb()

        self.logger.info("Obtaining skymap from GraceDB")

        if event_name is None:
            superevent_iterator = ligo_client.superevents("category: Production")
            superevent_ids = [
                superevent["superevent_id"] for superevent in superevent_iterator
            ]
            event_name = superevent_ids[0]

        try:
            res = ligo_client.voevents(event_name)
            if res.status_code == 200:
                voevents = res.json()["voevents"]
            else:
                raise EventNotFound(
                    f"The specified LIGO event, {event_name}, was not found on GraceDB. Please check that you entered the correct event name."
                )
        except HTTPError:
            raise EventNotFound(
                f"The specified LIGO event, {event_name}, was not found on GraceDB. Please check that you entered the correct event name."
            )

        if rev is None:
            rev = len(voevents)

        elif rev > len(voevents):
            raise Exception("Revision {0} not found".format(rev))

        self.rev = rev

        latest_voevent = voevents[rev - 1]
        self.logger.info(f"Found voevent {latest_voevent['filename']}")

        if "Retraction" in latest_voevent["filename"]:
            raise RetractionError(
                f"The specified LIGO event, {latest_voevent['filename']}, was retracted."
            )

        response = requests.get(latest_voevent["links"]["file"])

        root = lxml.etree.fromstring(response.content)
        params = {
            elem.attrib["name"]: elem.attrib["value"]
            for elem in root.iterfind(".//Param")
        }

        latest_skymap = params["skymap_fits"]

        self.logger.info(f"Latest skymap URL: {latest_skymap}")

        base_file_name = os.path.basename(latest_skymap)
        savepath = os.path.join(
            self.base_skymap_dir,
            f"{event_name}_{latest_voevent['N']}_{base_file_name}",
        )

        self.logger.info(f"Saving to: {savepath}")
        response = requests.get(latest_skymap)

        with open(savepath, "wb") as f:
            f.write(response.content)

        summary_path = f"{self.base_skymap_dir}/{event_name}_{latest_voevent['N']}_{self.prob_threshold}"

        return savepath, summary_path, event_name

    def read_map(
        self,
    ):
        """Read the skymap"""

        self.logger.info(f"Reading file: {self.skymap_path}")

        with fits.open(self.skymap_path) as hdul:
            data = None
            h = hdul[0].header

            dist = None
            dist_unc = None
            t_obs = None
            ordering = None

            for x in hdul:
                if data is None:
                    if x.data is not None:
                        data = np.array(x.data)

                if "DISTMEAN" in x.header:
                    dist = x.header["DISTMEAN"]

                if "DISTSTD" in x.header:
                    dist_unc = x.header["DISTSTD"]

                if "DATE-OBS" in x.header:
                    t_obs = x.header["DATE-OBS"]

                elif "EVENTMJD" in x.header:
                    t_obs_mjd = x.header["EVENTMJD"]
                    t_obs = Time(t_obs_mjd, format="mjd").isot

                if "ORDERING" in x.header:
                    ordering = x.header["ORDERING"]

        if "PROB" in data.dtype.names:
            key = "PROB"
        elif "PROBABILITY" in data.dtype.names:
            key = "PROB"
            prob = np.array(data["PROBABILITY"]).flatten()
            data = append_fields(data, "PROB", prob)
        elif "PROBDENSITY" in data.dtype.names:
            key = "PROB"
            prob = np.array(data["PROBDENSITY"])
            data = append_fields(data, "PROB", prob)
        elif "T" in data.dtype.names:  # weird IceCube format
            key = "PROB"
            prob = np.array(data["T"]).flatten()
            data = append_fields(data, "PROB", prob)
        else:
            raise Exception(
                f"No recognised probability key in map. This is probably a weird one, right? "
                f"Found the following keys: {data.dtype.names}"
            )

        if ordering == "NUNIQ":
            self.logger.info("Rasterising skymap to convert to nested format")
            data = data[list(["UNIQ", key])]
            probs = rasterize(data, order=7)
            data = np.array(probs, dtype=np.dtype([("PROB", float)]))

        if not isinstance(data["PROB"][0], float):
            self.logger.info("Flattening skymap")
            probs = np.array(data["PROB"]).flatten()
            data = np.array(probs, dtype=np.dtype([("PROB", float)]))

        if "NSIDE" not in h.keys():
            h["NSIDE"] = hp.npix2nside(len(data[key]))

        data["PROB"] /= np.sum(data["PROB"])

        self.logger.info(f"Summed probability is {100. * np.sum(data['PROB']):.1f}%")

        if ordering == "RING":
            data["PROB"] = hp.pixelfunc.reorder(data["PROB"], inp="RING", out="NESTED")

        if ordering is not None:
            h["ORDERING"] = "NESTED"
        else:
            raise Exception(
                f"Error parsing fits file, no ordering found. "
                f"Please enter the ordewring (NESTED/RING/NUNIQ)"
            )

        hpm = HEALPix(nside=h["NSIDE"], order=h["ORDERING"], frame="icrs")

        return data, t_obs, hpm, key, dist, dist_unc

    def find_pixel_threshold(self, data):
        """ """

        ranked_pixels = np.sort(data)[::-1]
        int_sum = 0.0
        pixel_threshold = 0.0

        for i, prob in enumerate(ranked_pixels):
            int_sum += prob
            if int_sum > self.prob_threshold:
                self.logger.info(
                    f"Threshold found! \n To reach {int_sum * 100.0:.2f}% of probability, pixels with probability greater than {prob} are included."
                )
                pixel_threshold = prob
                break

        return pixel_threshold

    def interpolate_map(self, ra_deg, dec_deg):
        """ """
        interpol_map = self.hpm.interpolate_bilinear_skycoord(
            SkyCoord(ra_deg * u.deg, dec_deg * u.deg), self.data[self.key]
        )
        return interpol_map

    def in_contour(self, ra_deg, dec_deg):
        """ """
        return self.interpolate_map(ra_deg, dec_deg) > self.pixel_threshold
