#!/usr/bin/env python3
# License: BSD-3-Clause

import unittest
import logging
from astropy.time import Time
import astropy.units as u
from nuztf.ampel_api import (
    ampel_api_name,
    ampel_api_cutout,
    ampel_api_cone,
    ampel_api_timerange,
    ampel_api_healpix,
)


class TestAPI(unittest.TestCase):

    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_api(self):

        self.logger.info("\n\n Testing API queries \n\n")

        ztf_id = "ZTF21abyonuw"

        t_min_jd = Time("2019-04-01T00:00:00.123456789", format="isot", scale="utc").jd
        t_max_jd_timerange = t_min_jd + 0.127163

        logging.info(f"Retrieving alerts for {ztf_id}")
        api_name = ampel_api_name(
            ztf_name=ztf_id, with_history=True, with_cutouts=False, logger=self.logger
        )
        self.assertEqual(len(api_name), 1)
        self.logger.info(f"Successfully retrieved the alert for {ztf_id}")

        candid = api_name[0]["candid"]

        self.logger.info(f"Retrieving cutouts for {ztf_id}")
        api_cutouts = ampel_api_cutout(candid=candid)
        nr_cutouts = len(api_cutouts)

        ref = 5

        self.logger.info(f"Retrieved {nr_cutouts}. Reference value is {ref}")

        self.assertEqual(nr_cutouts, ref)

        self.logger.info("Commencing API cone search")

        t_max_jd_cone = Time("2021-10-07", format="isot").jd

        api_cone = ampel_api_cone(
            ra=30, dec=30, radius=0.1, t_max_jd=t_max_jd_cone, logger=self.logger
        )

        ztf_ids = []
        for entry in api_cone:
            ztf_ids.append(entry["objectId"])
        ztf_ids = list(set(ztf_ids))

        nr_transients = len(ztf_ids)
        ref = 94

        self.logger.info(f"Found {nr_transients} transients. Reference value is {ref}")

        self.assertEqual(nr_transients, ref)

        self.logger.info("Commencing API time search")
        api_time = ampel_api_timerange(
            t_min_jd=t_min_jd,
            t_max_jd=t_max_jd_timerange,
            chunk_size=2000,
            logger=self.logger,
        )

        ztf_ids = []
        for entry in api_time:
            ztf_ids.append(entry["objectId"])
        ztf_ids = list(set(ztf_ids))

        nr_transients = len(ztf_ids)
        ref = 1887

        self.logger.info(f"Found {nr_transients} transients. Reference value is {ref}")

        self.assertEqual(nr_transients, ref)

        self.logger.info("Commencing API healpix search")

        t_min = Time("2019-07-01T00:00:00.0", format="isot", scale="utc")
        t_max = t_min + 5 * u.day

        api_healpix = ampel_api_healpix(
            ipix=1000,
            nside=64,
            t_min_jd=t_min.jd,
            t_max_jd=t_max.jd,
            with_history=True,
            with_cutouts=True,
            logger=self.logger,
        )

        nr_transients = len(api_healpix)
        ref = 5

        self.assertEqual(nr_transients, ref)

        self.logger.info(f"Found {nr_transients} transients. Reference value is {ref}")
