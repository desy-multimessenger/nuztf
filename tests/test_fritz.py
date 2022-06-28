#!/usr/bin/env python3
# License: BSD-3-Clause

import unittest
import logging
from astropy.time import Time
import astropy.units as u
from nuztf.fritz import save_source_to_group, delete_source_from_group


class TestFritz(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_fritz(self):

        self.logger.info("\n\n Testing Fritz queries \n\n")

        ztf_id = "ZTF21abyonuw"

        response = save_source_to_group(ztf_id, group_id=1430)
        self.assertEqual(response.status_code, 200)

        response = delete_source_from_group(ztf_id, group_id=1430)
        self.assertEqual(response.status_code, 200)
