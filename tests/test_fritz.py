#!/usr/bin/env python3
# License: BSD-3-Clause

import logging
import unittest

import astropy.units as u
import backoff
from astropy.time import Time

from nuztf.fritz import delete_source_from_group, save_source_to_group

ztf_id = "ZTF21abyonuw"
ztf_group = 1430


class TestFritz(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_fritz(self):
        self.logger.info("Deactivated for now, Fritz too unstable!")

    #     self.logger.info("\n\n Testing Fritz queries \n\n")

    #     self.logger.debug(f"Deleting {ztf_id} from group {ztf_group}")
    #     delete_source_from_group(ztf_id, group_id=ztf_group)

    #     self.logger.debug(f"Saving {ztf_id} to group {ztf_group}")
    #     response = save_source_to_group(ztf_id, group_id=ztf_group)
    #     self.assertEqual(response.status_code, 200)

    #     self.logger.debug(f"Deleting {ztf_id} from group {ztf_group}")
    #     response = delete_source_from_group(ztf_id, group_id=ztf_group)
    #     self.assertEqual(response.status_code, 200)

    # def tearDown(self):
    #     self.logger.debug(f"Deleting {ztf_id} from group {ztf_group}")
    #     delete_source_from_group(ztf_id, group_id=ztf_group)
