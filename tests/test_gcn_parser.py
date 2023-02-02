#!/usr/bin/env python
# coding: utf-8

import unittest, logging

from nuztf.parse_nu_gcn import get_latest_gcn, gcn_url, find_gcn_no


class TestNeutrinoScanner(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_latest(self):
        self.logger.info("\n\n Testing parsing of GCNs \n\n")
        no = get_latest_gcn()
        self.logger.info(f"Latest alert is {no}")
        url = gcn_url(gcn_number=no)
        self.logger.info(f"URL is {url}")

    def test_named(self):
        name = "IC200620A"

        num = int(find_gcn_no(base_nu_name=name))

        ref = 27997

        self.logger.info(f"GCN number for {name} was found to be {num}")
        self.logger.info(f"Reference value was {ref}")

        self.assertEqual(num, ref)

        fakename = "IC130921A"

        self.logger.info(f"Searching for fictional alert {fakename}")

        gcn_nr = find_gcn_no(base_nu_name=fakename)

        if gcn_nr:
            raise Exception(
                f"Somehow found a GCN ({gcn_nr}) matching "
                f"fictional neutrino alert {fakename}"
            )
        else:
            self.logger.info("No GCN found, as expected.")
