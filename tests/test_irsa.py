#!/usr/bin/env python
# coding: utf-8

import unittest
import logging
import nuztf
from nuztf.irsa import plot_irsa_lightcurve, load_irsa


class TestIrsa(unittest.TestCase):
    def setUp(self):
        logging.getLogger("nuztf.irsa").setLevel(logging.DEBUG)
        logging.getLogger("nuztf.observations").setLevel(logging.DEBUG)

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_lightcurve(self):
        self.logger.info("\n\n Testing IRSA \n\n")
        self.logger.info("Getting lightcurve from IPAC.")

        res = load_irsa(77.358185, 5.693148, 0.5, TIME=[59204.1, 59210.9])

        res = res[res["programid"] == 1]
        expected = 2

        self.logger.info(f"Found {len(res)} entries. Expected: {expected}")

        self.assertEqual(len(res), expected)

        src_names = ["PKS1502+106", "SN2021gpw"]
        nu_names = ["IC190730A", "IC211216B"]

        for i, src_name in enumerate(src_names):
            plot_irsa_lightcurve(
                source_name=src_name,
                nu_name=nu_names[i],
                check_obs=True,
                check_obs_lookback_weeks=1,
                query_irsa_for_logs=False,
            )
