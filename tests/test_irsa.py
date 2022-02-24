#!/usr/bin/env python
# coding: utf-8

import unittest
import logging
from nuztf.irsa import plot_irsa_lightcurve, load_irsa


class TestIrsa(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_lightcurve(self):
        self.logger.info("\n\n Testing IRSA \n\n")

        res = load_irsa(77.358185, 5.693148, 0.5, TIME=[58204.1, 59678.9])

        self.assertEqual(len(res), 436)

        src_names = ["PKS1502+106", "SN2021gpw"]
        nu_names = ["IC190730A", "IC211216B"]

        for i, src_name in enumerate(src_names):

            plot_irsa_lightcurve(
                source_name=src_name,
                nu_name=nu_names[i],
                check_obs=True,
                check_obs_lookback_weeks=1,
            )
