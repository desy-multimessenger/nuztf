#!/usr/bin/env python
# coding: utf-8

import unittest
import logging
from nuztf.irsa import plot_irsa_lightcurve


class TestIrsa(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_lightcurve(self):
        self.logger.info("\n\n Testing IRSA \n\n")

        src_name = "PKS1502+106"

        plot_irsa_lightcurve(
            source_name=src_name,
            nu_name="IC190730A",
            check_obs=True,
            check_obs_lookback_weeks=12,
        )
