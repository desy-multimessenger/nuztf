#!/usr/bin/env python
# coding: utf-8

import unittest
import logging
from astropy.time import Time
from astropy import units as u
from nuztf.observations import get_obs_summary, get_obs_summary_skyvision


class TestCoverage(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_lightcurve(self):
        self.logger.info("\n\n Testing observation log parsing \n\n")

        t_start = Time(2458865.96, format="jd")
        t_end = Time(2458866.96, format="jd")

        res = get_obs_summary(t_start, t_end)

        expected = {
            "obsid": 111223429.0,
            "field": 3.550000e02,
            "obsjd": 2458866.734294,
            "seeing": 3.4250149727,
            "limmag": 19.998298645,
            "exposure_time": 3.000000e01,
            "fid": 2.000000e00,
            "processed_fraction": 1.000000e00,
        }

        self.assertEqual(len(res.data), 211)

        for (name, val) in expected.items():
            self.assertEqual(res.data.iloc[0][name], val)

        # res2 = get_obs_summary_skyvision(
        #     t_start,
        #     t_end
        # )
        #
        # for name in ["obsjd"]:
        #     val = expected[name]
        #     self.assertAlmostEqual(res2.data.iloc[0][name], val)
        #
        # print("Lens", len(res2.data), len(res.data))

        # get_obs_summary(Time.now() - 1.0 * u.day, Time.now())
