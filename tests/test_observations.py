#!/usr/bin/env python
# coding: utf-8

import logging
import unittest

from astropy import units as u
from astropy.time import Time

from nuztf.observations import get_obs_summary_irsa, get_obs_summary_skyvision


class TestCoverage(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_lightcurve(self):
        self.logger.info("\n\n Testing observation log parsing \n\n")

        t_start = Time(2458865.96, format="jd")
        t_end = Time(2458866.96, format="jd")

        # res = get_obs_summary_irsa(t_start, t_end)

        # expected = {
        #     "obsid": 111223429.0,
        #     "field": 3.550000e02,
        #     "obsjd": 2458866.734294,
        #     "seeing": 3.4250149727,
        #     "limmag": 19.998298645,
        #     "exposure_time": 3.000000e01,
        #     "fid": 2.000000e00,
        #     "processed_fraction": 1.000000e00,
        # }

        # self.assertEqual(len(res.data), 211)

        # for name, val in expected.items():
        #     self.assertEqual(res.data.iloc[0][name], val)

        res2 = get_obs_summary_skyvision(t_start, t_end)

        expected_2 = {
            "datetime": "2020-01-18T05:37:54.356",
            "date": "2020-01-18",
            "exptime": 30.0,
            "totalexptime": 207.428,
            "fid": 2,
            "field": 355,
            "pid": 1,
            "ra": "+05:02:46.24",
            "dec": "-09:51:00",
            "totaltime": 207.428,
            "base_name": "ztf_20200118233438_000355_zr",
            "obsjd": 2458866.7346568983,
        }

        for name, val in expected_2.items():
            self.assertEqual(res2.data.iloc[0][name], val)
