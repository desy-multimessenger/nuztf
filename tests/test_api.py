#!/usr/bin/env python3
# License: BSD-3-Clause

import unittest
from astropy.time import Time
from nuztf.ampel_api import ampel_api_name, ampel_api_cutout, ampel_api_cone, ampel_api_timerange

from ampel.log.AmpelLogger import AmpelLogger

logger = AmpelLogger()

class TestAPI(unittest.TestCase):

    maxDiff = None
    
    ztf_id = "ZTF21abyonuw"

    t_min_jd = Time(
        "2019-04-01T00:00:00.123456789",
        format="isot",
        scale="utc"
    ).jd
    t_max_jd = t_min_jd+0.127163

    self.ztf_id = "ZTF21abyonuw"


    def test_query_ztfid(self):

        logger.info(f"Retrieving alerts for {self.ztf_id}")
        api_name = ampel_api_name(ztf_name=self.ztf_id, with_history=True, with_cutouts=False)
        self.assertEqual(
            len(api_name),
            1
        )

    def test_query_cone(self):

        logger.info("Commencing API cone search")
        api_cone = ampel_api_cone(ra=30, dec=30, radius=0.1)
        nr_transients = len(get_ztf_ids(api_cone))
        self.assertEqual(
            nr_transients,
            94
        )

    def test_query_time(self):
        logger.info("Commencing API time search")
        api_time = ampel_api_timerange(t_min_jd=t_min_jd, t_max_jd=t_max_jd, chunk_size=2000)
        nr_transients = len(get_ztf_ids(api_time))
        self.assertEqual(
            nr_transients,
            1887
        )

    def test_query_cutouts(self):
        logger.info(f"Retrieving cutouts for {self.ztf_id}")
        api_cutouts = ampel_api_cutout(candid=candid)
        nr_cutouts = len(api_cutouts)
        self.assertEqual(
            nr_cutouts,
            3
        )      

    @staticmethod
    def get_ztf_ids(query_result):
        ztf_ids = []
        for entry in query_result:
            ztf_ids.append(entry["objectId"])
        return list(set(ztf_ids))



