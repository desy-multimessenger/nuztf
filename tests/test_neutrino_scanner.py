import logging
import unittest

import numpy as np
from astropy.coordinates import Distance
from astropy.time import Time

from nuztf.ampel_api import ampel_api_catalog
from nuztf.base_scanner import cosmo
from nuztf.cat_match import query_ned_for_z
from nuztf.neutrino_scanner import NeutrinoScanner


class TestNeutrinoScanner(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.max_distance_diff_arcsec = 2

    def test_scan(self):
        self.logger.info("\n\n Testing Neutrino Scanner \n\n")

        neutrino_name = "IC200620A"
        expected_candidates = 2

        self.logger.info(f"scanning with neutrino {neutrino_name}")
        nu = NeutrinoScanner(nu_name=neutrino_name)

        t_max = nu.default_t_max - 8

        # nu.scan_cones(t_max=t_max)
        nu.scan_area(t_max=t_max)
        retrieved_candidates = len(nu.cache)

        self.logger.info(
            f"found {retrieved_candidates}, expected {expected_candidates}"
        )
        self.assertEqual(expected_candidates, retrieved_candidates)

        for name, res in sorted(nu.cache.items()):
            # Only use old data, so new detections do not change CI
            dets = [
                x
                for x in res["prv_candidates"]
                if ("isdiffpos" in x.keys()) & (x["jd"] < nu.default_t_max.jd)
            ]
            cand = dets[-1]
            res["candidate"] = cand
            res["prv_candidates"] = dets[:-1]

        nu.plot_overlap_with_observations(
            first_det_window_days=(t_max - nu.t_min).to("d").value
        )
        res = nu.draft_gcn()

        print(repr(res))

        # Update the true using repr(res)
        true_gcn = "Astronomer Name (Institute of Somewhere), ............. report,\n\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n\nAs part of the ZTF neutrino follow up program (Stein et al. 2023), we observed the localization region of the neutrino event IceCube-200620A (Santander et. al, GCN 27997) with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g- and r-band beginning at 2020-06-21 04:53 UTC, approximately 25.8 hours after event time. We covered 77.6% (1.3 sq deg) of the reported localization region. This estimate accounts for chip gaps. Each exposure was 300s with a typical depth of 21.0 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) . We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF18acvhwtf | AT2020ncs | 162.0678742 | +12.1264130 | r      | 20.55 | 0.11   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF20abgvabi | AT2020ncr | 162.5306820 | +12.1462203 | r      | 20.67 | 0.10   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \n\nZTF18acvhwtf was first detected on 2018-12-09. It has a spec-z of 0.291 [1500 Mpc] and an abs. mag of -20.7. Distance to SDSS galaxy is 0.39 arcsec. [MILLIQUAS: SDSS J104816.25+120734.7 - 'Q'-type source (0.06 arsec)] [TNS NAME=AT2020ncs]\nZTF20abgvabi was first detected on 2020-05-26. WISE DETECTION: W1-W2=0.04 (1.03 arsec) [TNS NAME=AT2020ncr]\n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the nuztf (Stein et al. 2021, https://github.com/desy-multimessenger/nuztf ).\n"

        self.assertEqual(res, true_gcn)

        # Test manually adding candidates

        nu.add_to_cache_by_names(
            ztf_ids=["ZTF18abteipt"],
        )

        # Check

        false_candidate = nu.check_ampel_filter("ZTF18abteipt")

        self.logger.info(
            f"For the false candidate, the pipeline bool is {false_candidate}"
        )

        self.assertFalse(false_candidate)

        true_candidate = nu.check_ampel_filter("ZTF20abgvabi")

        self.logger.info(
            f"For the true candidate, the pipeline bool is {true_candidate}"
        )

        self.assertTrue(true_candidate)
