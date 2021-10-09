import unittest
import logging
from nuztf import NeutrinoScanner


class TestNeutrinoScanner(unittest.TestCase):

    maxDiff = None

    def test_scan(self):
        logger = logging.getlogger(__name__)
        logger.setLevel(logging.DEBUG)

        logger.info('\n\n Testing Neutrino Scanner \n\n')
        
        neutrino_name = "IC200620A"
        expected_candidates = 2

        logging.info(f'scanning with neutrino {name}')
        nu = NeutrinoScanner(nu_name=neutrino_name, logger=logger)

        t_max = nu.default_t_max - 8

        nu.scan_cones(t_max=t_max)
        retrieved_candidates = len(nu.cache)

        logging.info(f"found {retrieved_candidates}, expected {expected_candidates}")
        self.assertEqual(expected_candidates, retrieved_candidates)

        nu.plot_overlap_with_observations(
            first_det_window_days=(t_max - nu.t_min).to("d").value
        )
        res = nu.draft_gcn()

        # Update the true using repr(res)
        true_gcn = "Astronomer Name (Institute of Somewhere), ............. report,\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \nWe observed the localization region of the neutrino event IceCube-200620A (Santander et. al, GCN 27997) with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g-band and r-band beginning at 2020-06-21 04:53 UTC, approximately 25.8 hours after event time. We covered 1.2 sq deg, corresponding to 77.7% of the reported localization region. This estimate accounts for chip gaps. Each exposure was 300s with a typical depth of 21.0 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) . We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF18acvhwtf | AT2020ncs | 162.0678527 | +12.1263986 | r      | 20.11 | 0.16   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF20abgvabi | AT2020ncr | 162.5306341 | +12.1461187 | g      | 20.58 | 0.19   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \nZTF18acvhwtf was first detected on 2458461.9815278. It has a spec-z of 0.291 [1548 Mpc] and an abs. mag of -20.8. Distance to SDSS galaxy is 0.52 arcsec. [MILLIQUAS: SDSS J104816.25+120734.7 - 'Q'-type source (0.55 arsec)]\nZTF20abgvabi was first detected on 2458995.6705903. \n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the AMPEL Follow-up Pipeline (Stein et al. 2021).\n"

        self.assertEqual(
            res,
            true_gcn
        )

        # Test manually adding candidates

        nu.add_to_cache_by_names("ZTF18abteipt",)

        # Check

        false_candidate = nu.check_ampel_filter("ZTF18abteipt")

        logging.info(f"For the false candidate, the pipeline bool is {false_candidate}")

        self.assertFalse(false_candidate)

        true_candidate = nu.check_ampel_filter("ZTF20abgvabi")

        logging.info(f"For the true candidate, the pipeline bool is {true_candidate}")

        self.assertTrue(true_candidate)
