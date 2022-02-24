import unittest
import logging
from nuztf.skymap_scanner import SkymapScanner


class TestSkymapScanner(unittest.TestCase):

    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_gw_scan(self):

        self.logger.info("\n\n Testing GW Scanning \n\n")

        gw_name = "S190814bv"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GW {gw_name}")

        scanner = SkymapScanner(
            event=gw_name,
            prob_threshold=prob_threshold,
            n_days=1,
            logger=self.logger,
        )

        scanner.plot_skymap()

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 1649

        self.logger.info(
            f"Retrieved {n_retrieved_alerts} alerts. {n_expected_alerts} alerts expected."
        )

        self.assertEqual(n_retrieved_alerts, n_expected_alerts)

        scanner.filter_alerts()

        n_retrieved_candidates = len(scanner.final_candidates)
        n_expected_candidates = 1

        self.logger.info(
            f"Retrieved {n_retrieved_candidates} candidates. {n_expected_candidates} candidates expected."
        )

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)

        scanner.create_overview_table()

        scanner.create_candidate_summary()

        fig, coverage_summary = scanner.plot_coverage()

        true_coverage_summary = "In total, 88.57 % of the contour was observed at least once.\nThis estimate includes 0.00 % of the contour at a galactic latitude <10 deg.\nIn total, 73.81 % of the contour was observed at least twice. \nIn total, 73.81 % of the contour was observed at least twice, and excluding low galactic latitudes.\nThese estimates account for chip gaps."

        self.assertEqual(coverage_summary, true_coverage_summary)

        tns_summary = scanner.tns_summary()

        true_tns_summary = "Candidate: ZTF19abpuhbh / RA=11.2569224 / Dec=-22.5161471 / First detection=2458710.9475231\nLast Upper Limit: None\nFirst Detection: 2458710.9475231 / band=r / mag=20.857 +/- 0.284\nFirst observed 13.56 hours after merger\n[2458710.9475231, 2458710.9948611]\n"

        self.assertEqual(tns_summary, true_tns_summary)

    def test_grb_scan(self):

        self.logger.info("\n\n Testing GRB Scanner \n\n")

        grb_name = "GRB210927A"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GRB {grb_name}")

        scanner = SkymapScanner(
            event=grb_name,
            prob_threshold=prob_threshold,
            n_days=0.1,
            logger=self.logger,
        )

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 6356

        self.logger.info(
            f"Retrieved {n_retrieved_alerts} alerts. {n_expected_alerts} alerts expected."
        )

        self.assertEqual(n_retrieved_alerts, n_expected_alerts)

        scanner.filter_alerts()

        n_retrieved_candidates = len(scanner.final_candidates)
        n_expected_candidates = 3

        self.logger.info(
            f"Retrieved {n_retrieved_candidates} candidates. {n_expected_candidates} candidates expected."
        )

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)
