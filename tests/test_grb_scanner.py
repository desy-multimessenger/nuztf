import unittest
import logging
from nuztf.skymap_scanner import SkymapScanner


class TestGRBScanner(unittest.TestCase):

    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_scan(self):

        self.logger.info("\n\n Testing GRB Scanner \n\n")

        grb_name = "GRB210927A"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GRB {grb_name}")

        scanner = SkymapScanner(
            event_name=grb_name,
            scan_mode="grb",
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
