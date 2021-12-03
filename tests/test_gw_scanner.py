import unittest
import logging
from nuztf.skymap_scanner import SkymapScanner


class TestGWScanner(unittest.TestCase):

    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

    def test_scan(self):

        self.logger.info("\n\n Testing GW Scanner \n\n")

        gw_name = "S190814bv"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GW {gw_name}")

        scanner = SkymapScanner(
            event_name=gw_name,
            scan_mode="gw",
            prob_threshold=prob_threshold,
            n_days=3,
            logger=self.logger,
        )

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 3474

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
