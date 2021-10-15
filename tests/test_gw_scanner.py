import unittest
import logging
from nuztf.skymap_scanner import SkymapScanner


class TestGWScanner(unittest.TestCase):

    maxDiff = None

    def test_scan(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)

        logger.info("\n\n Testing GW Scanner \n\n")

        gw_name = "S190426c"
        prob_threshold = 0.9

        logging.info(f"Scanning with GW {gw_name}")

        scanner = SkymapScanner(
            event_name=gw_name,
            scan_mode="gw",
            prob_threshold=prob_threshold,
            n_days=3,
            logger=logger,
        )

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 110332

        logging.info(
            f"Retrieved {n_retrieved_alerts} alerts. {n_expected_alerts} alerts expected."
        )

        self.assertEqual(n_retrieved_alerts, n_expected_alerts)

        scanner.filter_alerts()

        n_retrieved_candidates = len(scanner.final_candidates)
        n_expected_candidates = 41

        logging.info(
            f"Retrieved {n_retrieved_candidates} candidates. {n_expected_candidates} candidates expected."
        )

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)
