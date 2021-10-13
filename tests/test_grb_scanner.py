import unittest
import logging
from nuztf.skymap_scanner import SkymapScanner


class TestNeutrinoScanner(unittest.TestCase):

    maxDiff = None

    def test_scan(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)

        logger.info("\n\n Testing GRB Scanner \n\n")

        grb_name = "GRB210927A"
        prob_threshold = 0.9

        logging.info(f"Scanning with GRB {grb_name}")

        scanner = SkymapScanner(
            event_name=grb_name,
            scan_mode="grb",
            prob_threshold=prob_threshold,
            n_days=3,
            logger=logger,
        )

        scanner.get_alerts()

        retrieved_alerts = scanner.n_alerts
        expected_alerts = 35130

        logging.info(
            f"Retrieved {retrieved_alerts} alerts. {expected_alerts} alerts expected."
        )

        self.assertEqual(retrieved_alerts, expected_alerts)

        scanner.filter_alerts()

        retrieved_candidates = scanner.final_candidates
        expected_candidates = 98

        logging.info(
            f"Retrieved {retrieved_candidates} candidates. {expected_candidates} candidates expected."
        )

        self.assertEqual(retrieved_candidates, expected_candidates)
