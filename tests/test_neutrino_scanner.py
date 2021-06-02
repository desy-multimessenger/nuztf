import unittest

from neutrino_scanner import NeutrinoScanner
from ampel.log.AmpelLogger import AmpelLogger


logger = AmpelLogger()


class TestNeutrinoScanner(unittest.TestCase):

    def test_scan(self):
        logger.info('\n\n Testing Neutrino Scanner \n\n')
        name = "IC210510A"
        expected_candidates = 6

        logger.info(f'scanning with neutrino {name}')
        nu = NeutrinoScanner(name, logger=logger)
        nu.scan_cones()
        retrieved_candidates = len(nu.cache)

        logger.info(f"found {retrieved_candidates}, expected {expected_candidates}")
        self.assertEqual(expected_candidates, retrieved_candidates)