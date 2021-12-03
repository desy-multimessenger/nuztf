import unittest
import logging
from nuztf.irsa import plot_irsa_lightcurve


class TestIrsa(unittest.TestCase):

    def test_lightcurve(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)

        logger.info("\n\n Testing IRSA \n\n")

        src_name = "PKS1502+106"

        plot_irsa_lightcurve(
            source_name=src_name,
            nu_name="IC190730A",
            check_obs=False,
            # check_obs_lookback_weeks=12,
        )
