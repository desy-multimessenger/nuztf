import unittest

from nuztf.parse_nu_gcn import get_latest_gcn, gcn_url, find_gcn_no, ParsingError
from ampel.log.AmpelLogger import AmpelLogger


logger = AmpelLogger()


class TestNeutrinoScanner(unittest.TestCase):

    maxDiff = None

    def test_latest(self):
        logger.info("\n\n Testing parsing of GCNs \n\n")
        no = get_latest_gcn()
        logger.info(f"Latest alert is {no}")
        url = gcn_url(gcn_number=no)
        logger.info(f"URL is {url}")

    def test_named(self):

        name = "IC200620A"

        num = int(find_gcn_no(name))

        ref = 27997

        logger.info(f"GCN number for {name} was found to be {num}")
        logger.info(f"Reference value was {ref}")

        self.assertEqual(num, ref)

        fakename = "IC130921A"

        logger.info(f"Searching for fictional alert {fakename}")

        try:
            no = find_gcn_no(fakename)
            raise Exception(
                f"Somehow found a GCN ({no}) matching "
                f"fictional neutrino alert {fakename}"
            )
        except ParsingError:
            logger.info("No GCN found, as expected.")
            pass
