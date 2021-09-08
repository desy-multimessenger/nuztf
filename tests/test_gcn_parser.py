import unittest

from nuztf.parse_nu_gcn import get_latest_gcn, gcn_url, find_gcn_no
from ampel.log.AmpelLogger import AmpelLogger


logger = AmpelLogger()


class TestNeutrinoScanner(unittest.TestCase):

    maxDiff = None

    def test_latest(self):
        logger.info('\n\n Testing parsing of GCNs \n\n')
        no = get_latest_gcn()
        logger.info(f'\n\n Latest alert is {no} \n\n')
        url = gcn_url(gcn_number=no)

    def test_named(self):

        name = "IC200620A"

        num = int(find_gcn_no(name))

        ref = 27997

        logger.info(f'GCN number for {name} was found to be {num}')
        logger.info(f'Reference value was {ref}')

        self.assertEqual(num, ref)
