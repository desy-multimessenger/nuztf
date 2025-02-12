import unittest
import logging
from pathlib import Path
from nuztf.cli import main


logger = logging.getLogger(__name__)


class TestCLI(unittest.TestCase):
    def test_cli(self):
        neutrino_name = "IC200620A"
        true_gcn = "Astronomer Name (Institute of Somewhere), ............. report,\n\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n\nAs part of the ZTF neutrino follow up program (Stein et al. 2023), we observed the localization region of the neutrino event IceCube-200620A (Santander et. al, GCN 27997) with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g- and r-band beginning at 2020-06-21 04:53 UTC, approximately 25.8 hours after event time. We covered 77.6% (1.3 sq deg) of the reported localization region. This estimate accounts for chip gaps. Each exposure was 300s with a typical depth of 21.0 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) . We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF18acvhwtf | AT2020ncs | 162.0677272 | +12.1263357 | g      | 19.74 | 0.16   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF20abgvabi | AT2020ncr | 162.5306341 | +12.1461187 | g      | 20.58 | 0.19   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \n\nZTF18acvhwtf was first detected on 2018-12-09. It has a spec-z of 0.291 [1500 Mpc] and an abs. mag of -21.1. Distance to SDSS galaxy is 0.06 arcsec. [MILLIQUAS: SDSS J104816.25+120734.7 - 'Q'-type source (0.06 arsec)] [TNS NAME=AT2020ncs]\nZTF20abgvabi was first detected on 2020-05-26. WISE DETECTION: W1-W2=0.04 (1.03 arsec) [TNS NAME=AT2020ncr]\n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the nuztf (Stein et al. 2021, https://github.com/desy-multimessenger/nuztf ).\n"
        tmpfile = Path("tmpfile.txt")
        main(
            nu_name=neutrino_name,
            logging_level="DEBUG",
            gcn_filename=tmpfile,
            rich_handler=False,
            stream_handler=True,
        )
        with tmpfile.open("r") as f:
            gcn = f.read()
        self.assertEqual(gcn, true_gcn)
        tmpfile.unlink()
