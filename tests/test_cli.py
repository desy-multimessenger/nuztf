import logging
import unittest
from pathlib import Path

import requests
from typer.testing import CliRunner

from nuztf.cli import app
from nuztf.neutrino_kafka_scanner import save_alert
from tests.test_neutrino_scanner import EXAMPLE_KAFKA_ALERT

logger = logging.getLogger(__name__)


# EXAMPLE_KAFKA_ALERT = {
#     "$schema": "https://gcn.nasa.gov/schema/v6.0.0/gcn/notices/icecube/single_neutrino_alerts.schema.json",
#     "mission": "IceCube",
#     "instrument": "IC86",
#     "messenger": "Neutrino",
#     "pipeline": "Bronze Track Alert",
#     "record_number": 1,
#     "event_name": ["IceCube-230416A"],
#     "id": ["137840_57034692_0"],
#     "alert_datetime": "2023-04-16T05:42:00.0Z",
#     "alert_type": "initial",
#     "alert_tense": "current",
#     "alert_topology": "Track",
#     "number_of_events": 1,
#     "ra": 345.82,
#     "dec": 9.01,
#     "ra_dec_error": 0.5,
#     "containment_probability": 0.9,
#     "systematic_included": False,
#     "healpix_url": "https://roc.icecube.wisc.edu/public/alerts/example/run00140078.evt000030891383.example.skymap_nside_1024_probability.fits.gz_error",
#     "trigger_time": "2023-04-16T05:22:26.150574Z",
#     "nu_energy": 127.29,
#     "p_astro": 0.34064,
#     "far": 8.029e-8,
# }


class TestCLI(unittest.TestCase):
    def test_run_classic(self):
        # Test only creation of NeutrinoScanner and GCN draft generation
        # The scanner itself is tested elsewhere
        runner = CliRunner()
        neutrino_name = "IC200620A"
        tmpfile = Path("tmpfile.txt")
        res = runner.invoke(app, ["nu-classic", neutrino_name, "-f", str(tmpfile)])
        if res.exit_code != 0:
            raise res.exc_info[1]
        assert tmpfile.exists(), f"No GCN draft at {tmpfile}!"
        tmpfile.unlink()

    def test_run_saved_kafka(self):
        # this will only test entry into NeutrinoKafkaScanner
        # The scanner itself is tested elsewhere
        save_alert(EXAMPLE_KAFKA_ALERT, overwrite=True)
        runner = CliRunner()
        tmpfile = Path("tmpfile.txt")
        with self.assertRaisesRegex(
            requests.HTTPError, "404 Client Error: Not Found for url"
        ):
            res = runner.invoke(
                app,
                [
                    "nu-saved-kafka",
                    EXAMPLE_KAFKA_ALERT["event_name"][0],
                    "-f",
                    str(tmpfile),
                ],
            )
            if res.exit_code != 0:
                raise res.exc_info[1]
