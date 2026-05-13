import logging
import unittest
from pathlib import Path

import requests
from typer.testing import CliRunner

from nuztf.cli import app
from nuztf.neutrino_kafka_scanner import save_alert
from tests.test_neutrino_scanner import EXAMPLE_KAFKA_ALERT

logger = logging.getLogger(__name__)


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
