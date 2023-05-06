#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import logging
import os
from logging import Handler
from pathlib import Path

from astropy.time import Time  # type: ignore
from slack import WebClient  # type: ignore

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.skymap import EventNotFound
from nuztf.skymap_scanner import SkymapScanner

logging.basicConfig()


class SlackLogHandler(Handler):
    def __init__(self, channel: str, ts: str, webclient: WebClient):
        Handler.__init__(self)

        self.webclient = webclient
        self.channel = channel
        self.ts = ts

    def info(self, record) -> None:
        self.webclient.chat_postMessage(
            channel=self.channel,
            thread_ts=self.ts,
            text=record,
        )

    def debug(self, record) -> None:
        return None


class Slackbot:
    def __init__(
        self,
        channel: str,
        ts,
        name: str,
        event_type: str,
        dl_results: bool = False,
        do_gcn: bool = False,
        time_window: int | None = None,
    ):
        self.channel = channel
        self.ts = ts
        self.name = name
        self.event_type = event_type
        self.dl_results = dl_results
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.webclient = WebClient(token=os.environ.get("SLACK_TOKEN"))

        if time_window is None:
            if self.event_type == "nu":
                self.time_window = 10
            elif self.event_type == "gw":
                self.time_window = 2
        else:
            self.time_window = time_window

        self.scanner: NeutrinoScanner | SkymapScanner

        if self.event_type == "nu":
            self.scanner = NeutrinoScanner(self.name)
        elif self.event_type == "gw":
            try:
                self.scanner = SkymapScanner(self.name)
            except EventNotFound:
                self.post(
                    f"The specified LIGO event, {self.name}, was not found on GraceDB. Please check that you entered the correct event name."
                )
                return

        self.scanner.logger = SlackLogHandler(
            channel=self.channel, ts=self.ts, webclient=self.webclient
        )

        if self.event_type == "nu":
            self.scanner.scan_area(t_max=self.scanner.t_min + self.time_window)

        elif self.event_type == "gw":
            if self.dl_results:
                self.scanner.download_results()
            else:
                self.scanner.get_alerts()
                self.scanner.filter_alerts()

        scan_message = "Scanning done."

        if len(self.scanner.cache) > 0:
            scan_message += f" Found {len(self.scanner.cache)} candidates:\n\n"
            for entry in list(self.scanner.cache.keys()):
                scan_message += f"{entry}\n"
        else:
            scan_message += "\nNo candidates found."

        self.post(scan_message)

        if len(self.scanner.cache) > 0:
            pdf_overview_path = self.scanner.summary_path + ".pdf"
            self.post_file(pdf_overview_path, f"{self.name}_candidates.pdf")

            self.scanner.create_overview_table()
            csv_path = self.scanner.summary_path + ".csv"
            self.post_file(csv_path, f"{self.name}_candidates.csv")

        if do_gcn and len(self.scanner.cache) > 0:
            self.create_gcn()

    def create_gcn(self):
        self.scanner.plot_overlap_with_observations(
            first_det_window_days=self.time_window
        )
        self.post(f"*Observations*\n{self.scanner.observations}")
        gcn = self.scanner.draft_gcn()
        self.post(text=gcn)

    def post(self, text: str):
        """Post to Slack"""
        self.webclient.chat_postMessage(
            channel=self.channel,
            text=text,
            thread_ts=self.ts,
        )

    def post_file(self, filepath, filename):
        """Post a file to Slack"""
        with open(filepath, "rb") as file:
            self.webclient.files_upload(
                file=file,
                filename=filename,
                channels=self.channel,
                thread_ts=self.ts,
            )
