#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import logging
import os
from pathlib import Path

from astropy.time import Time  # type: ignore

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.skymap_scanner import SkymapScanner
from slack import WebClient  # type: ignore

logging.basicConfig()


class Slackbot:
    def __init__(
        self,
        channel: str,
        ts,
        name: str,
        event_type: str,
        do_gcn: bool = False,
        time_window: int | None = None,
    ):
        self.channel = channel
        self.ts = ts
        self.name = name
        self.event_type = event_type
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.webclient = WebClient(token=os.environ.get("SLACK_TOKEN"))

        if time_window is None:
            self.time_window = 10
        else:
            self.time_window = time_window

        self.nu = NeutrinoScanner(self.name)
        self.nu.scan_area(t_max=self.nu.t_min + self.time_window)

        scan_message = "Scanning done."
        if len(self.nu.cache) > 0:
            scan_message += f" Found {len(self.nu.cache)} candidates:\n\n"
            for entry in list(self.nu.cache.keys()):
                scan_message += f"{entry}\n"
        else:
            scan_message += "\nNo candidates found."

        self.post(scan_message)

        if len(self.nu.cache) > 0:
            pdf_overview_path = self.nu.summary_path + ".pdf"
            self.post_file(pdf_overview_path, f"{self.name} candidates")

        if do_gcn and len(self.nu.cache) > 0:
            self.create_gcn()

    def create_gcn(self):
        self.nu.plot_overlap_with_observations(first_det_window_days=self.time_window)
        self.post(f"observations:\n{self.nu.observations}")
        gcn = self.nu.draft_gcn()
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
