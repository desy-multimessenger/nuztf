#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import logging
import os

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

        self.post("Scan complete")

        self.nu.peak_mag_summary()

        self.post(self.nu.observations)

        if do_gcn:
            self.post("Creating GCN (obtaining and unpacking observations)")
            self.create_gcn()
            # self.nu.plot_overlap_with_observations(first_det_window_days=time_window)
            # self.gcn = nu.draft_gcn()
            # self.webclient.chat_postMessage(
            #     channel=channel,
            #     text=slack_bot.gcn,
            #     thread_ts=ts,
            # )

    def create_gcn(self):
        self.nu.plot_overlap_with_observations(first_det_window_days=self.time_window)
        gcn = self.nu.draft_gcn()
        self.post(text=gcn)

    def post(self, text: str):
        """Post to Slack"""
        self.webclient.chat_postMessage(
            channel=self.channel,
            text=text,
            thread_ts=self.ts,
        )
