#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)

import logging

from astropy.time import Time  # type: ignore

from nuztf.neutrino_scanner import NeutrinoScanner
from nuztf.skymap_scanner import SkymapScanner

logging.basicConfig()


class Slackbot:
    def __init__(
        self,
        channel: str,
        name: str,
        event_type: str,
        do_gcn: bool = False,
        time_window: int | None = None,
    ):
        self.channel = channel
        self.name = name
        self.event_type = event_type
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

        if time_window is None:
            time_window = 10

        nu = NeutrinoScanner(self.name)
        nu.scan_area(t_max=nu.t_min + time_window)
        nu.peak_mag_summary()

        if do_gcn:
            nu.plot_overlap_with_observations(first_det_window_days=time_window)
            self.gcn = nu.draft_gcn()
