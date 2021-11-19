#!/usr/bin/env python
# coding: utf-8
import time
import logging
import matplotlib.pyplot as plt
from nuztf.skymap_scanner import SkymapScanner


if __name__ == "__main__":

    event_name = "GRB210927A"
    prob_threshold = 0.9

    start = time.time()

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    scanner = SkymapScanner(
        event_name=event_name,
        # skymap_file=skymap_file,
        scan_mode="grb",
        prob_threshold=prob_threshold,
        n_days=3,
        logger=logger,
    )

    scanner.plot_skymap()
    scanner.get_alerts()
    scanner.filter_alerts()  # load_cachefile=True)
    scanner.create_overview_table()
    scanner.create_candidate_summary()
    scanner.plot_coverage()

    end = time.time()

    logger.info(f"The skript took {end-start} s.")
