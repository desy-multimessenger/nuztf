#!/usr/bin/env python
# coding: utf-8

import logging

import numpy as np
from astropy.time import Time
from nuztf.skymap_scanner import SkymapScanner

logger = logging.getLogger("nuztf")
logger.setLevel(logging.INFO)

prob_threshold = 0.9
time_window_det_days = 3


def draw_random_map_name():
    max_index_maps = 8257
    rand_index = np.random.randint(0, high=max_index_maps + 1)
    random_map_name = f"O4_simul_{rand_index}.fits"
    return random_map_name


def draw_random_time():
    lower_time_bound_jd = 2458300.5  # July 1, 2018
    higher_time_bound_jd = 2459761.5  # June 31, 2022
    rand_time = np.random.uniform(low=lower_time_bound_jd, high=higher_time_bound_jd)
    random_time = Time(rand_time, format="jd")
    return random_time


scanner = SkymapScanner(
    event=draw_random_map_name(),
    prob_threshold=prob_threshold,
    n_days=time_window_det_days,
    t_obs=draw_random_time(),
)

scanner.get_alerts()

n_retrieved_alerts = scanner.n_alerts

scanner.filter_alerts()

scanner.create_overview_table()

scanner.create_candidate_summary()

fig, coverage_summary = scanner.plot_coverage()

tns_summary = scanner.tns_summary()
res = scanner.draft_gcn()
