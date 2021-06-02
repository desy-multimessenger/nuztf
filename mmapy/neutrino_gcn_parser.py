#!/usr/bin/env python3

import time
import numpy as np
from neutrino_scanner import NeutrinoScanner

# nu_name = "IC200530A"
nu_name = "IC210510A"
n_days = 3

from astropy.time import Time

nu = NeutrinoScanner(nu_name=nu_name)

t_start = time.time()
nu.scan_cones()
t_end = time.time()
print(f"Querying took {t_end - t_start:.0f} seconds")

summary = nu.text_summary()
print(summary)

fig, message = nu.plot_overlap_with_observations(first_det_window_days=n_days)
fig.savefig("overlap_{}.png".format(nu_name), dpi=300)
print(nu.draft_gcn())
