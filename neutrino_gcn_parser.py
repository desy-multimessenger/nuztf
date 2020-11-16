#!/usr/bin/env python3

import pickle
import numpy as np
from neutrino_scanner import NeutrinoScanner
from ampel.log.AmpelLogger import AmpelLogger

logger = AmpelLogger()

nu_name = "IC201014A"
n_days = 4

# nu = NeutrinoScanner(nu_name=nu_name, logger=logger)
# nu.scan_cones()
# fig, message = nu.plot_overlap_with_observations()
# fig.savefig('overlap_{}.png'.format(nu_name), dpi=300)
# print(nu.draft_gcn())

from astropy.time import Time
nu = NeutrinoScanner(
	# manual_args=(
	#     "IC-201007A", 
	#     [265.17, +0.52, -0.52],
	#     [5.34, +0.32, -0.23],
	#     Time(
	#     	"2020-10-07T22:01:49.28", 
	#     	format='isot', 
	#     	scale='utc'
	#     ),
	# ),
	nu_name=nu_name,
	logger=logger,
)
nu.scan_cones()
fig, message = nu.plot_overlap_with_observations(first_det_window_days=n_days)
fig.savefig('overlap_{}.png'.format(nu_name), dpi=300)
print(nu.draft_gcn())
