#!/usr/bin/env python3

from neutrino_scanner import NeutrinoScanner
from ampel.log.AmpelLogger import AmpelLogger

logger = AmpelLogger()


nu_name = "IC-200916A"
n_days = 1

nu = NeutrinoScanner(nu_name, logger=logger)
nu.scan_cones()
fig, message = nu.plot_overlap_with_observations()
fig.savefig('overlap_{}.png'.format(nu_name), dpi=300)
print(nu.draft_gcn())
