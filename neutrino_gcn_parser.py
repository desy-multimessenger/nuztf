#!/usr/bin/env python3

import pickle
import backoff, time
import numpy as np
from neutrino_scanner import NeutrinoScanner
from ampel.log.AmpelLogger import AmpelLogger

logger = AmpelLogger()

# nu_name = "IC200530A"
nu_name = "IC210510A"
n_days = 3

# nu = NeutrinoScanner(nu_name=nu_name, logger=logger)
# nu.scan_cones()
# fig, message = nu.plot_overlap_with_observations()
# fig.savefig('overlap_{}.png'.format(nu_name), dpi=300)
# print(nu.draft_gcn())

from astropy.time import Time

nu = NeutrinoScanner(
    # manual_args=(
    #     "IC-201007A",
    #     [265.17, +0.01, -0.01],
    #     [5.34, +0.01, -0.01],
    #     Time(
    #       "2020-10-07T22:01:49.28",
    #       format='isot',
    #       scale='utc'
    #     ),
    # ),
    nu_name=nu_name,
    logger=logger,
)

# nu.add_to_cache_by_names("ZTF19aatubsj")
# lol = nu.get_avro_by_name("ZTF19aatubsj")
# import json
# with open('/Users/simeon/Desktop/avro_api.json', 'w') as outfile:
#     json.dump(lol, outfile)

# quit()
t_start = time.time()
nu.scan_cones()
t_end = time.time()
print(t_end - t_start)

quit()

fig, message = nu.plot_overlap_with_observations(first_det_window_days=n_days)
fig.savefig("overlap_{}.png".format(nu_name), dpi=300)
print(nu.draft_gcn())
