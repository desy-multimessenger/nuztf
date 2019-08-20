#!/usr/bin/env python
# coding: utf-8

from ampel.ztf.archive.ArchiveDB import ArchiveDB
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy as scp
import datetime
import ztfquery
import datetime
import re
from ztfquery import alert
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import csv
import os,io
import pickle
from astropy.coordinates import SkyCoord
from astropy import units as u
import getpass
import psycopg2
import sqlalchemy
import healpy as hp
from tqdm import tqdm
import zerorpc
from ampel.contrib.hu.t0.DecentFilter import DecentFilter
from ampel.pipeline.t0.DevAlertProcessor import DevAlertProcessor
from ampel.base.AmpelAlert import AmpelAlert


try:
    with open(".AMPEL_user.txt", "r") as f:
        username = f.read()
except FileNotFoundError:
    username = getpass.getpass(prompt='Username: ', stream=None)
    with open(".AMPEL_user.txt", "wb") as f:
        f.write(username.encode())
        
try:
    with open(".AMPEL_pass.txt", "r") as f:
        password = f.read()
except FileNotFoundError:
    password = getpass.getpass(prompt='Password: ', stream=None)
    with open(".AMPEL_pass.txt", "wb") as f:
        f.write(password.encode())


try:
    ampel_client = ArchiveDB('postgresql://{0}:{1}@localhost:5432/ztfarchive'.format(username, password))
except sqlalchemy.exc.OperationalError as e:
    print("You can't access the archive database without first opening the port.")
    print("Open a new terminal, and into that terminal, run the following command:")
    print("ssh -L5432:localhost:5433 ztf-wgs.zeuthen.desy.de")
    print("If that command doesn't work, you are either not a desy user or you have a problem in your ssh config.")
    raise e

class AmpelWizard:

    def __init__(self, run_config, t_min, logger=None, base_config=None, filter_class=DecentFilter, cone_nside=64,):
        self.cone_nside = cone_nside
        self.t_min = t_min

        if base_config is None:
            base_config = {'catsHTM.default': "tcp://127.0.0.1:27020"}

        self.ampel_filter_class = filter_class(set(), base_config=base_config,
                                               run_config=filter_class.RunConfig(**run_config),
                                               logger=logger)
        self.dap = DevAlertProcessor(self.ampel_filter_class)

        if not hasattr(self, "output_path"):
            self.output_path = None

        self.scanned_pixels = []
        self.cone_ids, self.cone_coords = self.find_cone_coords()
        self.cache = dict()
        self.default_t_max = Time.now()

    def filter_ampel(self, res):
        return self.ampel_filter_class.apply(AmpelAlert(res['objectId'], *self.dap._shape(res))) is not None

    def scan_cones(self, t_max=None, max_cones=None):

        if max_cones is None:
            max_cones = len(self.cone_ids)

        scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))
        print("Commencing Ampel queries!")
        print("Scan radius is", scan_radius)
        print("So far, {0} pixels out of {1} have already been scanned.".format(
            len(self.scanned_pixels), len(self.cone_ids)
        ))

        for i, cone_id in enumerate(tqdm(list(self.cone_ids)[:max_cones])):
            ra, dec = self.cone_coords[i]

            if cone_id not in self.scanned_pixels:
                res = self.query_ampel(np.degrees(ra), np.degrees(dec), scan_radius, t_max)
                #                 print(len(res))

                for res_alert in res:

                    if res_alert['objectId'] not in self.cache.keys():
                        self.cache[res_alert['objectId']] = res_alert
                    elif res_alert['candidate']["jd"] > self.cache[res_alert['objectId']]['candidate']["jd"]:
                        self.cache[res_alert['objectId']] = res_alert
                self.scanned_pixels.append(cone_id)

        print("Scanned {0} pixels".format(len(self.scanned_pixels)))
        print("Found {0} candidates".format(len(self.cache)))

        self.create_candidate_summary()

    def filter_f_no_prv(self, res):
        raise NotImplementedError

    def filter_f_history(self, res):
        raise NotImplementedError

    def find_cone_coords(self):
        raise NotImplementedError

    def query_ampel(self, ra, dec, rad, t_max=None):

        if t_max is None:
            t_max = self.default_t_max

        ztf_object = ampel_client.get_alerts_in_cone(
            ra, dec, rad, self.t_min.jd, t_max.jd, with_history=False)
        query_res = [i for i in ztf_object]

        candids = []
        for res in query_res:
            if self.filter_f_no_prv(res):
                if self.filter_ampel(res) is not None:
                    candids.append(res["candid"])

        ztf_object = ampel_client.get_alerts(candids, with_history=True)
        query_res = [i for i in ztf_object]
        final_res = []

        for res in query_res:
            if self.filter_f_history(res):
                final_res.append(res)

        return final_res

    @staticmethod
    def reassemble_alert(mock_alert):
        cutouts = ampel_client.get_cutout(mock_alert["candid"])
        for k in cutouts:
            mock_alert['cutout{}'.format(k.title())] = {'stampData': cutouts[k], 'fileName': 'dunno'}
        mock_alert['schemavsn'] = 'dunno'
        mock_alert['publisher'] = 'dunno'
        for pp in [mock_alert['candidate']] + mock_alert['prv_candidates']:
            # if pp['isdiffpos'] is not None:
            # pp['isdiffpos'] = ['f', 't'][pp['isdiffpos']]
            pp['pdiffimfilename'] = 'dunno'
            pp['programpi'] = 'dunno'
            pp['ssnamenr'] = 'dunno'

        return mock_alert

    @staticmethod
    def extract_ra_dec(nside, index):
        (colat, ra) = hp.pix2ang(nside, index, nest=True)
        dec = np.pi / 2. - colat
        return (ra, dec)

    @staticmethod
    def extract_npix(nside, ra, dec):
        colat = np.pi / 2. - dec
        return hp.ang2pix(nside, colat, ra, nest=True)

    def create_candidate_summary(self):

        print("Saving to:", self.output_path)

        with PdfPages(self.output_path) as pdf:
            for (name, old_alert) in tqdm(sorted(self.cache.items())):
                mock_alert = self.reassemble_alert(old_alert)
                fig = alert.display_alert(mock_alert)
                fig.text(0, 0, name)
                pdf.savefig()
                plt.close()