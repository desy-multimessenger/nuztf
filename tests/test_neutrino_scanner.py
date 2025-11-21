import logging
import unittest
from collections import OrderedDict
from pathlib import Path

import healpy as hp
import ligo.skymap.plot as lkp
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from astropy.visualization.wcsaxes import Quadrangle
from matplotlib.patches import Polygon
from ztfquery.fields import get_field_vertices

from nuztf.neutrino_kafka_scanner import NeutrinoKafkaScanner
from nuztf.neutrino_scanner import NeutrinoScanner

EXAMPLE_KAFKA_ALERT = {
    "$schema": "https://gcn.nasa.gov/schema/v6.0.0/gcn/notices/icecube/single_neutrino_alerts.schema.json",
    "mission": "IceCube",
    "instrument": "IC86",
    "messenger": "Neutrino",
    "pipeline": "Bronze Track Alert",
    "record_number": 1,
    "event_name": ["IceCube-230416A"],
    "id": ["137840_57034692_0"],
    "alert_datetime": "2023-04-16T05:42:00.0Z",
    "alert_type": "initial",
    "alert_tense": "current",
    "alert_topology": "Track",
    "number_of_events": 1,
    "ra": 345.82,
    "dec": 9.01,
    "ra_dec_error": 0.5,
    "containment_probability": 0.9,
    "systematic_included": False,
    "healpix_url": "https://roc.icecube.wisc.edu/public/alerts/example/run00140078.evt000030891383.example.skymap_nside_1024_probability.fits.gz",
    "trigger_time": "2023-04-16T05:22:26.150574Z",
    "nu_energy": 127.29,
    "p_astro": 0.34064,
    "far": 8.029e-8,
}


class TestNeutrinoScanner(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.max_distance_diff_arcsec = 2
        self.maxDiff = None

        self.neutrino_name = "IC200620A"
        self.expected_candidates = 2

    @pytest.fixture(autouse=True)
    def initdir(self, tmp_path):
        self.tmp_path = tmp_path

    def test_classical_scan(self):
        self.logger.info("\n\n Testing Neutrino Scanner \n\n")

        self.logger.info(f"scanning with neutrino {self.neutrino_name}")
        nu = NeutrinoScanner(nu_name=self.neutrino_name)
        self.scan(nu)

    def test_kafka_scan(self):
        self.logger.info("\n\n Testing Neutrino Kafka Scanner \n\n")

        # using classical scanner to get GCN infos
        gcn_info = {
            "author": "Santander",
            "name": "IceCube-200620A",
            "time": Time("2020-06-20T03:03:32.280"),
            "ra": 162.11,
            "ra_err": [0.64, -0.95],
            "dec": 11.95,
            "dec_err": [0.63, -0.48],
        }
        nu = NeutrinoScanner(nu_name=self.neutrino_name)

        # make a mock skymap
        nside = 1024
        npix = hp.nside2npix(nside)

        ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)

        # -------- rectangle bounds --------
        ra_min = gcn_info["ra"] + gcn_info["ra_err"][1]
        ra_max = gcn_info["ra"] + gcn_info["ra_err"][0]
        dec_min = gcn_info["dec"] + gcn_info["dec_err"][1]
        dec_max = gcn_info["dec"] + gcn_info["dec_err"][0]

        # -------- mask pixels inside rectangle --------
        # RA may wrap across 0/360 → handle both cases
        if ra_min <= ra_max:
            inside_ra = (ra >= ra_min) & (ra <= ra_max)
        else:
            inside_ra = (ra >= ra_min) | (ra <= ra_max)

        inside_dec = (dec >= dec_min) & (dec <= dec_max)

        mask = inside_ra & inside_dec
        skymap = np.zeros(npix)
        num_inside = np.sum(mask)
        num_outside = npix - num_inside
        skymap[mask] = 0.9 / num_inside
        skymap[~mask] = 0.1 / num_outside

        # normalize (just in case)
        skymap /= skymap.sum()

        meta = OrderedDict(
            [
                ("PIXTYPE", "HEALPIX"),
                ("ORDERING", "RING"),
                ("COORDSYS", "C"),
                ("EXTNAME", "xtension"),
                ("NSIDE", nside),
                ("FIRSTPIX", 0),
                ("LASTPIX", 12582911),
                ("INDXSCHM", "IMPLICIT"),
                ("OBJECT", "FULLSKY"),
                ("RUNID", 140078),
                ("EVENTID", 30891383),
                ("SENDER", "IceCube Collaboration"),
                ("DATE-OBS", "2024-11-13T00:22:20.682"),
                ("MJD-OBS", nu.t_min.mjd),
                ("I3TYPE", "EHE"),
                ("RA", gcn_info["ra"]),
                ("DEC", gcn_info["dec"]),
                ("RA_ERR_PLUS", gcn_info["ra_err"][0]),
                ("RA_ERR_MINUS", gcn_info["ra_err"][1]),
                ("DEC_ERR_PLUS", gcn_info["dec_err"][0]),
                ("DEC_ERR_MINUS", gcn_info["dec_err"][1]),
                (
                    "COMMENTS",
                    "90% uncertainty location => Highest posterior density 90% credible region",
                ),
                (
                    "NOTE",
                    "Please ignore pixels with infinite or NaN values. They are rare cases of the minimizer failing to converge",
                ),
            ]
        )

        alert = dict(EXAMPLE_KAFKA_ALERT)
        alert["trigger_time"] = gcn_info["time"].isot
        alert["event_name"] = ["IceCube-" + self.neutrino_name.replace("IC", "")]

        filename = self.tmp_path / "skymap.fits"
        Table([skymap], meta=meta, names=["PROB"], units=["1 / pix"]).write(
            filename, format="fits"
        )
        nu_kafka = NeutrinoKafkaScanner(alert=alert, map_path=filename)

        assert len(nu_kafka.pixel_nos) == len(nu.pixel_nos)
        assert (nu_kafka.pixel_nos == nu.pixel_nos).all()

        # fake GCN circular info for compatibility
        nu_kafka.author = "Santander"
        nu_kafka.gcn_no = 27997

        center_ra = meta["RA"] * u.deg
        center_dec = meta["DEC"] * u.deg
        ra_err_minus = meta["RA_ERR_MINUS"] * u.deg
        ra_err_plus = meta["RA_ERR_PLUS"] * u.deg
        dec_err_minus = meta["DEC_ERR_MINUS"] * u.deg
        dec_err_plus = meta["DEC_ERR_PLUS"] * u.deg
        left_lower_corner_ra = center_ra + ra_err_minus
        left_lower_corner_dec = center_dec + dec_err_minus
        dra = ra_err_plus - ra_err_minus
        ddec = dec_err_plus - dec_err_minus

        fig = plt.figure()
        ax = plt.axes(
            projection="astro degrees zoom",
            center=f"{gcn_info['ra']}d {gcn_info['dec']}d",
            radius="2 deg",
        )
        _t = ax.get_transform("world")
        ax.imshow_hpx(str(filename), label="skymap")
        ax.contour_hpx(nu_kafka.credible_levels, levels=[0.9], colors="C0", nested=True)
        ax.plot([], [], color="C0", label="90% contour from skymap", ls="-")
        gcn_rect = Quadrangle(
            [left_lower_corner_ra, left_lower_corner_dec],
            dra,
            ddec,
            edgecolor="C1",
            facecolor="none",
            transform=_t,
            ls="--",
        )
        ax.add_patch(gcn_rect)
        ax.plot([], [], color="C1", label="90% rectangle from GCN", ls="--")

        fields = [522]
        verts = np.squeeze(get_field_vertices(fields, squeeze=False), axis=1)
        for i, vert in enumerate(verts):
            c = f"C{2 + i}"
            ax.add_patch(
                Polygon(vert, transform=_t, edgecolor=c, facecolor=c, ls="-", alpha=0.5)
            )
            ax.plot([], [], color=c, label=f"Field {fields[i]}")
        ax.legend(loc="upper right", ncol=2)

        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
        fig_fn = self.tmp_path / "skymap.pdf"
        self.logger.info(f"Saving skymap figure to {fig_fn}")
        fig.savefig(fig_fn)
        plt.close()

        assert nu_kafka.cone_nside == nu.cone_nside
        assert nu_kafka.cone_ids == nu.cone_ids
        self.scan(scanner=nu_kafka)

    def scan(self, scanner: NeutrinoScanner | NeutrinoKafkaScanner):
        t_max = scanner.default_t_max - 8

        scanner.scan_area(t_max=t_max)
        retrieved_candidates = scanner.n_candidates

        self.logger.info(
            f"found {retrieved_candidates}, expected {self.expected_candidates}"
        )
        self.assertEqual(self.expected_candidates, retrieved_candidates)

        for name, res in sorted(scanner.cache_candidates.items()):
            # Only use old data, so new detections do not change CI
            dets = [
                x
                for x in res["prv_candidates"]
                if ("isdiffpos" in x.keys()) & (x["jd"] < scanner.default_t_max.jd)
            ]
            dets = [
                x for x in dets if str(x["isdiffpos"]).lower() in ["t", "true", "1"]
            ]
            cand = dets[-1]
            res["candidate"] = cand
            res["prv_candidates"] = dets[:-1]

        scanner.plot_overlap_with_observations(
            first_det_window_days=(t_max - scanner.t_min).to("d").value
        )
        res = scanner.draft_gcn()

        print(repr(res))

        # Update the true using repr(res)
        true_gcn = "Astronomer Name (Institute of Somewhere), ............. report,\n\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n\nAs part of the ZTF neutrino follow up program (Stein et al. 2023), we observed the localization region of the neutrino event IceCube-200620A (Santander et. al, GCN 27997) with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g- and r-band beginning at 2020-06-21 04:53 UTC, approximately 25.8 hours after event time. We covered 77.6% (1.3 sq deg) of the reported localization region. This estimate accounts for chip gaps. Each exposure was 300s with a typical depth of 21.0 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) . We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF18acvhwtf | AT2020ncs | 162.0678742 | +12.1264130 | r      | 20.55 | 0.11   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF20abgvabi | AT2020ncr | 162.5306820 | +12.1462203 | r      | 20.67 | 0.10   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \n\nZTF18acvhwtf was first detected on 2018-12-09. It has a spec-z of 0.291 [1500 Mpc] and an abs. mag of -20.4. Distance to SDSS galaxy is 0.09 arcsec. [MILLIQUAS: SDSS J104816.25+120734.7 - 'Q'-type source (0.06 arsec)] [TNS NAME=AT2020ncs]\nZTF20abgvabi was first detected on 2020-05-26. WISE DETECTION: W1-W2=0.04 (1.03 arsec) [TNS NAME=AT2020ncr]\n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the nuztf (Stein et al. 2021, https://github.com/desy-multimessenger/nuztf ).\n"

        self.assertEqual(res, true_gcn)

        # Test manually adding candidates

        scanner.add_to_cache_by_names(
            ztf_ids=["ZTF18abteipt"],
        )

        # Check

        false_candidate = scanner.check_ampel_filter("ZTF18abteipt")

        self.logger.info(
            f"For the false candidate, the pipeline bool is {false_candidate}"
        )

        self.assertFalse(false_candidate)

        true_candidate = scanner.check_ampel_filter("ZTF20abgvabi")

        self.logger.info(
            f"For the true candidate, the pipeline bool is {true_candidate}"
        )

        self.assertTrue(true_candidate)
