import unittest
import logging
import numpy as np

from astropy.coordinates import Distance
from astropy.time import Time

from nuztf import NeutrinoScanner
from nuztf.cat_match import query_ned_for_z
from nuztf.ampel_api import ampel_api_catalog
from nuztf.base_scanner import cosmo


class TestNeutrinoScanner(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.max_distance_diff_arcsec = 2

    def test_scan(self):
        self.logger.info("\n\n Testing Neutrino Scanner \n\n")

        neutrino_name = "IC200620A"
        expected_candidates = 2

        self.logger.info(f"scanning with neutrino {neutrino_name}")
        nu = NeutrinoScanner(nu_name=neutrino_name)

        t_max = nu.default_t_max - 8

        # nu.scan_cones(t_max=t_max)
        nu.scan_area(t_max=t_max)
        retrieved_candidates = len(nu.cache)

        self.logger.info(
            f"found {retrieved_candidates}, expected {expected_candidates}"
        )
        self.assertEqual(expected_candidates, retrieved_candidates)

        hist_and_new_values = {
            "ZTF18acvhwtf": {"ned_dist_hist": 0.52, "milliquas_dist_hist": 0.55},
            "ZTF20abgvabi": {},
        }

        for name, res in sorted(nu.cache.items()):
            latest = [
                x
                for x in res["prv_candidates"] + [res["candidate"]]
                if "isdiffpos" in x.keys()
            ][-1]

            if name == "ZTF18acvhwtf":
                old_flag = ""
                jds = [x["jd"] for x in res["prv_candidates"]]
                second_det = [x for x in jds if x > min(jds) + 0.01]
                if len(second_det) > 0:
                    if Time.now().jd - second_det[0] > 1.0:
                        old_flag = "(MORE THAN ONE DAY SINCE SECOND DETECTION)"

            hist_and_new_values[name]["ra"] = latest["ra"]
            hist_and_new_values[name]["dec"] = latest["dec"]
            hist_and_new_values[name]["mag"] = latest["magpsf"]
            hist_and_new_values[name]["mag_err"] = latest["sigmapsf"]

            ned_z, ned_dist = query_ned_for_z(
                ra_deg=latest["ra"],
                dec_deg=latest["dec"],
                searchradius_arcsec=20,
                logger=self.logger,
            )

            if ned_z:
                delta_to_historic_value = np.abs(
                    hist_and_new_values[name]["ned_dist_hist"] - ned_dist
                )

                self.assertTrue(delta_to_historic_value < self.max_distance_diff_arcsec)

                absmag = nu.calculate_abs_mag(latest["magpsf"], ned_z)
                z_dist = Distance(z=ned_z, cosmology=cosmo).value

                hist_and_new_values[name]["ned_dist_new"] = ned_dist
                hist_and_new_values[name]["absmag"] = absmag
                hist_and_new_values[name]["z_dist"] = z_dist

            milliquas_res = ampel_api_catalog(
                catalog="milliquas",
                catalog_type="extcats",
                ra_deg=res["candidate"]["ra"],
                dec_deg=res["candidate"]["dec"],
                search_radius_arcsec=1.5,
                logger=self.logger,
            )

            if milliquas_res:
                milliquas_dist = milliquas_res[0]["dist_arcsec"]
                delta_to_historic_value = np.abs(
                    hist_and_new_values[name]["milliquas_dist_hist"] - milliquas_dist
                )

                self.assertTrue(delta_to_historic_value < self.max_distance_diff_arcsec)

                hist_and_new_values[name]["milliquas_dist_new"] = milliquas_dist

        nu.plot_overlap_with_observations(
            first_det_window_days=(t_max - nu.t_min).to("d").value
        )
        res = nu.draft_gcn()

        print(repr(res))

        # Update the true using repr(res)
        true_gcn = f"Astronomer Name (Institute of Somewhere), ............. report,\n\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n\nAs part of the ZTF neutrino follow up program (Stein et al. 2022), we observed the localization region of the neutrino event IceCube-200620A (Santander et. al, GCN 27997) with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g- and r-band beginning at 2020-06-21 04:53 UTC, approximately 25.8 hours after event time. We covered 77.7% (1.2 sq deg) of the reported localization region. This estimate accounts for chip gaps. Each exposure was 300s with a typical depth of 21.0 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) . We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF18acvhwtf | AT2020ncs | {hist_and_new_values['ZTF18acvhwtf']['ra']:011.7f} | {hist_and_new_values['ZTF18acvhwtf']['dec']:+011.7f} | r      | {hist_and_new_values['ZTF18acvhwtf']['mag']:.2f} | {hist_and_new_values['ZTF18acvhwtf']['mag_err']:.2f}   | {old_flag} \n| ZTF20abgvabi | AT2020ncr | 162.5306341 | +12.1461187 | g      | 20.58 | 0.19   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \n\nZTF18acvhwtf was first detected on 2018-12-09. It has a spec-z of 0.291 [{hist_and_new_values['ZTF18acvhwtf']['z_dist']:.0f} Mpc] and an abs. mag of {hist_and_new_values['ZTF18acvhwtf']['absmag']:.1f}. Distance to SDSS galaxy is {hist_and_new_values['ZTF18acvhwtf']['ned_dist_new']:.2f} arcsec. [MILLIQUAS: SDSS J104816.25+120734.7 - 'Q'-type source ({hist_and_new_values['ZTF18acvhwtf']['milliquas_dist_new']:.2f} arsec)] [TNS NAME=AT2020ncs]\nZTF20abgvabi was first detected on 2020-05-26. WISEA J105007.28+120846.1 ['G'-type source (0.00 arsec)] [TNS NAME=AT2020ncr]\n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the nuztf (Stein et al. 2021, https://github.com/desy-multimessenger/nuztf).\n"

        self.assertEqual(res, true_gcn)

        # Test manually adding candidates

        nu.add_to_cache_by_names(
            "ZTF18abteipt",
        )

        # Check

        false_candidate = nu.check_ampel_filter("ZTF18abteipt")

        self.logger.info(
            f"For the false candidate, the pipeline bool is {false_candidate}"
        )

        self.assertFalse(false_candidate)

        true_candidate = nu.check_ampel_filter("ZTF20abgvabi")

        self.logger.info(
            f"For the true candidate, the pipeline bool is {true_candidate}"
        )

        self.assertTrue(true_candidate)
