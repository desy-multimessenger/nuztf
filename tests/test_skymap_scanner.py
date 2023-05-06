import logging
import unittest

from nuztf.skymap_scanner import SkymapScanner


class TestSkymapScanner(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_gw_scan(self):
        self.logger.info("\n\n Testing GW Scanning \n\n")

        gw_name = "S190814bv"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GW {gw_name}")

        scanner = SkymapScanner(
            event=gw_name,
            prob_threshold=prob_threshold,
            n_days=1,
        )

        scanner.plot_skymap()

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 1649

        self.logger.info(
            f"Retrieved {n_retrieved_alerts} alerts. {n_expected_alerts} alerts expected."
        )

        self.assertEqual(n_retrieved_alerts, n_expected_alerts)

        scanner.filter_alerts()

        n_retrieved_candidates = len(scanner.final_candidates)
        n_expected_candidates = 1

        self.logger.info(
            f"Retrieved {n_retrieved_candidates} candidates. {n_expected_candidates} candidates expected."
        )

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)

        self.logger.info("Creating overview table")

        scanner.create_overview_table()

        self.logger.info("Creating candidate summary")

        scanner.create_candidate_summary()

        fig, coverage_summary = scanner.plot_coverage()

        true_coverage_summary = "In total, 88.57 % of the contour was observed at least once.\nThis estimate includes 0.00 % of the contour at a galactic latitude <10 deg.\nIn total, 73.81 % of the contour was observed at least twice. \nIn total, 73.81 % of the contour was observed at least twice, and excluding low galactic latitudes.\nThese estimates account for chip gaps."

        self.assertEqual(coverage_summary, true_coverage_summary)

        tns_summary = scanner.tns_summary()

        true_tns_summary = "Candidate: ZTF19abpuhbh / RA=11.2569224 / Dec=-22.5161471 / First detection=2458710.9475231\nLast Upper Limit: None\nFirst Detection: 2458710.9475231 / band=r / mag=20.857 +/- 0.284\nFirst observed 13.56 hours after merger\n[2458710.9475231, 2458710.9948611]\n"

        self.assertEqual(tns_summary, true_tns_summary)

    def test_grb_scan(self):
        self.logger.info("\n\n Testing GRB Scanner \n\n")

        grb_name = "GRB210927A"
        prob_threshold = 0.9
        t_window_d = 0.1

        self.logger.info(f"Scanning with GRB {grb_name}")

        scanner = SkymapScanner(
            event=grb_name,
            prob_threshold=prob_threshold,
            n_days=t_window_d,
        )

        scanner.get_alerts()

        n_retrieved_alerts = scanner.n_alerts
        n_expected_alerts = 6356

        self.logger.info(
            f"Retrieved {n_retrieved_alerts} alerts. {n_expected_alerts} alerts expected."
        )

        self.assertEqual(n_retrieved_alerts, n_expected_alerts)

        scanner.filter_alerts()

        n_retrieved_candidates = len(scanner.final_candidates)
        n_expected_candidates = 3

        self.logger.info(
            f"Retrieved {n_retrieved_candidates} candidates. {n_expected_candidates} candidates expected."
        )

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)

        scanner.plot_coverage()
        res = scanner.draft_gcn()

        print(repr(res))

        # Update the true using repr(res)
        true_gcn = "Astronomer Name (Institute of Somewhere), ............. report,\n\nOn behalf of the Zwicky Transient Facility (ZTF) and Global Relay of Observatories Watching Transients Happen (GROWTH) collaborations: \n\nAs part of the ZTF neutrino follow up program (Stein et al. 2022), we observed the localization region of the GRB210927A with the Palomar 48-inch telescope, equipped with the 47 square degree ZTF camera (Bellm et al. 2019, Graham et al. 2019). We started observations in the g- and r-band beginning at 2021-09-27 02:44 UTC, approximately 2.1 hours after event time. We covered 0.1% (1.0 sq deg) of the reported localization region. This estimate accounts for chip gaps. Each exposure was 30s with a typical depth of 20.5 mag. \n \nThe images were processed in real-time through the ZTF reduction and image subtraction pipelines at IPAC to search for potential counterparts (Masci et al. 2019). AMPEL (Nordin et al. 2019, Stein et al. 2021) was used to search the alerts database for candidates. We reject stellar sources (Tachibana and Miller 2018) and moving objects, and apply machine learning algorithms (Mahabal et al. 2019) , and removing candidates with history of variability prior to the merger time. We are left with the following high-significance transient candidates by our pipeline, all lying within the 90.0% localization of the skymap.\n\n+--------------------------------------------------------------------------------+\n| ZTF Name     | IAU Name  | RA (deg)    | DEC (deg)   | Filter | Mag   | MagErr |\n+--------------------------------------------------------------------------------+\n| ZTF21acdvtxc |  -------  | 250.2336698 | +05.3908972 | g      | 21.80 | 0.21   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF21acdvtxp |  -------  | 250.4636648 | +01.8436867 | g      | 21.33 | 0.18   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n| ZTF21acdvuzf |  -------  | 241.8979602 | +19.0755373 | g      | 20.82 | 0.17   | (MORE THAN ONE DAY SINCE SECOND DETECTION) \n+--------------------------------------------------------------------------------+\n\n \n\nAmongst our candidates, \n\nZTF21acdvtxc had upper limit problems. PLEASE FILL IN NUMBERS BY HAND!!! WISEA J164056.10+052327.1 ['UvS'-type source (0.03 arsec)]\nZTF21acdvtxp had upper limit problems. PLEASE FILL IN NUMBERS BY HAND!!! [MILLIQUAS: SDSS J164151.27+015037.0 - Likely QSO (prob = 95.0%) (0.21 arsec)]\nZTF21acdvuzf had upper limit problems. PLEASE FILL IN NUMBERS BY HAND!!! \n\n\nZTF and GROWTH are worldwide collaborations comprising Caltech, USA; IPAC, USA; WIS, Israel; OKC, Sweden; JSI/UMd, USA; DESY, Germany; TANGO, Taiwan; UW Milwaukee, USA; LANL, USA; TCD, Ireland; IN2P3, France.\n\nGROWTH acknowledges generous support of the NSF under PIRE Grant No 1545949.\nAlert distribution service provided by DIRAC@UW (Patterson et al. 2019).\nAlert database searches are done by AMPEL (Nordin et al. 2019).\nAlert filtering is performed with the nuztf (Stein et al. 2021, https://github.com/desy-multimessenger/nuztf).\n"

        self.assertEqual(res, true_gcn)

    def test_gw_scan_desy_download(self):
        self.logger.info("\n\n Testing GW result download from DESY \n\n")

        gw_name = "S200115j"
        prob_threshold = 0.9

        self.logger.info(f"Scanning with GW {gw_name}")

        scanner = SkymapScanner(
            event=gw_name,
            prob_threshold=prob_threshold,
            n_days=2,
        )

        scanner.download_results()

        n_expected_candidates = 122

        n_retrieved_candidates = len(scanner.final_candidates)

        self.assertEqual(n_retrieved_candidates, n_expected_candidates)
