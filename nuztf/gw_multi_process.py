import pickle
import argparse
from multiprocessing import JoinableQueue, Process, freeze_support
from nuztf.gw_scanner import GravWaveScanner
import numpy as np
import healpy as hp
import os
from tqdm import tqdm
from astropy.time import Time
from pathlib import Path
from nuztf.ampel_api import ampel_api_cone, ampel_api_timerange, ampel_api_name
import logging

ligo_candidate_cache = os.path.join(Path(__file__).resolve().parents[1], "LIGO_cache")


class MultiGwProcessor(GravWaveScanner):
    queue = None
    results = dict()

    def __init__(
        self,
        n_cpu: int = os.cpu_count() - 1,
        mp_id: int = 0,
        n_days=None,
        verbose: bool = False,
        logger=None,
        *args,
        **kwargs,
    ):

        self.mp_id = mp_id
        self.n_cpu = n_cpu
        self.verbose = verbose

        GravWaveScanner.__init__(self, n_days=n_days, verbose=verbose, *args, **kwargs)

        print(self.t_min.jd)
        print(self.default_t_max)

        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)

        if verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        self.logger.info(f"Running on {self.n_cpu} CPUs")

        self.logger.info("Now the queue will be filled")

        self.fill_queue, self.n_sky, self.scan_method = self.optimize_scan_method()

        self.cache_dir = os.path.join(
            ligo_candidate_cache,
            os.path.splitext(os.path.basename(self.output_path))[0],
        )

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

        pass

    def run(self):

        self.scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))
        self.queue = JoinableQueue()

        kwargs = {}
        for i in range(int(self.n_cpu)):
            kwargs.update({"mp_id": i + 1})

        self.processes = [
            Process(target=self.scan_wrapper, kwargs={"mp_id": i + 1}) for i in [0]
        ]  # range(int(self.n_cpu))]

        self.obj_names = []

        for p in self.processes:
            print(p)
            p.start()

    def add_to_queue(self, item):
        self.queue.put(item)

    def scan_wrapper(self, **kwargs):

        self.mp_id = kwargs["mp_id"]

        while True:
            item = self.queue.get()
            if item is None:
                break

            (j, mts, query_res) = item

            # print("{0} of {1} queries: Staring {2} alerts".format(j, mts, len(query_res)))
            print("now checking")
            res = self.filter(query_res)

            # print("{0} of {1} queries: {2} accepted out of {3} alerts".format(j, mts, len(res), len(query_res)))

            self.obj_names += [x["objectId"] for x in res]

            if len(res) > 0:
                print("Dumping cache")
                self.dump_cache()

            # self.dump_cache(res)
            self.queue.task_done()

    def filter(self, query_res):

        indices = []

        for i, res in enumerate(query_res):
            if self.filter_f_no_prv(res):
                indices.append(i)

        return [query_res[i] for i in indices]

    def filter_f_no_prv(self, res):

        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            if self.verbose:
                print(f"{res['objectId']}: Transient is too old")
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > self.default_t_max.jd:
            if self.verbose:
                print(f"{res['objectId']}: Transient is too new")
            return False

        # Positive detection
        if res["candidate"]["isdiffpos"] not in ["t", "1"]:
            if self.verbose:
                print(f"{res['objectId']}: Negative subtraction")
            return False

        try:
            if res["candidate"]["drb"] < 0.3:
                if self.verbose:
                    print(f"{res['objectId']}: DRB too low")
                return False
        except KeyError:
            pass
        except TypeError:
            pass

        # Check contour
        print(f"ra={res['candidate']['ra']}")
        print(f"dec={res['candidate']['dec']}")
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            if self.verbose:
                print(f"{res['objectId']}: Not in contour")
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            if self.verbose:
                print(
                    f"{res['objectId']}: Does not have 2 detections separated  by >15 mins"
                )
            return False

        print(f"{res['objectId']}: Passed first filtering stage")
        self.logger.info(f"{res['objectId']}: Passed first filtering stage")

        return True

    def dump_cache(self):

        path = os.path.join(self.cache_dir, "{0}.pkl".format(self.mp_id))

        with open(path, "wb") as f:
            pickle.dump(self.obj_names, f)

    def optimize_scan_method(self, t_max=None):

        if t_max is None:
            t_max = Time.now()

        length_window = t_max.jd - self.t_min.jd

        survey_start_jd = 2458000.0

        full_ztf = t_max.jd - survey_start_jd

        n_skys_time = length_window * 1.0 / 4.0

        n_skys_space = full_ztf * self.pixel_area / 1000.0

        n_skys_time = 2323123123123

        if n_skys_time < n_skys_space:
            f = self.fill_queue_time
            method = "time"
        else:
            f = self.fill_queue_space
            method = "space"

        n_sky = min([n_skys_space, n_skys_time])

        self.logger.info(f"Scanning method: {method} \n N_sky: {n_sky}")

        return f, n_sky, method
        # return n_sky, method

    def fill_queue_time(self, t_max=None):
        """query the AMPEL API based on time-range search"""
        if t_max is None:
            t_max = Time(self.default_t_max.jd + 4.0, format="jd")

        # time_steps = np.arange(self.t_min.jd, t_max.jd, step=0.005)
        time_steps = np.arange(self.t_min.jd, t_max.jd, step=1)
        mts = len(time_steps)

        n_tot = 0

        self.logger.info(f"Scanning between {time_steps[0]}JD and {time_steps[-1]} JD")

        for j, t_start in enumerate(tqdm(list(time_steps[:-1]))):

            jd_min = Time(t_start, format="jd").jd
            jd_max = Time(time_steps[j + 1], format="jd").jd

            query_res = ampel_api_timerange(
                t_min_jd=jd_min,
                t_max_jd=jd_max,
                with_history=False,
                chunk_size=1000,
                logger=self.logger,
            )

            n_tot += len(query_res)
            self.add_to_queue((j, mts, query_res))
            self.scanned_pixels.append(j)

        self.logger.info(f"Added {n_tot} candidates since {time_steps[0]}")

    def fill_queue_space(self, t_max=None):
        """query the AMPEL API based on cone search"""
        if t_max is None:
            t_max = Time(self.default_t_max.jd, format="jd")

        mts = len(list(self.cone_ids))
        n_tot = 0

        for j, cone_id in enumerate(tqdm(list(self.cone_ids))):
            ra, dec = self.cone_coords[j]

            if cone_id not in self.scanned_pixels:

                query_res = ampel_api_cone(
                    ra=ra,
                    dec=dec,
                    radius=self.scan_radius,
                    t_min_jd=self.t_min.jd,
                    t_max_jd=t_max.jd,
                    logger=self.logger,
                )
                print(f"{len(query_res)} alerts found.")

                n_tot += len(query_res)
                self.add_to_queue((j, mts, query_res))
                self.scanned_pixels.append(j)

        self.logger.info(f"Added {n_tot} candidates since {self.t_min.jd}")

    def terminate(self):
        """wait until queue is empty and terminate processes"""
        self.queue.join()
        for p in self.processes:
            p.terminate()

    def combine_cache(self):
        """read the pickled result from first filtering stage and cut more"""

        for name in self.get_cache_file():
            print(name)
            try:
                with open(os.path.join(self.cache_dir, name), "rb") as f:
                    self.obj_names += pickle.load(f)
            except:
                pass

        self.obj_names = list(set(self.obj_names))

        self.logger.info(f"Scanned {len(self.scanned_pixels)} pixels")
        self.logger.info(
            f"Found {len(self.obj_names)} candidates passing the first filtering stage."
        )

        self.logger.info(f"Now checking: {self.obj_names}")

        all_results = []

        for ztf_name in self.obj_names:

            query_res = ampel_api_name(ztf_name=ztf_name, with_history=True)

            for res in tqdm(query_res):
                if self.filter_f_history(res):
                    if self.filter_ampel(res):
                        self.logger.info(f"{res['objectId']}: Passed all cuts")
                        self.cache[res["objectId"]] = res
                    else:
                        self.logger.info(f"{res['objectId']}: Failed Ampel")
                else:
                    self.logger.info(f"{res['objectId']}: Failed History")

        self.logger.info(
            f"Found {len(self.cache)} candidates passing the final filtering stage."
        )

        self.create_candidate_summary()

    def get_cache_file(self):
        return [
            os.path.join(self.cache_dir, x)
            for x in os.listdir(self.cache_dir)
            if ".pkl" in x
        ]

    def clean_cache(self):

        for name in self.get_cache_file():
            os.remove(name)

        self.logger.info("Cache cleaned!")


if __name__ == "__main__":

    import os
    import logging
    from ampel.pipeline.logging.ExtraLogFormatter import ExtraLogFormatter

    logging.basicConfig()
    logging.root.handlers[0].setFormatter(ExtraLogFormatter())
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_cpu", default=min(24, max(1, os.cpu_count() - 1)))
    parser.add_argument("-p", "--prob_threshold", default=0.9, type=float)
    parser.add_argument("-n", "--name", default=None)
    cfg = parser.parse_args()

    logger.info("N CPU available", os.cpu_count())
    logger.info("Using {0} CPUs".format(cfg.n_cpu))

    gw = MultiGwProcessor(
        gw_name=cfg.name,
        logger=logger,
        prob_threshold=cfg.prob_threshold,
        n_cpu=cfg.n_cpu,
        verbose=True,
    )

    gw.clean_cache()
    gw.fill_queue()
    gw.terminate()
    gw.combine_cache()
    gw.clean_cache()
    # gw.plot_overlap_with_observations()
