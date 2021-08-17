import pickle
import argparse
from multiprocessing import JoinableQueue, Process
from nuztf.gw_scanner import GravWaveScanner
import numpy as np
import healpy as hp
import os
from tqdm import tqdm
from astropy.time import Time
from pathlib import Path

ligo_candidate_cache = os.path.join(Path(__file__).resolve().parents[1], "LIGO_cache")

class MultiGwProcessor(GravWaveScanner):
    queue = None
    results = dict()

    def __init__(self, n_cpu=os.cpu_count()-1, mp_id=0, n_days=None, *args, **kwargs):

        print("Running on {0} CPUs".format(n_cpu))

        GravWaveScanner.__init__(self, n_days=n_days, *args, **kwargs)
        self.fill_queue, self.n_sky, self.scan_method = self.optimise_scan_method()

        self.cache_dir = os.path.join(
            ligo_candidate_cache,
            os.path.splitext(os.path.basename(self.output_path))[0]
        )
        try:
            os.makedirs(self.cache_dir)
        except OSError:
            pass
        self.scan_radius = np.degrees(hp.max_pixrad(self.cone_nside))
        self.queue = JoinableQueue()
        self.processes = [Process(target=self.scan_wrapper, kwargs={"mp_id": i+1}) for i in range(int(n_cpu))]

        self.obj_names = []
        self.mp_id = mp_id

        for p in self.processes:
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

            res = self.filter(query_res)

            # print("{0} of {1} queries: {2} accepted out of {3} alerts".format(j, mts, len(res), len(query_res)))

            self.obj_names += [x["objectId"] for x in res]
            if len(res) > 0:
                self.dump_cache()

            # self.dump_cache(res)
            self.queue.task_done()

    def filter(self, query_res):

        indexes = []

        for i, res in enumerate(query_res):
            if self.filter_f_no_prv(res):
                indexes.append(i)

        return [query_res[i] for i in indexes]

    def filter_f_no_prv(self, res):

        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            logging.debug("Transient is too old")
            return False

        # Veto new transients
        if res["candidate"]["jdstarthist"] > self.default_t_max.jd:
            logging.debug("Transient is too new")
            return False

        # Positive detection
        if res['candidate']['isdiffpos'] not in ["t", "1"]:
            logging.debug("Negative subtraction")
            return False

        try:
            if res['candidate']['drb'] < 0.3:
                logging.debug("DRB too low")
                return False
        except KeyError:
            pass
        except TypeError:
            pass

        # Check contour
        if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
            logging.debug("Not in contour")
            return False

        # Require 2 detections separated by 15 mins
        if (res["candidate"]["jdendhist"] - res["candidate"]["jdstarthist"]) < 0.01:
            logging.debug("Does not have 2 detections separated  by >15 mins")
            return False

        return True

    def dump_cache(self):

        path = os.path.join(self.cache_dir, "{0}.pkl".format(self.mp_id))

        with open(path, "wb") as f:
            pickle.dump(self.obj_names, f)

    def optimise_scan_method(self, t_max=None):

        if t_max is None:
            t_max = Time.now()

        length_window = t_max.jd - self.t_min.jd

        survey_start_jd = 2458000.

        full_ztf = t_max.jd - survey_start_jd

        n_skys_time = length_window * 1./4.

        n_skys_space = full_ztf * self.pixel_area/1000.

        if n_skys_time < n_skys_space:
            f = self.fill_queue_time
            method = "time"
        else:
            f = self.fill_queue_space
            method  = "space"

        n_sky = min([n_skys_space, n_skys_time])

        print("Scanning method: {0} \n N_sky: {1}".format(method, n_sky))

        return f, n_sky, method

    def fill_queue_time(self, t_max=None):

        if t_max is None:
            t_max = Time(self.default_t_max.jd + 7., format="jd")

        time_steps = np.arange(self.t_min.jd, t_max.jd, step=0.005)
        mts = len(time_steps)

        n_tot = 0

        print("Scanning between {0}JD and {1}JD".format(time_steps[0], time_steps[-1]))

        for j, t_start in enumerate(tqdm(list(time_steps[:-1]))):

            ztf_object = ampel_client.get_alerts_in_time_range(
                jd_min=t_start, jd_max=time_steps[j+1], with_history=False)
            query_res = [x for x in ztf_object]
            n_tot += len(query_res)
            self.add_to_queue((j, mts, query_res))
            self.scanned_pixels.append(j)

        print("Added {0} candidates since {1}".format(n_tot, time_steps[0]))

    def fill_queue_space(self, t_max=None):

        if t_max is None:
            t_max = Time.now()

        mts = len(list(self.cone_ids))
        n_tot = 0

        for j, cone_id in enumerate(tqdm(list(self.cone_ids))):
            ra, dec = self.cone_coords[j]

            if cone_id not in self.scanned_pixels:
                ztf_object = ampel_client.get_alerts_in_cone(
                    ra, dec, self.scan_radius, self.t_min.jd, t_max.jd, with_history=False)
                query_res = [x for x in ztf_object]
                n_tot += len(query_res)
                self.add_to_queue((j, mts, query_res))
                self.scanned_pixels.append(j)

        print("Added {0} candidates since {1}".format(n_tot, self.t_min.jd))

    def terminate(self):
        """ wait until queue is empty and terminate processes """
        self.queue.join()
        for p in self.processes:
            p.terminate()

    def combine_cache(self):

        for name in self.get_cache_file():
            try:
                with open(os.path.join(self.cache_dir, name), "rb") as f:
                    self.obj_names += pickle.load(f)
            except:
                pass

        self.obj_names = list(set(self.obj_names))

        print("Scanned {0} pixels".format(len(self.scanned_pixels)))
        print("Found {0} candidates passing the first filtering stage.".format(len(self.obj_names)))
        print(self.obj_names)

        ztf_object = ampel_client.get_alerts_for_object(self.obj_names, with_history=True)

        query_res = [i for i in ztf_object]

        query_res = self.merge_alerts(query_res)

        self.cache = dict()

        for res in tqdm(query_res):
            if self.filter_f_history(res):
                if self.filter_ampel(res):
                    self.cache[res["objectId"]] = res
                else:
                    print("Failed Ampel")
            else:
                print("Failed History")

        print("Found {0} candidates passing the final filtering stage.".format(len(self.cache)))

        self.create_candidate_summary()

    def get_cache_file(self):
        return [os.path.join(self.cache_dir, x) for x in os.listdir(self.cache_dir) if ".pkl" in x]

    def clean_cache(self):

        for name in self.get_cache_file():
            os.remove(name)

        print("Cache cleaned!")

if __name__ == '__main__':
    import os
    import logging
    from ampel.pipeline.logging.ExtraLogFormatter import ExtraLogFormatter

    logging.basicConfig()
    logging.root.handlers[0].setFormatter(ExtraLogFormatter())
    logger = logging.getLogger()
    logger.setLevel(logging.ERROR)

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_cpu", default=min(24, max(1, os.cpu_count()-1)))
    parser.add_argument("-p", "--prob_threshold", default=0.9, type=float)
    parser.add_argument("-n", "--name", default=None)
    cfg = parser.parse_args()

    print("N CPU available", os.cpu_count())
    print("Using {0} CPUs".format(cfg.n_cpu))

    gw = MultiGwProcessor(gw_name=cfg.name, logger=logger, prob_threshold=cfg.prob_threshold,
                          n_cpu=cfg.n_cpu, fast_query=True)
    gw.clean_cache()
    gw.fill_queue()
    gw.terminate()
    gw.combine_cache()
    gw.clean_cache()
    # gw.plot_overlap_with_observations()