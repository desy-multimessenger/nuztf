import pickle
import argparse
from multiprocessing import JoinableQueue, Process
from gw_scanner import GravWaveScanner
import random
import numpy as np
import healpy as hp
import os
from pathlib import Path
from ampel_magic import ampel_client
from tqdm import tqdm
import pickle as pickle


ligo_candidate_cache = os.path.join(Path().absolute(), "LIGO_cache")

class MultiGwProcessor(GravWaveScanner):
    queue = None
    results = dict()

    def __init__(self, n_cpu, **kwargs):
        GravWaveScanner.__init__(self, **kwargs)
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
        self.processes = [Process(target=self.scan_wrapper) for _ in range(int(n_cpu))]

        for p in self.processes:
            p.start()

    def add_to_queue(self, item):
        self.queue.put(item)

    def scan_wrapper(self):
        while True:
            ztf_object = self.queue.get()
            if ztf_object is None:
                break

            query_res = [x for x in ztf_object]
            res = self.filter(query_res)

            self.dump_cache(res)
            self.queue.task_done()

    def filter(self, query_res):

        indexes = []

        for i, res in enumerate(query_res):
            if self.fast_filter_f_no_prv(res):
                if self.filter_ampel(res) is not None:
                    indexes.append(i)

        return [query_res[i] for i in indexes]

    def fast_filter_f_no_prv(self, res):

        # Positive detection
        if res['candidate']['isdiffpos'] not in ["t", "1"]:
            return False

        # Veto old transients
        if res["candidate"]["jdstarthist"] < self.t_min.jd:
            return False

        # # Check contour
        # if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
        #     return False

        return True

    def dump_cache(self, res):
        for obj in res:
            path = os.path.join(self.cache_dir, obj["objectId"]+".pkl")
            with open(path, "wb") as f:
                pickle.dump(obj, f)

    def fill_queue(self):

        t_max = self.default_t_max

        for i, cone_id in enumerate(tqdm(list(self.cone_ids))):
            self.scanned_pixels.append(cone_id)
            ra, dec = self.cone_coords[i]
            ztf_object = ampel_client.get_alerts_in_cone(
                ra, dec, self.scan_radius, self.t_min.jd, t_max.jd, with_history=False)
            query_res = [x for x in ztf_object]
            r.add_to_queue(query_res)

    def terminate(self):
        """ wait until queue is empty and terminate processes """
        self.queue.join()
        for p in self.processes:
            p.terminate()

    def combine_cache(self):
        self.cache = dict()

        for name in self.get_cache_file():
            try:
                with open(os.path.join(self.cache_dir, name), "rb") as f:
                    obj = pickle.load(f)
                    self.cache[obj["objectId"]] = obj
            except:
                pass

        print("Scanned {0} pixels".format(len(self.scanned_pixels)))
        print("Found {0} candidates".format(len(self.cache)))

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

    logger = logging.getLogger("quiet_logger")
    logger.setLevel(logging.ERROR)

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--n_cpu", default=max(1, os.cpu_count()-1))
    cfg = parser.parse_args()

    print("N CPU available", os.cpu_count())

    r = MultiGwProcessor(n_cpu=cfg.n_cpu, logger=logger, fast_query=True)
    r.clean_cache()
    r.fill_queue()
    r.terminate()
    r.combine_cache()
    r.clean_cache()