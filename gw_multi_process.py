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

    def __init__(self, n_cpu, id=0, **kwargs):
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
        self.processes = [Process(target=self.scan_wrapper, kwargs={"id": i+1}) for i in range(int(n_cpu))]

        self.obj_names = []
        self.mp_id = id

        for p in self.processes:
            p.start()

    def add_to_queue(self, item):
        self.queue.put(item)

    def scan_wrapper(self, **kwargs):
        self.mp_id = kwargs["id"]

        while True:
            ztf_object = self.queue.get()
            if ztf_object is None:
                break

            query_res = [x for x in ztf_object]
            res = self.filter(query_res)

            self.obj_names += [x["objectId"] for x in res]

            # self.dump_cache(res)
            self.queue.task_done()

        self.dump_cache()

    def filter(self, query_res):

        indexes = []

        for i, res in enumerate(query_res):
            # print(i)
            if self.filter_f_no_prv(res):
                if self.filter_ampel(res) is not None:
                    indexes.append(i)

        return [query_res[i] for i in indexes]

    def filter_f_no_prv(self, res):

        # Positive detection
        if res['candidate']['isdiffpos'] not in ["t", "1"]:
            return False

        # Veto old transients
        # if res["candidate"]["jdstarthist"] < self.t_min.jd:
        #     return False

        # # Check contour
        # if not self.in_contour(res["candidate"]["ra"], res["candidate"]["dec"]):
        #     return False

        return True

    # def dump_cache(self, res):
    #     for obj in res:
    #         path = os.path.join(self.cache_dir, obj["objectId"] + ".pkl")
    #         with open(path, "wb") as f:
    #             pickle.dump(obj, f)
    def dump_cache(self):

        path = os.path.join(self.cache_dir, "{0}.pkl".format(self.mp_id))

        with open(path, "rb") as f:
            pickle.dump(self.obj_names, f)

    def fill_queue(self, max_cones=None):

        if max_cones is None:
            max_cones = len(self.cone_ids)

        t_max = self.default_t_max

        for k, cone_id in enumerate(tqdm(list(self.cone_ids)[:max_cones])):
            if cone_id not in self.scanned_pixels:
                ra, dec = self.cone_coords[k]
                ztf_object = ampel_client.get_alerts_in_cone(
                    np.degrees(ra), np.degrees(dec), self.scan_radius, self.t_min.jd-10, t_max.jd, with_history=False)
                query_res = [x for x in ztf_object]
                # print(query_res, ra, dec)
                r.add_to_queue(query_res)
                self.scanned_pixels.append(cone_id)

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
        print("Found {0} candidates".format(len(self.obj_names)))

        # self.create_candidate_summary()

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
    parser.add_argument("-n", "--n_cpu", default=min(24, max(1, os.cpu_count()-1)))
    parser.add_argument("-p", "--prob_threshold", default=0.9)
    cfg = parser.parse_args()

    print("N CPU available", os.cpu_count())

    r = MultiGwProcessor(n_cpu=cfg.n_cpu, logger=logger, prob_threshold=cfg.prob_threshold, fast_query=True)
    r.clean_cache()
    r.fill_queue(max_cones=100)
    r.terminate()
    r.combine_cache()
    # r.clean_cache()