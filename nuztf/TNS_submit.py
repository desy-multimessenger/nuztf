#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, argparse, logging, re

# from ampel.ztf.archive.ArchiveDB import ArchiveDB
import ampel.contrib.hu.t3.TNSTalker as TNSTalker

# from ampel.ztf.utils.ZIAlertUtils import ZIAlertUtils
from ampel.ztf.dev.ZTFAlert import ZTFAlert
from ampel.ztf import legacy_utils

# from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils
from ampel.contrib.hu.t3.ampel_tns import sendTNSreports

from nuztf.credentials import load_credentials


def TNSSubmit(ztf_ids: list, reporter: str = None, sandbox: bool = True, logger=None):
    """ """
    if logger is None:
        logger = logging.getLogger(__name__)

    _, tns_api_token = load_credentials("tns_api_token")

    # Run config for TNSTalker
    ztf_default_values = {
        "flux_units": "1",
        "instrument_value": "196",
        "Observer": "Robot",
    }

    RUN_CONFIG = {
        "submit_tns": True,
        "sandbox": True,
        "resubmit_tns_nonztf": False,
        "resubmit_tns_ztf": False,
        "submit_unless_journal": True,
        "lc_filters": [],
        "ztf_tns_at": ztf_default_values,
        "max_age": 5,
        "tns_api_key": tns_api_token,
    }

    # Connect to TNSTalker instance
    talker = TNSTalker.TNSTalker(logger=logger, run_config=RUN_CONFIG)
    talker.run_config = talker.RunConfig(**RUN_CONFIG)


# # Connect to TNSTalker instance
# TNS = TNSTalker.TNSTalker(logger=logger, run_config=RUN_CONFIG)
# TNS.run_config = TNS.RunConfig(**RUN_CONFIG)
# PORT = 5432

# # Connect to AMPEL instance
# AMPEL_CLIENT = ArchiveDB(
#     "postgresql://{0}:{1}@localhost:{2}/ztfarchive".format(username, password, PORT)
# )


# def is_ztf_name(name):
#     """ """
#     return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


# def use_if_ztf(file_or_name):
#     """ """
#     if is_ztf_name(file_or_name):
#         object_list = [file_or_name]
#     else:
#         object_list = []
#         try:
#             file = open(f"{file_or_name}", "r")
#             lines = file.read().splitlines()
#             for line in lines:
#                 if is_ztf_name(line):
#                     object_list.append(line)
#         except FileNotFoundError as error:
#             print(
#                 "\nYou have to provide either a ZTF name or a file containing ZTF names (1 per line).\n"
#             )
#             raise error
#         assert (
#             object_list[0][:3] == "ZTF" and len(object_list[0]) == 12
#         ), "You have to provide either a ZTF name or a file containing ZTF names (1 per line)"

#     return object_list


# def TNS_submit(ztf_names, reporter=None, sandbox=True):
#     """
#     Note: IF YOU WANT TO SUBMIT, SET sandbox=False
#     """
#     atreports = {}

#     for ztf_name in ztf_names:
#         ztf_object = AMPEL_CLIENT.get_alerts_for_object(ztf_name, with_history=True)
#         query_res = [i for i in ztf_object]
#         last_alert = query_res[-1]
#         tv = ZIAlertUtils.to_transientview(
#             file_path=None, content=last_alert, science_records=None
#         )
#         ztf_id = tv.tran_id
#         ampel_id = ZTFUtils.to_ampel_id(ztf_id)
#         lc_id = tv.lightcurves[0].id
#         object.__setattr__(tv, "tran_id", ampel_id)
#         object.__setattr__(tv, "latest_state", lc_id)
#         atreport = TNS.create_atreport(tv)
#         if reporter is None:
#             atreport["reporter"] = "The Zwicky Transient Facility (ZTF) Collaboration"
#         else:
#             atreport[
#                 "reporter"
#             ] = f"{reporter} for the Zwicky Transient Facility (ZTF) Collaboration"
#         atreport["exptime"] = "300"
#         atreports.update({ztf_name: atreport})

#     # Reformat the report (code by AMPEL!)
#     atreportlist = []
#     atreport = {}
#     k = 0
#     for tranid in atreports.keys():
#         if len(atreport) > 90:
#             self.logger.info("adding another report to TNS submit")
#             atreportlist.append({"at_report": atreport})
#             atreport = {}
#             k = 0
#         atreport[int(k)] = atreports[tranid]
#         k += 1
#     atreportlist.append({"at_report": atreport})
#     tnsreplies = sendTNSreports(
#         atreportlist, TNS.run_config.tns_api_key, logger, sandbox=sandbox
#     )


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Submit ZTF object(s) to TNS")
#     parser.add_argument(
#         "name",
#         type=str,
#         help='Provide a ZTF name (e.g. "ZTF19aaelulu") or a .txt-file containing a list of ZTF names',
#     )
#     parser.add_argument(
#         "--commit",
#         "-commit",
#         action="store_false",
#         help="Run an actual commit (instead of just using the TNS sandbox",
#     )
#     parser.add_argument(
#         "--reporter",
#         "-reporter",
#         type=str,
#         default=None,
#         help="Provide an author for the submission",
#     )

#     commandline_args = parser.parse_args()
#     name = commandline_args.name
#     sandbox = commandline_args.commit
#     reporter = commandline_args.reporter
#     object_list = use_if_ztf(name)

#     if sandbox:
#         logger.info(
#             f"Submitting {object_list} to the TNS - SANDBOX MODE ONLY (to change, use -commit)"
#         )
#     else:
#         logger.info(f"Submitting {object_list} to the TNS")

#     TNS_submit(object_list, reporter, sandbox=sandbox)
