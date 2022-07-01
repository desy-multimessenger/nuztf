#!/usr/bin/env python
# coding: utf-8

import os
import warnings
import logging
from ztfquery import io

# Manage ztfquery logins from environment variables


def load_credentials(name: str, token_based: bool = False):
    """ZTFquery wrapper for loading credentials."""
    return io._load_id_(name, token_based=token_based)


try:
    io.set_account(
        "irsa", username=os.environ["IRSA_USER"], password=os.environ["IRSA_PASSWORD"]
    )
    logging.info('Set up "irsa" credentials')

except KeyError:
    logging.info(
        'No Credentials for "irsa" found in environment' "Assuming they are set."
    )

try:
    io.set_account(
        "skyvision",
        username=os.environ["SKYVISION_USER"],
        password=os.environ["SKYVISION_PASSWORD"],
    )
    logging.info('Set up "skyvision" credentials')

except KeyError:
    logging.info(
        'No Credentials for "skyvision" found in environment' "Assuming they are set."
    )

try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        io.set_account(
            "ampel_api_archive_token",
            username=os.environ["AMPEL_API_ARCHIVE_TOKEN_USER"],
            password=os.environ["AMPEL_API_ARCHIVE_TOKEN_PASSWORD"],
        )
        logging.info('Set up "ampel_api_archive_token" credentials')

except KeyError:
    logging.info("No Token for AMPEL API found in environment" "Assume they are set.")

try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        io.set_account(
            "tns_api_token", token=os.environ["TNS_API_TOKEN"], token_based=True
        )
        logging.info('Set up "tns_api_token" credentials')

except KeyError:
    logging.info("No Token for TNS API found in environment" "Assume they are set.")

try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        io.set_account("fritz", token=os.environ["FRITZ_TOKEN"], token_based=True)
        logging.info('Set up "fritz" credentials')

except KeyError:
    logging.info("No Token for Fritz API found in environment" "Assume they are set.")
