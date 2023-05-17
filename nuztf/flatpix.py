import logging
import os
import pickle
from pathlib import Path

from gwemopt.ztf_tiling import get_quadrant_ipix
from tqdm import tqdm
from ztfquery.fields import FIELD_DATAFRAME

from nuztf.paths import FLATPIX_CACHE_DIR


def get_flatpix_path(nside: int) -> Path:
    """
    Flatpix lookup table for nside

    :param nside: nside of healpix
    :return: path to flatpix lookup table
    """
    return FLATPIX_CACHE_DIR.joinpath(f"ztf_fields_ipix_nside={nside}.pickle")


def get_nested_pix_path(nside: int) -> Path:
    """
    Nested pix lookup table for nside

    :param nside: nside of healpix
    :return: path to nested pix lookup table
    """

    outfile = FLATPIX_CACHE_DIR.joinpath(f"ztf_fields_nested_ipix_nside={nside}.pickle")
    return outfile


def generate_flatpix_file(nside: int, logger=logging.getLogger(__name__)):
    """
    Generate and save the fields-healpix lookup table

    :param nside: nside of healpix
    :param logger: logger
    :return: None
    """

    logger.info(f"Generating field-healpix lookup table for nside={nside}")

    field_dataframe = FIELD_DATAFRAME.reset_index()

    fields = field_dataframe["ID"].values
    ras = field_dataframe["RA"].values
    decs = field_dataframe["Dec"].values

    flat_pix_dict = dict()
    nested_pix_dict = dict()

    for i, field in tqdm(enumerate(fields), total=len(fields)):
        ra = ras[i]
        dec = decs[i]
        pix = get_quadrant_ipix(nside, ra, dec)

        flat_pix = []

        for sub_list in pix:
            for p in sub_list:
                flat_pix.append(p)

        flat_pix = list(set(flat_pix))
        flat_pix_dict[field] = flat_pix
        nested_pix_dict[field] = pix

    outfile = get_flatpix_path(nside=nside)

    logger.info(f"Saving to {outfile}")

    with open(outfile, "wb") as f:
        pickle.dump(flat_pix_dict, f)

    outfile = get_nested_pix_path(nside=nside)

    logger.info(f"Saving to {outfile}")

    with open(outfile, "wb") as f:
        pickle.dump(nested_pix_dict, f)


def get_flatpix(nside: int, logger=logging.getLogger(__name__)):
    """
    Get the fields-healpix lookup table

    :param nside: nside of healpix
    :param logger: logger
    :return: flatpix lookup table
    """
    infile = get_flatpix_path(nside=nside)

    # Generate a lookup table for field healpix
    # if none exists (because this is computationally costly)
    if not infile.exists():
        generate_flatpix_file(nside=nside, logger=logger)

    logger.info(f"Loading from {infile}")

    with open(infile, "rb") as f:
        field_pix = pickle.load(f)

    return field_pix


def get_nested_pix(nside: int, logger=logging.getLogger(__name__)):
    """
    Nested Flatpix lookup table for nside

    :param nside: nside of healpix
    :param logger: logger
    :return: flatpix lookup table
    """
    infile = get_nested_pix_path(nside=nside)

    # Generate a lookup table for field healpix
    # if none exists (because this is computationally costly)
    if not os.path.isfile(infile):
        generate_flatpix_file(nside=nside, logger=logger)

    logger.info(f"Loading from {infile}")

    with open(infile, "rb") as f:
        field_pix = pickle.load(f)

    return field_pix
