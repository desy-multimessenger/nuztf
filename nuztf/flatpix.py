import logging
import os
import pickle

from gwemopt.ztf_tiling import get_quadrant_ipix
from tqdm import tqdm
from ztfquery.fields import FIELD_DATAFRAME


def get_flatpix_path(nside: int):
    outdir = os.path.join(os.path.dirname(__file__), "data")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile = os.path.join(outdir, f"ztf_fields_ipix_nside={nside}.pickle")
    return outfile


def generate_flatpix_file(nside: int, logger=logging.getLogger(__name__)):
    """
    Generate and save the fields-healpix lookup table
    """

    logger.info(f"Generating field-healpix lookup table for nside={nside}")

    field_dataframe = FIELD_DATAFRAME.reset_index()

    fields = field_dataframe["ID"].values
    ras = field_dataframe["RA"].values
    decs = field_dataframe["Dec"].values

    flat_pix_dict = dict()

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

    outfile = get_flatpix_path(nside=nside)

    logger.info(f"Saving to {outfile}")

    with open(outfile, "wb") as f:
        pickle.dump(flat_pix_dict, f)


def get_flatpix(nside: int, logger=logging.getLogger(__name__)):
    infile = get_flatpix_path(nside=nside)

    # Generate a lookup table for field healpix
    # if none exists (because this is computationally costly)
    if not os.path.isfile(infile):
        generate_flatpix_file(nside=nside, logger=logger)

    logger.info(f"Loading from {infile}")

    with open(infile, "rb") as f:
        field_pix = pickle.load(f)

    return field_pix
