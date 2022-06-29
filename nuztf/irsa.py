#!/usr/bin/env python3
# coding: utf-8

import os
import logging

import numpy as np
import pandas as pd

from ztfquery.lightcurve import LCQuery
from astropy.table import Table
from ztfquery.io import LOCALSOURCE

from nuztf.style import plot_dir
from nuztf.plot import lightcurve_from_science_image
from nuztf.cat_match import context_from_name

logger = logging.getLogger(__name__)


def load_irsa(ra_deg: float, dec_deg: float, radius_arcsec: float = 0.5, **kwargs):
    df = LCQuery.from_position(ra_deg, dec_deg, radius_arcsec, **kwargs).data

    mask = df.catflags > 0

    flags = list(set(df.catflags))

    logger.info(
        f"Found {len(df)} datapoints, masking {np.sum(mask)} datapoints with bad flags."
    )

    for flag in sorted(flags):
        logger.debug(f"{np.sum(df.catflags == flag)} datapoints with flag {flag}")

    df = df.drop(df[mask].index)
    return df


def get_irsa_path(
    source_name: str,
    cache_dir: str = os.path.join(LOCALSOURCE, "cache/"),
) -> str:
    return os.path.join(cache_dir, f'{source_name.replace(" ", "")}.csv')


def plot_irsa_lightcurve(
    source_name: str,
    nu_name: list = None,
    source_coords: list = None,
    source_redshift: float = None,
    plot_mag: bool = False,
    atel: bool = True,
    plot_folder: str = plot_dir,
    extra_folder: str = None,
    check_obs: bool = False,
    check_obs_lookback_weeks: float = 4,
    from_cache: bool = False,
    cache_dir: str = os.path.join(LOCALSOURCE, "cache/"),
    expanded_labels: bool = True,
    ylim: tuple = None,
    radius_arcsec: float = 0.5,
):

    # Query IRSA, or load from cache

    try:
        os.makedirs(cache_dir)
    except OSError:
        pass

    cache_path = get_irsa_path(source_name, cache_dir=cache_dir)

    if from_cache:
        logger.debug(f"Loading from {cache_path}")
        df = pd.read_csv(cache_path)

    else:

        source_coords, source_redshift, _ = context_from_name(
            source_name=source_name,
            source_coords=source_coords,
            source_redshift=source_redshift,
        )

        df = load_irsa(source_coords[0], source_coords[1], radius_arcsec=radius_arcsec)

        logger.debug(f"Saving to {cache_path}")
        df.to_csv(cache_path)

    data = Table.from_pandas(df)

    logger.info(f"There are a total of {len(data)} detections for {source_name}")

    lightcurve_from_science_image(
        source_name=source_name,
        data=data,
        nu_name=nu_name,
        source_coords=source_coords,
        source_redshift=source_redshift,
        plot_mag=plot_mag,
        atel=atel,
        plot_folder=plot_folder,
        extra_folder=extra_folder,
        check_obs=check_obs,
        check_obs_lookback_weeks=check_obs_lookback_weeks,
        expanded_labels=expanded_labels,
        ylim=ylim,
    )
