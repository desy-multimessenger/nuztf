#!/usr/bin/env python3
# coding: utf-8

import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from ztfquery.io import LOCALSOURCE


logger = logging.getLogger(__name__)

sns.set_style("white")

# Use latex if available
try:
    subprocess.check_output(["which", "latex"])
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{romanbar}")
except subprocess.CalledProcessError:
    logger.warning(
        "No Latex installation found. Proceeding without, but plots may look weird."
    )

plt.rcParams["font.family"] = "sans-serif"

plot_dir = os.path.join(LOCALSOURCE, "plots")
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

dpi = 300

fontsize = 7.0
big_fontsize = 10.0

golden_ratio = 1.618

base_width = 4.0
base_height = base_width / golden_ratio

margin_width = 0.5 * base_width
margin_height = margin_width / golden_ratio

full_width = 1.5 * base_width
full_height_landscape = full_width / golden_ratio
full_height_a4 = 11.75 / 8.25 * full_width

cmap = "rocket"
