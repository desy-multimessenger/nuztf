#!/usr/bin/env python3
# coding: utf-8

import os
import matplotlib.pyplot as plt
import seaborn as sns

from ztfquery.io import LOCALSOURCE

sns.set_style("white")
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{romanbar}")
plt.rcParams["font.family"] = "sans-serif"

output_folder = os.path.join(LOCALSOURCE, "IRSA")
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

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
