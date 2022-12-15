import os, pathlib, json
import numpy as np
import pandas

import sys
import matplotlib.pyplot as plt

from pathlib import Path
here = Path(__file__).parent
import sys
sys.path.insert(0, str(here.parent))

from plots.definitions import groups, intervals
from plots.helpers import get_data_in_intervals, sim_at_mri_times
import argparse
import scipy
import scipy.integrate

pats = groups["all"]
pats.remove("091")
# pats = [pats[0]]


dpath = pathlib.Path("/home/basti/Dropbox (UiO)/Sleep/")

roi = "avg"

for pat in pats:

    # path = dpath / pat / "diffusion_reaction" / "avgDTIavgT1"

    path = dpath / pat / "TOTALCLEARANCE"
    subfolders = sorted(os.listdir(path), key=lambda x: float(x[1:]), reverse=True)

    subfolder = path / subfolders[0]

    data = json.load(open(subfolder / "params"))

    volumes = json.load(open(subfolder / "region_volumes.json"))

    removed_tracer = pandas.read_csv(subfolder / "rm.csv")

    # print(removed_tracer)
    # exit()

    boundary_tracer = pandas.read_csv(subfolder / "boundaryconcentration.csv")

    totalremoved = scipy.integrate.trapezoid(y=removed_tracer[roi], x=removed_tracer["t"])

    totalremoved *= (volumes[roi] / 1e6)

    

    plt.plot(removed_tracer["t"] / 3600, removed_tracer[roi] * volumes[roi] / 1e6)
    # plt.plot(boundary_tracer["t"] / 3600, boundary_tracer["outflux"] * volumes[roi] / 1e6)

    

    print(pat, totalremoved, data["reaction_rate"])

plt.show()