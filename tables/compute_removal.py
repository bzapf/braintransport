import os, pathlib, json
import numpy as np
import pandas

import sys

from pathlib import Path
here = Path(__file__).parent
import sys
sys.path.insert(0, str(here.parent))
import h5py

from plots.definitions import groups, intervals
from plots.helpers import get_data_in_intervals, sim_at_mri_times
import argparse


pats = groups["all"]

# pats = [pats[0]]

if "cluster" in os.getcwd():

    dpath = pathlib.Path("/cluster/projects/nn9279k/bastian/SleepData/")

os.chdir("diffusion_reaction_simulations/")

pats = ["215"]

for pat in pats:

    path = dpath / pat / "diffusion_reaction" / "avgDTIavgT1"

    # path = dpath / pat / "TOTALCLEARANCE"

    subfolders = sorted(os.listdir(path), key=lambda x: float(x[1:]), reverse=True)

    subfolder = path / subfolders[0]

    data = json.load(open(subfolder / "params"))

    # rm = pandas.read_csv(subfolder / "rm.csv")

    # print(rm)

    # rm = pandas.read_csv(subfolder / "absrm.csv")

    # print(rm)


    # exit()

    # f = h5py.File(subfolder / "finalstate.hdf", "r")
    
    os.system("sbatch reaction.slurm " + pat + " 288 nodti avgt1 " + str(data["alpha_final"]) + " " + str(data["r_d_final"]))

    # os.system("sbatch reaction.slurm " + pat + " 144 nodti avgt1 1 0")