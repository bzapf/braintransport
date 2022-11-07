import os
import argparse
import matplotlib
from definitions import groups, resultfoldername
from barplots import make_barplot

parser = argparse.ArgumentParser()

default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/simulations/"

parser.add_argument("--inputpath", type=str, default=default_input_path)
parser.add_argument("--outputpath", type=str, default=default_output_path)

parserargs = parser.parse_args()
argparse_dict = vars(parserargs)


path_to_files = argparse_dict["inputpath"]
plotpath = argparse_dict["outputpath"]

if path_to_files != default_input_path:
    assert plotpath != default_output_path

os.chdir(path_to_files)

group = "all"
pats = groups[group]

alphas = [1, 2, 3, 4, 5]

imageformat = ".png"

FS = 35

dpi = 400
figsize = (10, 10)

matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["axes.labelsize"] = FS # "xx-large"
matplotlib.rcParams["grid.linewidth"] = 1
matplotlib.rcParams["xtick.labelsize"] = FS # "xx-large"
matplotlib.rcParams["ytick.labelsize"] = FS # "xx-large"
matplotlib.rcParams["legend.fontsize"] = FS # "xx-large"
matplotlib.rcParams["font.size"] = FS

fs = None

paperformat = True


for region in ["avg", "white", "gray"]:
    print(region)

    make_barplot(region, pats, alphas, paperformat, resultfoldername, path_to_files, 
                savepath=plotpath + "alpha_barpolots_" + region + imageformat, fs=fs, figsize=figsize, dpi=dpi)

