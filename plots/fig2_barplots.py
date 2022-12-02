import os
import argparse
if "cluster" not in os.getcwd():
    import matplotlib
from definitions import groups, ylabels # , resultfoldername
from barplots import make_barplot

parser = argparse.ArgumentParser()

default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/simulations/"

parser.add_argument("--inputpath", type=str, default=default_input_path)
parser.add_argument("--outputpath", type=str, default=default_output_path)
parser.add_argument("--printAvg", action="store_true", default=False, help="print info and make plots using average tracer (mmol/L)")

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

FS = 36

dpi = 400
figsize = (12, 8)

if "cluster" not in os.getcwd():
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


GREY_WHITE=False

plotname = ""

if GREY_WHITE:
    plotname = "Grey-White"


print_format=".4f"

width = 0.8 + 0.1

def resultfoldername(pat, alpha):

    return str(pat) + "/alphatests/alpha" + str(alpha) + "/"

for region in ["avg", 
                "white", 
                "gray"]:

    ylabel = ylabels[region]

    print("Region=", region)

    print("Total tracer (mmol)", "." * 50)

    make_barplot(region, pats, alphas, paperformat, resultfoldername, path_to_files, GREY_WHITE=GREY_WHITE, width=width, ylabel=ylabel, print_format=print_format,
                savepath=plotpath + plotname + "alpha_barpolots_" + region + imageformat, fs=fs, figsize=figsize, dpi=dpi, average_tracer=False)

    print("Average tracer (mmol / L)", "." * 50)
    if not argparse_dict["printAvg"]:
        continue
    
    make_barplot(region, pats, alphas, paperformat, resultfoldername, path_to_files, GREY_WHITE=GREY_WHITE, width=width, ylabel=ylabel, print_format=print_format,
                savepath=plotpath + plotname + "alpha_barpolots__avg_" + region + imageformat, fs=fs, figsize=figsize, dpi=dpi, average_tracer=True)

    # # exit()

    print("---------------------------------------------------------------------------------------------------------------")
    print("---------------------------------------------------------------------------------------------------------------")