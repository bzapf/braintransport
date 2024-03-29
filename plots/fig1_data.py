import pandas as pd
import numpy as np
import os
import pathlib, json
import matplotlib.pyplot as plt
import argparse
import scipy
from scipy import stats
from definitions import groups, labels, intervals
from helpers import get_data_in_intervals, significance_bar



parser = argparse.ArgumentParser()
parser.add_argument("--ylim", default=None, type=float)

default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/data/"

parser.add_argument("--inputpath", type=str, default=default_input_path)
parser.add_argument("--outputpath", type=str, default=default_output_path)


parserargs = parser.parse_args()
argparse_dict = vars(parserargs)


path_to_files = argparse_dict["inputpath"]
plotpath = argparse_dict["outputpath"]

if path_to_files != default_input_path:
    assert plotpath != default_output_path


params = {'font.size': 10}
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams.update(params) 

fac = 1
dpi = 300 * fac


only_significant = True

latex_textwidth = 7.02352416667 # inches
figwidth = latex_textwidth / 3
figheight = figwidth * 0.75

size_inches = (figwidth, figheight)

width = 0.4  # the width of the bars
capsize_ = 2
fontsize = None


def stderr(x, axis): 
    return np.nanstd(x, axis=axis)
    # return stats.sem(x, axis=axis, nan_policy="omit")



group = "all"
pats = groups[group]
sleepers = groups["sleep"]
nonsleep = groups["sleepdep"]


def resultfoldername(pat):
    return "alphatests"


os.makedirs(plotpath, exist_ok=True)

for roi, ylabel in zip(["avg", "gray", "white", "avgds"], ['tracer in brain (mmol)', 'tracer in gray (mmol)', 'tracer in white (mmol)', r'tracer at surface (mmol/m)']):


    # fig1 = plt.figure()
    # ax1 = fig1.gca()

    # fig2 = plt.figure()
    # ax2 = fig2.gca()

    # fig3 = plt.figure()
    # ax3 = fig3.gca()

    tracer_dict = {}
    avg_tracer_dict = {}

    for pat in pats:

        pat = str(pat)

        excelpaths = []

        subfolder = pathlib.Path(path_to_files) / str(pat) / resultfoldername(pat) / "alpha1"

        v = 1

        if roi == "avgds":
            filename = 'region_areas.json'

        else:                
            filename = 'region_volumes.json'

        with open(subfolder / filename) as f:
            region_volumes = json.load(f)
            roi_volume = region_volumes[roi] / 1e6

        # if roi != "avgds":
        #     with open(path_to_files + str(pat) + '/region_volumes' + str(32), 'rb') as f:
        #         region_volumes = pickle.load(f)
        #     for r in roi:
        #         if "ds" in roi:
        #             raise KeyError
        # else:
        #     print("Loading surface areas")
        #     assert roi == "avgds"
        #     with open(path_to_files + str(pat) + '/region_areas' + str(32), 'rb') as f:
        #         region_volumes = pickle.load(f)

        # v = region_volumes[roi] * 1e-6

        try:
            exceltable = pd.read_csv(subfolder / "experimental_data.csv")
        except pd.errors.EmptyDataError:
            breakpoint()
            print(subfolder + "experimental_data.csv")

        t = exceltable["t"]

        assert max(t) < 2.5 * 24 * 3600

        total_tracer = exceltable[roi] * v
        average_tracer = exceltable[roi]

        if pat in sleepers:
            c = "blue"

        else:
            c = "red"

        ts, tracer_at_times = get_data_in_intervals(pat, stored_times=t, stored_data=total_tracer, intervals=intervals)
        _, avg_tracer_at_times = get_data_in_intervals(pat, stored_times=t, stored_data=average_tracer, intervals=intervals)

        tracer_dict[pat] = tracer_at_times.tolist()

        if pat in sleepers:
            tracer_dict[pat] = tracer_dict[pat] + ["sleep"]
        else:
            tracer_dict[pat] = tracer_dict[pat] + ["no sleep"]

        avg_tracer_dict[pat] = avg_tracer_at_times.tolist()

        if pat in sleepers:
            avg_tracer_dict[pat] = avg_tracer_dict[pat] + ["sleep"]
        else:
            avg_tracer_dict[pat] = avg_tracer_dict[pat] + ["no sleep"]

    patdf = pd.DataFrame.from_dict(tracer_dict).transpose()

    avg_tracer_dict= pd.DataFrame.from_dict(avg_tracer_dict).transpose()

    if type(patdf.iloc[0, patdf.columns[-1]]) is str:
        patdf = patdf.loc[:, :(patdf.columns[-1]-1)]
    
    if type(avg_tracer_dict.iloc[0, avg_tracer_dict.columns[-1]]) is str:
        avg_tracer_dict = avg_tracer_dict.loc[:, :(avg_tracer_dict.columns[-1]-1)]

    print("------------------------------------------------------------------------------")
    print("Total amount of tracer in " + roi + ", averaged over all subjects:")
    for i in range(1, 4):
        y = patdf[i]

        avg_tracer_at_i = avg_tracer_dict[i]

        mean, std = np.nanmean(y), np.nanstd(y)

        print("time", intervals[i], "total tracer  ", format(mean, ".4f"),  "pm", format(std, ".4f"), "mmol", 
            "(", format(mean * 100 / 0.5, ".0f"),  "pm", format(std * 100 / 0.5, ".0f"), " percent)")
    
        mean, std = np.nanmean(avg_tracer_at_i), np.nanstd(avg_tracer_at_i)
        print("time", intervals[i], "average tracer", format(mean, ".2f"),  "pm", format(std, ".2f"), "mmol / L", )
    
    print("------------------------------------------------------------------------------")

    sleep_means = np.mean(patdf.loc[sleepers, :], axis=0)
    nosleep_means = np.mean(patdf.loc[nonsleep, :], axis=0)
    sleep_arr = np.array(patdf.loc[sleepers, :]).astype(float)
    no_sleep_arr = np.array(patdf.loc[nonsleep, :]).astype(float)

    nan_policy = "omit"

    print("-------------------------------------------------------------------------")
    print("ROI =", roi)
    stat, pval = scipy.stats.ttest_ind(sleep_arr, no_sleep_arr, axis=0, equal_var=True, nan_policy="omit")
    print("Students's t-test, p=", pval)

    stat, pval2 = scipy.stats.ttest_ind(sleep_arr[:, -1], no_sleep_arr[:, -1], axis=0, equal_var=True, nan_policy="omit")
    print("Students's t-test, p=", pval2, "(48 h)")

    print("-------------------------------------------------------------------------")

    sleep_std = stderr(np.array(patdf.loc[sleepers, :]).astype(float), axis=0)
    nosleep_std = stderr(np.array(patdf.loc[nonsleep, :]).astype(float), axis=0)

    cs = "tab:blue"
    cns = "tab:red"


    x = np.arange(len(labels))

    fig, ax = plt.subplots(dpi=dpi)

    fig.set_size_inches(size_inches[0], size_inches[1])
    rects1 = ax.bar(x - width / 2, sleep_means, width, color=cs, label='Sleep', linewidth=0.5, 
                    yerr=sleep_std, capsize=capsize_,)

    rects2 = ax.bar(x + width / 2, nosleep_means, width, color=cns, label='Sleep deprivation',
                    yerr=nosleep_std, capsize=capsize_,)

    displaystring = ""
    for i, p in enumerate(pval):
        
        if p >= 0.05 and only_significant:
            displaystring = "n.s."
            continue

        elif p < 0.0001:
            displaystring = r'***'
        elif p < 0.001:
            displaystring = r'**'
        else:
            displaystring = r'*'

        markersize = 5

        if roi == "avgds":
            dh = 0.005
            text_dh = 0.003

        elif roi == "avg":
            markersize = 4
            dh = 0.01
            text_dh = 0.01

        elif roi == "white":
            # dh = 0.004
            # text_dh = 0.004
            markersize = 4
            dh = 0.01
            text_dh = 0.01
        
        elif roi == "gray":
            # dh = 0.008
            # text_dh = 0.015
            dh = 0.005
            text_dh = 0.01

        height = 1 * (dh + 1 * \
            max(sleep_means[i] + sleep_std[i],
                nosleep_means[i] + nosleep_std[i]))

        bar_centers = x[i] + np.array([-width / 2., width / 2])
            
        significance_bar(bar_centers[0], bar_centers[1], height, fontsize=fontsize, markersize=markersize, displaystring=displaystring, text_dh=text_dh)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.set_ylabel(ylabel, fontsize=fontsize, loc="top")
    ax.set_xlabel("time (hours)", fontsize=fontsize)
    ax.set_xticks(x, labels, fontsize=fontsize)
    ax.tick_params(axis='y', which='major', labelsize=fontsize)

    if "surface" in ylabel or "ds" in ylabel:
        ax.set_ylim(0, 0.07)
    else:
        ax.set_ylim(0, 0.18)

    if argparse_dict["ylim"] is not None:
        ax.set_ylim(0, argparse_dict["ylim"])
        print("Overwriting default ylim")

    if roi == "white":
        ax.legend(fontsize=9)

    plt.locator_params(axis='y', nbins=4)

    fig.savefig(plotpath + roi + group + ".pdf")
    fig.savefig(plotpath + roi + group + ".png", dpi=600)

    statistics_path = plotpath
    patdf.to_csv(statistics_path + roi + group + ".csv", index=True)
    patdf.loc[sleepers, :].to_csv(statistics_path + roi + "_sleep.csv", index=True)
    patdf.loc[nonsleep, :].to_csv(statistics_path + roi + "_nonsleep.csv", index=True)

    plt.close()


plt.close()

