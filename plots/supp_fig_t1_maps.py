import numpy as np
import os
import matplotlib.pyplot as plt
import json
import argparse
import pickle
import scipy
from scipy import stats
import pandas

from correlationplot import correlationplot
from helpers import get_data_in_intervals
from definitions import labels as data_labels, groups, intervals
from helpers import load_experimental_data, load_simulation_data

def stderr(x, axis): 
    return stats.sem(x, axis=axis, nan_policy="omit")

def method_barplots(parameters, roi, pats, path_to_files, resultfoldernames, labels=None, savepath=None, 
        savedpi=300, figsize=(6, 5), legendfs=12, 
        title=False, ylabel="", ylim=None, legend_outside=True, 
        ):

    n = 100

    average = True

    if not average:
        assert "avg" not in ylabel



    y_sim = np.zeros((n, len(parameters), len(pats))) - 1e16
    
    y_data = np.zeros((len(intervals), len(pats), len(parameters)))
    datatimes = np.zeros((len(intervals), len(pats), len(parameters)))

    print_from = 1
    late_times = intervals[print_from:]

    fontsize = 16
    legendfs = min(fontsize, legendfs)

    ratios = np.zeros((3, len(pats), len(parameters))) + np.nan

    for pat_counter, pat in enumerate(pats):

        pat_parameters = []
        paths_to_plot = []


        if roi != "avgds":
            with open(os.path.join(path_to_files, pat, 'region_volumes32'), 'rb') as f:
                region_volumes = pickle.load(f)
            for r in roi:
                if "ds" in roi:
                    raise KeyError
        else:
            print("Loading surface areas")
            assert roi == "avgds"
            with open(os.path.join(path_to_files, pat, 'region_areas32'), 'rb') as f:
                region_volumes = pickle.load(f)

        v = region_volumes[roi] * 1e-6



        for param_counter, p in enumerate(parameters):

            subfolder = os.path.join(path_to_files, pat, resultfoldernames[param_counter], "k" + str(iterk), "")

            paths_to_plot.append(subfolder)


        for param_counter, excelpath in enumerate(paths_to_plot):

            data_times, data = load_experimental_data(pat, excelpath, roi, intervals)

            data *= v

            y_data[:, pat_counter, param_counter] = data
            datatimes[:, pat_counter, param_counter] = days * np.array(data_times) / max(data_times)

            if "avgds" == roi:
                continue

            simtimes, simdata = load_simulation_data(pat, excelpath=os.path.join(excelpath, "concs_plain.csv"), roi=roi)
            
            simdata *= v

            y_sim[:, param_counter, pat_counter] = simdata
            
            simtimes = np.array(simtimes)

            _, measured_tracer_at_times = get_data_in_intervals(pat, stored_times=data_times, stored_data=data, intervals=intervals)
            _, simulated_tracer_at_times = get_data_in_intervals(pat, stored_times=simtimes, stored_data=simdata, intervals=intervals)

            ratios[:, pat_counter, param_counter] = (simulated_tracer_at_times / measured_tracer_at_times)[print_from:]
            
            assert len(simdata.shape) == 1
            assert simdata.shape[0] == y_sim.shape[0]
            del simtimes, simdata

    del paths_to_plot, pat_parameters

    print("Ratios of simulated / measured tracer")

    ratios = np.nanmean(ratios, axis=1)
    print("T1 maps are", labels)
    for idt, _ in enumerate(late_times):
        print("time interval", late_times[idt], "ratios", ratios[idt, :])

    mean_times = np.nanmean(datatimes, axis=1,)
    yd = np.nanmean(y_data, axis=1)
    yerr = stderr(y_data, axis=1)
    xerr = stderr(datatimes, axis=1)


    for tid in range(y_data.shape[0]):
        c1 = y_data[tid, :, 0]
        c2 = y_data[tid, :, 1]
        c3 = y_data[tid, :, 2]

        _, pval = scipy.stats.ttest_ind(c1, c2, axis=0, equal_var=True, nan_policy="omit")
        if pval < 0.05:
            print("c1 vs c2, p=", pval, "(", roi[0], ")")

        _, pval = scipy.stats.ttest_ind(c2, c3, axis=0, equal_var=True, nan_policy="omit")
        if pval < 0.05:
            print("c2 vs c3, p=", pval, "(", roi[0], ")")

        _, pval = scipy.stats.ttest_ind(c1, c3, axis=0, equal_var=True, nan_policy="omit")
        if pval < 0.05:
            print("c1 vs c3, p=", pval, "(", roi[0], ")")


    ts = np.array(list(range(mean_times.shape[0])))

    width = 1 / (len(labels) + 1) # the width of the bars

    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    plot_from = 0

    assert len(ts[plot_from:]) == len(data_labels)

    for i in range(len(parameters)):

        rects1 = ax.bar(ts[plot_from:] - width + i * width, yd[plot_from:, i], width, color=colors[i], linewidth=0.5, 
                        yerr=yerr[plot_from:, i], capsize=2, label=labels[i])

    maxdiff12, maxdiff23 = 0, 0
    maxstd12, maxstd23 = 0, 0

    for id1 in range(yd.shape[0]):
        
        c1 = y_data[id1, :, 0]
        c2 = y_data[id1, :, 1]
        c3 = y_data[id1, :, 2]

        print("time point", id1, )
        diff12 = np.nanmean(np.abs(c1-c2) / np.abs(c1))
        std12 = np.nanstd(np.abs(c1-c2) / np.abs(c1))
        
        print("diff between 1 and 2", diff12, "+-", std12)
        
        if diff12 > maxdiff12:
            maxdiff12 = diff12
            maxstd12 = std12
            maxid12 = id1
        diff23 = np.nanmean(np.abs(c2-c3) / np.abs(c2))
        std23 = np.nanstd(np.abs(c3-c2) / np.abs(c2))
        
        print("diff between 2 and 3", diff23, "+-", std23)
        
        if diff23 > maxdiff23:
            maxdiff23 = diff23
            maxstd23 = std23
            maxid23 = id1

    print("1-2: maximum difference at tid",  maxid12, maxdiff12, "+-", maxstd12)
    print("2-3: maximum difference at tid",  maxid23, maxdiff23, "+-", maxstd23)

    if "ds" in roi[0]:
        ax.set_ylim((0, 0.08))
        if "avg" in roi[0]:
            ax.set_ylim((0, 0.5))
    else:
        ax.set_ylim((0, 0.2))

    plt.locator_params(axis='y', nbins=4)
    ax.set_ylabel(ylabel, fontsize=fontsize, loc="top")
    ax.set_xlabel("time (hours)", fontsize=fontsize)
    ax.set_xticks(ts[plot_from:], data_labels, fontsize=fontsize)
    
    ax.tick_params(axis='x', which='major', length=0, width=0)
    ax.tick_params(axis='y', which='major', labelsize=fontsize)
    ax.legend(fontsize=fontsize,
        bbox_to_anchor=(0.01, 0.7 ), loc='lower left',
                        ncol=1, # mode="expand", 
                        borderaxespad=0)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.tight_layout()
    plt.savefig(plotpath + roi[0] + "concentrations.png")

    # if "ds" in regions[0]:
    return y_data, y_sim, datatimes



def plot_lines(y_data, y_sim, datatimes, savepath, savedpi, figsize, parameters):
    
    x = np.linspace(0, days, y_sim.shape[0])
    
    mean_times = np.nanmean(datatimes, axis=1,)
    yd = np.nanmean(y_data, axis=1)
    yerr = stderr(y_data, axis=1)
    xerr = stderr(datatimes, axis=1)

    fig = plt.figure(figsize=figsize, dpi=200)
    
    for param_counter in range(len(parameters)):

        plt.errorbar(mean_times[..., param_counter], yd[..., param_counter], yerr=yerr[..., param_counter], xerr=xerr[..., param_counter],
            elinewidth=1, capsize=1,
            marker="o", linewidth=0, markersize=12, color=colors[param_counter])

    assert y_sim.min() >= - 1e12

    y = np.mean(y_sim, axis=2)
    yerr = stderr(y_sim, axis=2)

    assert np.sum(np.isnan(y_sim)) == 0

    for idx in range(len(parameters)):

        _label = ""

        lw = 2

        if labels is not None:
            _label = labels[idx]
              
        plt.fill_between(x=x, y1=(y - yerr)[:, idx], y2=(y + yerr)[:, idx], alpha=0.5, color=colors[idx])
        plt.plot(x, y[:, idx], linewidth=lw, label=_label, color=colors[idx])

        print("Ratio sim/data after 6 h", y[-1, idx] / yd[-1, idx])
        print("Ratio sim/data after 48 h", y[-1, idx] / yd[-1, idx])

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlabel("time (days)", fontsize=fontsize)
    plt.locator_params(axis='x', nbins=3)
    plt.locator_params(axis='y', nbins=5)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=savedpi)





def load_sirs(roi):

    sirs = np.zeros((4, len(pats), 1))

    # sirs[:, :, :-1] = y_data

    for patid, pat in enumerate(pats):
        path2sir = os.path.join(data_folder, pat, "averaged_SIR", "experimental_data.csv")

        assert os.path.isfile(path2sir)

        exceltable = pandas.read_csv(path2sir)

        _, avgsir = get_data_in_intervals(pat, stored_times=exceltable["t"], stored_data=exceltable[roi], intervals=intervals)

        if roi != "avgds":
            with open(os.path.join(data_folder, pat, 'region_volumes' + str(32)), 'rb') as f:
                region_volumes = pickle.load(f)
            for r in roi:
                if "ds" in roi:
                    raise KeyError
        else:
            print("Loading surface areas")
            assert roi == "avgds"
            with open(os.path.join(data_folder, pat, 'region_areas' + str(32)), 'rb') as f:
                region_volumes = pickle.load(f)

        v = region_volumes[roi] * 1e-6

        # breakpoint()
        sirs[:, patid, -1] = avgsir * v

    return sirs


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
    default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
    default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/t1maps/"

    parser.add_argument("--inputpath", type=str, default=default_input_path)
    parser.add_argument("--outputpath", type=str, default=default_output_path)

    parserargs = parser.parse_args()
    argparse_dict = vars(parserargs)


    data_folder = argparse_dict["inputpath"]
    plotpath = argparse_dict["outputpath"]

    if data_folder != default_input_path:
        assert plotpath != default_output_path
    imageformat = ".png"

    fontsize = 16

    days = 2

    group = "t1map"
    pats = groups[group]
    sleepers = groups["sleep"]
    nonsleep = groups["sleepdep"]

    sleepers = [str(x) for x in sleepers if x in pats]
    nonsleep = [str(x) for x in nonsleep if x in pats]
    pats = nonsleep + sleepers
    pats = sorted(pats)

    groupids = {"sleep": sleepers, "nosleep": nonsleep}

    iterk = 144

    colors = ["tab:orange", "tab:green", "tab:blue", "tab:red"]

    labels=[r"raw $T_1$ map",      
            r"filtered $T_1$ map", 
            r"$\overline{T}_1$ in w/g/b",
            ]

    resultfoldernames = ["diffusion/avgDTI_T1", "diffusion/avgDTI_filteredT1", "diffusion/avgDTIavgT1"]

    regionnames = {"avg": "in brain", "avgds": "on surface"}


    for roi, ylabel in zip(["avg", "avgds"], 
                    ["avg. tracer in brain (mmol/L)", 
                    "avg. tracer on surface (mmol / L)"
                    ]):
        
        print("----------------------------------------------------------------------------------------------------------------")
        print(roi, ylabel)
        print("----------------------------------------------------------------------------------------------------------------")


        y_data, y_sim, datatimes = method_barplots(parameters=[(1, 0, 0) for _ in resultfoldernames], labels=labels, resultfoldernames=resultfoldernames,
                roi=roi, pats=pats, legend_outside=False, ylabel=ylabel,
                savepath=plotpath + "compare" + roi + ".png", path_to_files=data_folder)
        
        if roi == "avg":

            plot_lines(y_data, y_sim, datatimes, savepath=plotpath + "lines_" + roi + ".png", savedpi=600, figsize=(6, 5), parameters=[(1, 0, 0) for _ in resultfoldernames])


        sirs = load_sirs(roi)

        
        aspect = 0.8
        figsize_scale = 1
        dpi = 400
        fs = 22
        figsize = [6.4 * figsize_scale, figsize_scale * 4.8]

        ylabel = "mmol"
        
        figurename = plotpath + roi + "t1-correlation-t"
        xlabel = "tracer (mmol) \n(computed with " + labels[0] + ")"
        
        xlim = (0, 1.2 * max(np.max(np.nan_to_num(y_data)), np.max(np.nan_to_num(y_sim))))
        ylim = xlim

        correlationplot(
            x_data=np.expand_dims(y_data[..., 0], -1),
            y_data=y_data[..., 1:], # np.expand_dims(y_sim[..., -1], -1), 
            title="tracer in brain at ",
            fontsize=fs, labelsize=fs, legendsize=fs, aspect=aspect,
            labels=["" for r in labels], xlim=xlim, ylim=ylim, dpi=dpi, figsize=figsize,
            ylabel=ylabel, colors=colors[1:], xlabel=xlabel, roi=roi, figurepath=figurename)


        
        figurename = plotpath + roi + "sir-correlation-t"
        xlabel = r"avg. signal increase " + regionnames[roi] + r" (%)"
        xlim = (0, 1.2 * max(np.max(np.nan_to_num(sirs)), np.max(np.nan_to_num(sirs))))
        ylim = (0, 1.2 * max(np.max(np.nan_to_num(y_data)), np.max(np.nan_to_num(y_sim))))


        correlationplot(
            x_data=sirs, 
            y_data=y_data, # np.expand_dims(y_sim[..., -1], -1), 
            title="tracer in brain at ",
            fontsize=fs, labelsize=fs, legendsize=fs, aspect=aspect,
            labels=["" for r in labels], xlim=xlim, ylim=ylim, dpi=dpi, figsize=figsize,
            ylabel=ylabel, colors=colors, xlabel=xlabel, roi=roi, figurepath=figurename)