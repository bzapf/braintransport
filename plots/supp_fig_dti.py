import numpy as np
import os
import matplotlib.pyplot as plt
import json
import argparse
# from scripts.plots.read_data import extract_data
# from scripts.plot_helpers.analysis import lperrors
from definitions import intervals, groups
from helpers import load_experimental_data, load_simulation_data
from correlationplot import correlationplot
from scipy import stats

def stderr(x, axis): 
    return stats.sem(x, axis=axis, nan_policy="omit")


def plot_all(parameters, regions, pats, path_to_files, resultfoldernames, labels=None, savepath=None, 
        figsize=None, legendfs=12, fs=16, dpi=300, aspect=1,
        title=False, ylabel="", ylim=None, legend_outside=True, 
        ):

    n = 100

    average = True

    days = 2

    x = np.linspace(0, days, n)

    y_sim = np.zeros((n, len(parameters), len(pats))) - 1e16
    
    y_data = np.zeros((len(intervals), len(pats), len(parameters)))
    datatimes = np.zeros((len(intervals), len(pats), len(parameters)))

    legendfs = min(fs, legendfs)

    fig = plt.figure(figsize=figsize, dpi=dpi)

    ratios = np.zeros((len(intervals), len(pats), len(parameters))) + np.nan

    sims_at_measurement = np.zeros((len(intervals), len(pats), len(parameters))) + np.nan


    for pat_counter, pat in enumerate(pats):

        pat_parameters = []
        paths_to_plot = []

        for param_counter, p in enumerate(parameters):

            subfolder = os.path.join(path_to_files, pat, resultfoldernames[param_counter], "k" + str(iterk), "")

            paths_to_plot.append(subfolder)



        for param_counter, excelpath in enumerate(paths_to_plot):

            data_times, data = load_experimental_data(pat, excelpath, roi, intervals=intervals)

            y_data[:, pat_counter, param_counter] = data
            datatimes[:, pat_counter, param_counter] = days * np.array(data_times) / max(data_times)

            simtimes, simdata = load_simulation_data(pat, excelpath=os.path.join(excelpath, "concs_plain.csv"), roi=roi)
            
            assert max(simtimes) > 32
            assert max(simtimes)  < 60
            assert max(data_times) > 32
            assert max(data_times) < 60

            with open(excelpath + 'params') as data_file:
                hyperparameters = json.load(data_file)
            
            y_sim[:, param_counter, pat_counter] = simdata
            
            simtimes = np.array(simtimes)
            late_times = intervals[0:]
            for idt, (t1, t2) in enumerate(late_times):
                for idx, t in enumerate(data_times):
                    if t1 <= t and t <= t2:
                        # print("match", t1, t2, t)
                        argmin = np.argmin(np.abs(t-simtimes))

                        ratio = simdata[argmin] / data[idx]

                        sims_at_measurement[idt, pat_counter, param_counter] = simdata[argmin]

                        ratios[idt, pat_counter, param_counter] = ratio

            assert len(simdata.shape) == 1
            assert simdata.shape[0] == y_sim.shape[0]
            del simtimes, simdata


    del paths_to_plot, pat_parameters

    print("Ratios of simulated / measured tracer")
    # First axis is time, last axis is different T1 maps
    # breakpoint()
    ratios = np.nanmean(ratios, axis=1)
    print("T1 maps are", labels)
    for idt in range(ratios.shape[0]):
        print("time interval", late_times[idt], "ratios", ratios[idt, :])
    

    xd = np.nanmean(datatimes, axis=1,)
    yd = np.nanmean(y_data, axis=1)
    yerr = stderr(y_data, axis=1)
    xerr = stderr(datatimes, axis=1)

    for param_counter in range(len(parameters)):

        plt.errorbar(xd[..., param_counter], yd[..., param_counter], yerr=yerr[..., param_counter], xerr=xerr[..., param_counter],
            elinewidth=1, capsize=1, 
            marker="o", linewidth=0, markersize=12, color="k")

    # breakpoint()
    assert y_sim.min() >= - 1e12

    y = np.mean(y_sim, axis=2)
    yerr = stderr(y_sim, axis=2)

    assert np.sum(np.isnan(y_sim)) == 0

    lss =iter([":", "-"])

    for idx in range(len(parameters)):

        _label = ""

        lw = 2

        ls = next(lss)

        if labels is not None:
            _label = labels[idx]
              
        # plt.fill_between(x=x, y1=(y - yerr)[:, idx], y2=(y + yerr)[:, idx], alpha=0.5, color=colors[idx])
        plt.plot(x, y[:, idx], linestyle=ls, linewidth=lw, label=_label, color=colors[idx])

    # plot_night(fs=fs, max_t=50 / 2)

    print("Rel diff between lines", np.mean(np.abs(y[:, 0] - y[:, 1])) / np.mean(np.abs(y[:, 0])) )

    plt.plot([], [], label="measured", marker="o", linewidth=0, markersize=12, color="k")
    plt.legend(loc="upper left", fontsize=legendfs,)

    # if legend_outside is True:
    #     plt.legend(loc="upper left", fontsize=legendfs, bbox_to_anchor=(1.02, 0.8))
    # elif legend_outside == "no":
    #     pass  
    # else:
    #     plt.legend(loc="lower right", fontsize=fs, bbox_to_anchor=(0.5, 0.1, 0.5, 0.5))

    # plt.ylim(0, 0.16)

    plt.title("tracer in brain", fontsize=fs)
    
    plt.gca().set_box_aspect(aspect=aspect)

    plt.ylim(0, 0.3)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlabel("time (days)", fontsize=fs)
    plt.locator_params(axis='x', nbins=3)
    plt.locator_params(axis='y', nbins=5)
    plt.ylabel(ylabel, fontsize=fs)
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)

    return y_data, sims_at_measurement





if __name__ == "__main__":

    """

    Generate comparison plots for simulations based on dti and nodti

    """


    imageformat = ".png"

    parser = argparse.ArgumentParser()

    default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
    default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
    default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/dti/"

    parser.add_argument("--inputpath", type=str, default=default_input_path)
    parser.add_argument("--outputpath", type=str, default=default_output_path)

    parserargs = parser.parse_args()
    argparse_dict = vars(parserargs)


    data_folder = argparse_dict["inputpath"]
    plotpath = argparse_dict["outputpath"]

    if data_folder != default_input_path:
        assert plotpath != default_output_path

    os.chdir(data_folder)

    group = "all"
    pats = groups[group]

    group = "dti"
    pats = set(groups[group]).intersection(set(groups["t1map"]))

    pats = [str(x) for x in pats]

    print(len(pats), "have both t1map and dti")
    
    pats = sorted(pats)


    colors = ["tab:purple", "indigo", ]


    labels = ["simulated (DTI)", r"simulated ($\overline{D}$)"]

    resultfoldernames = ["diffusion/DTI_filteredT1", "diffusion/avgDTI_filteredT1"]

    iterk = 144

    roi = "avg" 
    ylabel = "concentration (mmol / L)"

    aspect = 0.8

    figsize_scale = 1

    dpi = 400

    fs = 22

    figsize = [6.4 * figsize_scale, figsize_scale * 4.8]

    y_data, y_sim = plot_all(parameters=[(1, 0, 0), (1, 0, 0), ], labels=labels, resultfoldernames=resultfoldernames,
            fs=fs, legendfs=fs - 2, dpi=dpi, aspect=aspect, figsize=figsize,
            regions=[roi], pats=pats, legend_outside=False, ylabel=ylabel,
            savepath=plotpath + "compare" + roi + ".png", path_to_files=data_folder)


    ylabel = "simulated (mmol / L)"
    xlabel = "measured (mmol / L)"
    
    figurename = plotpath + roi + "data-dti-avgd-correlation"

    labels = ["DTI ", r"$\overline{D}$ "]

    xlim = (0, 1.2 * max(np.max(np.nan_to_num(y_data)), np.max(np.nan_to_num(y_sim))))
    ylim = xlim

    correlationplot(
        x_data=np.expand_dims(y_data[..., -1], -1), y_data=y_sim, title="tracer in brain at ",
        fontsize=fs, labelsize=fs, legendsize=fs, aspect=aspect,
        labels=labels, xlim=xlim, ylim=ylim, dpi=dpi, figsize=figsize,
        ylabel=ylabel, colors=colors, xlabel=xlabel, roi=roi, figurepath=figurename)

    figurename = plotpath + roi + "dti-avgd-correlation"
    ylabel = "DTI simulation (mmol / L)"
    xlabel = "avg. D simulation (mmol / L)"

    correlationplot(
        x_data=np.expand_dims(y_sim[..., 0], -1), y_data=np.expand_dims(y_sim[..., -1], -1), title="tracer in brain at ",
        fontsize=fs, labelsize=fs, legendsize=fs, aspect=aspect,
        labels=labels, xlim=xlim, ylim=ylim, dpi=dpi, figsize=figsize,
        ylabel=ylabel, colors=["k"], xlabel=xlabel, roi=roi, figurepath=figurename)