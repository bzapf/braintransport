import numpy as np
import os
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import json

import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap
from typing import Callable, Union


# from helpers import get_data_in_intervals, sim_at_mri_times
from definitions import intervals # , reaction_resultfolder
from helpers import extract_data


def make_figs(region, pats, alphas, data_folder, average_tracer=False):
    n_t = 4  # number of time points
    n_a = len(alphas)  # number of alphas
    n_p = len(pats)  # number of patients
    conc_simulated = np.zeros((n_t, n_p, n_a)) + np.nan
    conc_experimental = np.zeros((n_t, n_p))  + np.nan

    
    # cmap_name = 'jet'
    # cmap = LinearSegmentedColormap.from_list(cmap_name, colorlist, N=100)

    cmap = matplotlib.colormaps["jet"].resampled(len(pats))(np.linspace(0, 1, len(pats)))
    cmap = iter([cmap[i, :] for i in range(len(pats))])    

    for pat in pats:

        fig, ax = plt.subplots()

        col = next(cmap)

        for alpha_idx, alpha in enumerate(alphas):

            resolution = 32
        
            # if "alphatest" in resultfoldername(pat) and not np.isnan(alpha):
            #     folder = data_folder + str(pat) + "/" + resultfoldername(pat) + "/alpha" + str(alpha) + "/"
            
            # else:
            #     assert np.isnan(alpha)
            raise NotImplementedError # TODO check and use correct folder
            folder = reaction_resultfolder(pat, best=True)
            assert folder is not None
            assert os.path.isdir(folder)

            print(folder)


            params = json.load(open(folder + "params"))

            assert params["concfolder"] == "FIGURES_CONC"

            if region == "avgds":
                with open(data_folder + str(pat) + '/region_areas' + str(resolution), 'rb') as f:
                    region_volumes = pickle.load(f)
                    roi_volume = region_volumes[region] / 1e6
            else:                
                with open(data_folder + str(pat) + '/region_volumes' + str(resolution), 'rb') as f:
                    region_volumes = pickle.load(f)
                    roi_volume = region_volumes[region] / 1e6

            if average_tracer:
                roi_volume = 1e-6

            experimental = pd.read_csv(folder + 'experimental_data.csv')

            if np.isnan(alpha):
                concs = pd.read_csv(folder + 'concs.csv')
            else:
                concs = pd.read_csv(folder + 'concs_plain.csv') 
            
            simulation_times = concs['t'] # / 3600
            simulated_tracer = concs[region] * roi_volume

            assert max(np.abs(simulated_tracer)) < 1
            
            experimental_times = experimental['t'] # / 3600
            measured_tracer = experimental[region] * roi_volume

            ls = "-"
            if np.isnan(alpha):
                ls = "--"

            plt.title(pat + " " + region)
            ax.plot(simulation_times / 3600, simulated_tracer, linestyle=ls, color=col)
            ax.plot(experimental_times / 3600, measured_tracer, linestyle="-", linewidth="0", marker="x", color=col)

            _, measured_tracer_at_times = get_data_in_intervals(pat, simtimes=experimental_times, simdata=measured_tracer, intervals=intervals)
            _, simulated_tracer_at_times = get_data_in_intervals(pat, simtimes=simulation_times, simdata=simulated_tracer, intervals=intervals)


    plt.show()
    exit()









def make_barplot(region, pats, alphas, paperformat, resultfoldername: Callable, data_folder, savepath: Union[None, str], fs, figsize, dpi, ylabel, 
                GREY_WHITE=False, width=None, print_format=".2f",
                average_tracer=False):


    fig, ax = plt.subplots()

    conc_experimental, conc_simulated = extract_data(alphas, pats, data_folder, resultfoldername, region, average_tracer, intervals)

    mridata, mridata_standarddev = np.nanmean(conc_experimental, axis=1), np.nanstd(conc_experimental, axis=1)

    simdata, simdata_standarddev = np.nanmean(conc_simulated, axis=1), np.nanstd(conc_simulated, axis=1)

    ratios = conc_simulated / np.expand_dims(conc_experimental, axis=-1)

    ratios, ratios_std = np.nanmean(ratios, axis=1), np.nanstd(ratios, axis=1)

    print("Simulated measured tracer level")
    print("time  ", "      2h   ", "            6h     ", "            24h       ", "         48h   ")
    print("data   ", [format(x, print_format) + r""" pm """ + format(stdx, print_format) for x, stdx in zip(mridata, mridata_standarddev)]) #, mridata_standarddev / np.sqrt(n_p))

    for i, alpha in enumerate(alphas):


        print("alpha", alpha,  [format(x, print_format) + " pm " + format(stdx, print_format) for x, stdx in zip(simdata[:, i], simdata_standarddev[:, i])])

    print("Ratios between simulation and data")
   
    print("time                         ", "      2h   ", "            6h     ", "        24h    ", "    48h   ")    
    for i, alpha in enumerate(alphas):
        print("alpha=", alpha, "Ratios sim/data =",  [format(x, print_format) + " pm " + format(stdx, print_format) for x, stdx in zip(ratios[:, i], ratios_std[:, i])])

    if "none" in savepath.lower():
        print("'none' found in <savepath>, only printing info and exiting make_barplot()")
        return


    # BZ Generalization to more alphas, populate list and skip 2 (to plot exp data and have space between groups)
    k = 0
    kk = 0
    x_sim_ = []
    x_e_s = []
    xticks_ = [kk]

    for j in range(simdata.shape[0]):
        if j > 0:
            k += 2
            kk += 1
        x_e_s.append(k)
        for i in range(len(alphas)):
            k += 1
            kk += 1
            x_sim_.append(k)
            xticks_.append(kk)

        if j > 0:
            kk += 1
            xticks_.append(kk)
    
    x_sim = np.array(x_sim_)
    y_sim = simdata.flatten()
    simdata_standarddev = simdata_standarddev.flatten()

    x_e = np.array(x_e_s)
    y_e = mridata.flatten()
    mridata_standarddev = mridata_standarddev.flatten()


    capsize_ = 2
    k = 2

    if region == "avg":
        colors = ["k" for _ in range(10)]

    elif region == "avgds":
        colors = ["magenta" for _ in range(10)]

    elif region == "white":
        colors = np.array([[0.45176471, 0.76708958, 0.46120723, 1.        ],
                            [0.25259516, 0.66812764, 0.36286044, 1.        ],
                            [0.13402537, 0.54232987, 0.26828143, 1.        ],
                            [0.        , 0.42303729, 0.17071895, 1.        ],
                            [0.        , 0.26666667, 0.10588235, 1.        ]])

        if GREY_WHITE:
            colors = ["white" for x in range(5)]

    elif region == "gray":
        colors = np.array([[0.41708574, 0.68063053, 0.83823145, 1.        ],
                            [0.25628604, 0.57001153, 0.7751634 , 1.        ],
                            [0.12710496, 0.44018454, 0.70749712, 1.        ],
                            [0.03137255, 0.31409458, 0.60648981, 1.        ],
                            [0.03137255, 0.18823529, 0.41960784, 1.        ]])

        if GREY_WHITE:
            colors = ["gainsboro", "silver", "darkgrey", "grey", "dimgrey"]

    datacolor = colors[0]

    colorlist = colors

    colors = iter(colors)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # if np.nan in alphas:
    label="data"
    # else:
    #     label = None

    datahatch = '...'
    edgecolor = None

    if GREY_WHITE:
        datahatch='...'
        edgecolor="k"

    ax.bar(x_e, y_e, yerr=mridata_standarddev / 1, # np.sqrt(n_p), 
            capsize=capsize_, width=width, color=datacolor, label=label, hatch=datahatch, edgecolor=edgecolor)

    for i, alpha in enumerate(alphas):    
        hatch = None
        label = None

        if len(alphas) == 1: # or (i == 0 and len(alphas) == 5):
            label = "simulation"
        
        elif len(alphas) == 2 and alpha == 1:
            label = "extracellular diffusion"
        
        elif np.isnan(alpha):
            hatch = "xx"
            label = "enhanced diffusion\n& local clearance"
            
        color = next(colors)

        ax.bar(x_sim[i::(len(alphas))], y_sim[i::(len(alphas))], yerr=simdata_standarddev[i::(len(alphas))] / 1, #np.sqrt(n_p), 
                edgecolor=edgecolor, width=width, 
                capsize=capsize_, label=label, color=color, hatch=hatch)


    ax.set_xticks([int(len(alphas)/2) + (len(alphas) + 2) * x for x in range(simdata.shape[0])],)
    ax.set_xticklabels([r"$\sim 2\,$", r"$\sim 6\,$", r"$\sim 24\,$", r"$\sim 48\,$",], fontsize=fs)

    ax.tick_params(axis='x', width=0)
    ax.tick_params(axis='y', labelsize=fs)

    ax.set_xlabel("time (hours)", fontsize=fs)
    
    
    if np.nan in alphas:
        plt.ylim(-0., 0.24)
        ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21])

    else:
        plt.ylim(-0.00, 0.20)
        ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18])
    
    if np.nan in alphas:
        plt.legend(fontsize=matplotlib.rcParams["legend.fontsize"]-2)



    if paperformat:
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

    if len(alphas) > 3:
        # x = np.linspace(0, 3 * np.pi, 500)
        # y = np.sin(x)
        # dydx = np.cos(0.5 * (x[:-1] + x[1:]))
        # print(dydx.shape)
        # points = np.array([x, y]).T.reshape(-1, 1, 2)
        # segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # # print(segments.shape)

        plt.legend(fontsize=matplotlib.rcParams["legend.fontsize"]-2)
        
        dydx = np.zeros(10)
        segments = np.zeros((10, 2, 2)) + np.nan

        cmap_name = 'my_list'
        cmap = LinearSegmentedColormap.from_list(cmap_name, colorlist, N=100)

        norm = plt.Normalize(min(alphas), max(alphas))
        lc = LineCollection([], cmap=cmap, norm=norm)
        # Set the values used for colormapping
        lc.set_array(dydx)
        lc.set_linewidth(0)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax, shrink=0.5)
        cbar.ax.set_xlabel(#"dispersion\nfactor\n"+
                            r"$\alpha$", rotation=0)
        cbar.ax.set_ylabel("simulation", rotation=270, labelpad=40)

        cbar.ax.set_yticks([1, 3, 5])

    if region == "white":
        ylabel = ylabel + "   "
        loc = "center"
    else:
        loc = "center"

    ax.set_ylabel(ylabel, # "tracer in " + region + " (mmol)", 
                    fontsize=fs, loc=loc)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)