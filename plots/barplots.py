import numpy as np
import os
import pandas as pd
import pickle
import itertools
import json

if "cluster" not in os.getcwd():

    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.colors import LinearSegmentedColormap
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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









def make_barplot(region, pats, alphas, paperformat, resultfoldername: Callable, data_folder, savepath: Union[None, str], FS, figsize, dpi, ylabel, 
                GREY_WHITE=False, width=None, print_format=".2f",
                average_tracer=False):


    

    conc_experimental, conc_simulated = extract_data(alphas, pats, data_folder, resultfoldername, region, average_tracer, intervals)

    mridata, mridata_standarddev = np.nanmean(conc_experimental, axis=1), np.nanstd(conc_experimental, axis=1)

    simdata, simdata_standarddev = np.nanmean(conc_simulated, axis=1), np.nanstd(conc_simulated, axis=1)

    ratios = conc_simulated / np.expand_dims(conc_experimental, axis=-1)

    ratios, ratios_std = np.nanmean(ratios, axis=1), np.nanstd(ratios, axis=1)

    if average_tracer:
        unit = "(mmol/L)"
    else:
        unit = "(mmol)"
    df = {}

    print("Simulated measured tracer level")
    columns = [""] + ["alpha"] + list(itertools.chain.from_iterable([["mean (" + x + ")", "std. (" + x + ")"] for x in ["2h", "6h", "24h", "48h"]]))
    df[0] = ["data " + unit] + ["-"] + list(itertools.chain.from_iterable([[format(x, print_format), format(stdx, print_format)] for x, stdx in zip(mridata, mridata_standarddev)]))
    
    df[1] = ["data (%)"] + ["-"] + list(itertools.chain.from_iterable([[format(100 * x / 0.5, ".1f"), format(100 * stdx / 0.5, ".1f")] for x, stdx in zip(mridata, mridata_standarddev)]))

    for i, alpha in enumerate(alphas):
        m = max(df.keys())

        df[m + 1] = ["" for x in df[m]]
        
        m = max(df.keys())

        df[m + 1] = ["simulated " + unit] + [alpha] 
        df[m + 1] += list(itertools.chain.from_iterable([[format(x, print_format), format(stdx, print_format)] for x, stdx in zip(simdata[:, i], simdata_standarddev[:, i])]))
        
        df[m + 2] = ["simulated (%)"] + [alpha] 
        df[m + 2] += list(itertools.chain.from_iterable([[format(100 * x / 0.5, ".1f"), format(100 * stdx / 0.5, ".1f")] for x, stdx in zip(simdata[:, i], simdata_standarddev[:, i])]))
        
        df[m + 3] = ["ratio sim/data"] + [alpha]
        df[m + 3] += list(itertools.chain.from_iterable([[format(100 * x, ".1f"), format(100 * stdx, ".1f")] for x, stdx in zip(ratios[:, i], ratios_std[:, i])]))



    df = pd.DataFrame.from_dict(df, orient="index", columns=columns)

    print("-" * 120)
    print("Region:", region, "average tracer =", average_tracer)
    print(df)
    print("-" * 120)

    if "none" in savepath.lower():
        print("'none' found in <savepath>, only printing info and exiting make_barplot()")
        return

    
    # exit()

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
    # y_sim = simdata.flatten()
    # simdata_standarddev = simdata_standarddev.flatten()

    x_e = np.array(x_e_s)
    y_e = mridata.flatten()
    mridata_standarddev = mridata_standarddev.flatten()


    capsize_ = 4
    k = 2

    if region == "avg":
        colors = ["k" for _ in range(10)]
        hatches = [None for x in range(5)]

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
            hatches = ["."*(i-1) for i in range(1, 6)]

    elif region == "gray":
        colors = np.array([[0.41708574, 0.68063053, 0.83823145, 1.        ],
                            [0.25628604, 0.57001153, 0.7751634 , 1.        ],
                            [0.12710496, 0.44018454, 0.70749712, 1.        ],
                            [0.03137255, 0.31409458, 0.60648981, 1.        ],
                            [0.03137255, 0.18823529, 0.41960784, 1.        ]])

        if GREY_WHITE:
            colors = ["gainsboro", "silver", "darkgrey", "grey", "dimgrey"]
            hatches = [None for x in range(5)]

    datacolor = colors[0]

    colorlist = colors

    colors = iter(colors)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # if np.nan in alphas:
    label="data"
    # else:
    #     label = None

    datahatch = 'x'
    edgecolor = None

    if GREY_WHITE:
        datahatch = 'x'
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
        
        elif len(alphas) > 3:
            hatch=hatches[i]

        color = next(colors)


        ax.bar(x_sim[i::(len(alphas))], 
                simdata[:, i], yerr=simdata_standarddev[:, i],
                # y_sim[i::(len(alphas))], yerr=simdata_standarddev[i::(len(alphas))] / 1, #np.sqrt(n_p), 
                edgecolor=edgecolor, width=width, 
                capsize=capsize_, label=label, color=color, hatch=hatch)

    # if region == "white":
    #     breakpoint()

    ax.set_xticks([int(len(alphas)/2) + (len(alphas) + 2) * x for x in range(simdata.shape[0])],)
    ax.set_xticklabels([r"$\sim 2\,$", r"$\sim 6\,$", r"$\sim 24\,$", r"$\sim 48\,$",], fontsize=None)

    ax.tick_params(axis='x', width=0)
    ax.tick_params(axis='y', labelsize=None)

    ax.set_xlabel("time (hours)", fontsize=None)
    
    
    if np.nan in alphas:
        if region == "gray":
            plt.ylim(-0., 0.24)
            ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21])
        if region == "white":
            plt.ylim(-0., 0.12)
            ax.set_yticks([0., 0.03, 0.06, 0.09, ])
            
    else:
        if region == "gray":
            plt.ylim(-0.00, 0.20)
            ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18])
    
        if region == "white":
            plt.ylim(-0.00, 0.11)
            ax.set_yticks([0., 0.03, 0.06, 0.09])


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

        
        
        if region == "avg":

            plt.legend(fontsize=matplotlib.rcParams["legend.fontsize"]-2, frameon=True)
            # old version:
            dydx = np.zeros(10)
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
        
        # New version:

        else:

            plt.legend(fontsize=matplotlib.rcParams["legend.fontsize"]-2, frameon=False, loc="upper left")
            
            # if region == "white":
            cbaxes = inset_axes(ax, width="40%", height="4%", 
                loc="upper left",
                bbox_to_anchor=(0.03, 0., 1, 0.84),
                bbox_transform=ax.transAxes,
                borderpad=0,
                
                )
            # elif region == "gray":
            #     cbaxes = inset_axes(ax, width="4%", height="40%", loc="center left")



            dydx = np.zeros(10)

            cmap_name = 'my_list'
            cmap = LinearSegmentedColormap.from_list(cmap_name, colorlist, N=100)

            norm = plt.Normalize(min(alphas), max(alphas))
            lc = LineCollection([], cmap=cmap, norm=norm)
            
            # Set the values used for colormapping
            lc.set_array(dydx)
            lc.set_linewidth(0)

            line = cbaxes.add_collection(lc)

            if region == "white":
                cbar = plt.colorbar(line, cax=cbaxes, shrink=0.5, orientation='horizontal')
                # print(ax.get_xlim())
                # print(alphas)
                # print([int(len(alphas)/2) + (len(alphas) + 2) * x for x in range(simdata.shape[0])])
                # breakpoint()
                # cbaxes.bar([0.2, 0.23], [0,2], width=0.5, color="k")
                # cbaxes.plot(np.linspace(2, 4, 10), np.linspace(0.05, 0.06, 10),color="red", markersize=100, marker="x")

                # nx, ny = 4, 4
                # xx, yy = np.meshgrid(np.linspace(0,1, nx), np.linspace(0,1, ny))                

                # xy = np.reshape(np.stack((xx, yy)), (2, nx*ny))

                # xy[1, :] += 1


                # cbaxes.scatter(xy[0, :], xy[1, :])
                # breakpoint()
                for idx, h in enumerate(hatches):
                    cbaxes.bar(height=1, width=4/5, x=1.5 + 4/5 * idx, # bottom=1 + idx, 
                            hatch=h, color="white")



                # pass

            elif region == "gray":
                cbar = plt.colorbar(line, cax=cbaxes, shrink=0.5, orientation='horizontal')
                # cbar.ax.set_xlabel(r"$\alpha$", rotation=0)
                # cbar.ax.set_ylabel("simulation", rotation=0, labelpad=40)

                # cbar.ax.set_yticks([1, 3, 5])

            cbar.ax.set_xticks([1.5 + 4/5 * idx for idx in range(5)], range(1,6) )


            cbar.ax.set_xlabel(r"$\alpha$", rotation=0)
            cbar.ax.set_title("simulation", rotation=0, fontsize=FS) # , labelpad=40)
            
        # cbar.ax.set_yticks([1, 3, 5])

    if region == "white":
        ylabel = ylabel + "   "
        loc = "center"
    else:
        loc = "center"

    ax.set_ylabel(ylabel, # "tracer in " + region + " (mmol)", 
                    fontsize=None, loc=loc)

    # ax.set_title(ylabel, # "tracer in " + region + " (mmol)", 
    #                 fontsize=None, loc=loc)

    

    plt.tight_layout()


    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)