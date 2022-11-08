import numpy as np
import os
import pandas as pd
import pickle
from definitions import intervals, reaction_resultfolder
import matplotlib.pyplot as plt
import json
from helpers import get_data_in_intervals



def make_barplot(region, pats, alphas, paperformat, resultfoldername, data_folder, savepath, fs, figsize, dpi, average_tracer=False):
    n_t = 4  # number of time points
    n_a = len(alphas)  # number of alphas
    n_p = len(pats)  # number of patients
    conc_simulated = np.zeros((n_t, n_p, n_a)) + np.nan
    conc_experimental = np.zeros((n_t, n_p))  + np.nan

    for alpha_idx, alpha in enumerate(alphas):
        
        pat_no = 0

        for pat in pats:

            resolution = 32
        
            if "alphatest" in resultfoldername(pat) and not np.isnan(alpha):
                folder = data_folder + str(pat) + "/" + resultfoldername(pat) + "/alpha" + str(alpha) + "/"
            
            else:
                assert np.isnan(alpha)
                folder = reaction_resultfolder(pat, best=True)
                assert folder is not None
                assert os.path.isdir(folder)

            params = json.load(open(folder + "params"))

            assert params["concfolder"] == "FIGURES_CONC"

            with open(data_folder + str(pat) + '/region_volumes' + str(resolution), 'rb') as f:
                region_volumes = pickle.load(f)
                roi_volume = region_volumes[region] / 1e6

            if average_tracer:
                roi_volume = 1e-6

            experimental = pd.read_csv(folder + 'experimental_data.csv')
            concs = pd.read_csv(folder + 'concs.csv') 
            
            simulation_times = concs['t'] # / 3600
            simulated_tracer = concs[region] * roi_volume

            assert max(np.abs(simulated_tracer)) < 1
            
            experimental_times = experimental['t'] # / 3600
            measured_tracer = experimental[region] * roi_volume

            _, measured_tracer_at_times = get_data_in_intervals(pat, stored_times=experimental_times, stored_data=measured_tracer, intervals=intervals)
            _, simulated_tracer_at_times = get_data_in_intervals(pat, stored_times=simulation_times, stored_data=simulated_tracer, intervals=intervals)

            conc_experimental[:, pat_no] = measured_tracer_at_times
            conc_simulated[:, pat_no, alpha_idx] = simulated_tracer_at_times   

            pat_no += 1

    e, st_e = np.nanmean(conc_experimental, axis=1), np.nanstd(conc_experimental, axis=1)

    sim, st_sim = np.nanmean(conc_simulated, axis=1), np.nanstd(conc_simulated, axis=1)

    # BZ Generalization to more alphas, populate list and skip 2 (to plot exp data and have space between groups)
    k = 0
    kk = 0
    x_sim_ = []
    x_e_s = []
    xticks_ = [kk]

    for j in range(n_t):
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
    y_sim = sim.flatten()
    st_sim = st_sim.flatten()

    x_e = np.array(x_e_s)
    y_e = e.flatten()
    st_e = st_e.flatten()

    names = ['exp'] + [r'$\alpha$ = %d' % i for i in range(1, max(alphas) + 1)]
    names = ['exp'] + [r'%d' % i for i in range(1, max(alphas) + 1)]

    capsize_ = 2
    k = 2

    if region == "avg":
        colors = ["k" for _ in range(10)]

    elif region == "white":
        colors = np.array([[0.45176471, 0.76708958, 0.46120723, 1.        ],
                            [0.25259516, 0.66812764, 0.36286044, 1.        ],
                            [0.13402537, 0.54232987, 0.26828143, 1.        ],
                            [0.        , 0.42303729, 0.17071895, 1.        ],
                            [0.        , 0.26666667, 0.10588235, 1.        ]])

    elif region == "gray":
        colors = np.array([[0.41708574, 0.68063053, 0.83823145, 1.        ],
                            [0.25628604, 0.57001153, 0.7751634 , 1.        ],
                            [0.12710496, 0.44018454, 0.70749712, 1.        ],
                            [0.03137255, 0.31409458, 0.60648981, 1.        ],
                            [0.03137255, 0.18823529, 0.41960784, 1.        ]])

    datacolor = colors[0]

    colors = iter(colors)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if np.nan in alphas:
        label="data"
    else:
        label = None 

    ax.bar(x_e, y_e, yerr=st_e / np.sqrt(n_p), capsize=capsize_, color=datacolor, label=label, 
        # hatch='xx',
        hatch='...'
        )
    print("time  ", "   2h   ", "       6h     ", "        24h    ", "    48h   ")
    print("data   ", [format(x, ".4f") + "+-" + format(stdx, ".4f") for x, stdx in zip(y_e, st_e)]) #, st_e / np.sqrt(n_p))

    for i, alpha in enumerate(alphas):


        print("alpha", alpha,  [format(x, ".2f") + "+-" + format(stdx, ".2f") for x, stdx in zip(y_sim[i::(len(alphas))], st_sim[i::(len(alphas))])])
        
        if len(alphas) == 1:
            label = "simulation"
        
        if len(alphas) == 2 and alpha == 1:
            label = "diffusion"
        
        if np.isnan(alpha):
            hatch = "xx"
            label = "diffusion-reaction"
        else:
            hatch = None
        color = next(colors)

        ax.bar(x_sim[i::(len(alphas))], y_sim[i::(len(alphas))], yerr=st_sim[i::(len(alphas))] / np.sqrt(n_p),
                capsize=capsize_, label=label, color=color, hatch=hatch)

    ax.set_xticks([int(len(alphas)/2) + (len(alphas) + 2) * x for x in range(n_t)],)
    ax.set_xticklabels([r"$\sim 2\,$", r"$\sim 6\,$", r"$\sim 24\,$", r"$\sim 48\,$",], fontsize=fs)

    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylabel("tracer in " + region + " (mmol)", fontsize=fs)
    ax.set_xlabel("time (hours)", fontsize=fs)
    
    plt.ylim(-0.00, 0.25)
    
    if np.nan in alphas:
        plt.ylim(-0., 0.17)
    
    
    # ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15])
    # plt.locator_params(axis='y', nbins=6, tight=True)
    
    if np.nan in alphas:
        plt.legend(fontsize=fs)
    

    if paperformat:
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)