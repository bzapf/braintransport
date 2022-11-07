import numpy as np
import os
import pandas as pd
import pickle
from definitions import intervals, reaction_resultfolder
import matplotlib.pyplot as plt
import json




def make_barplot(region, pats, alphas, paperformat, resultfoldername, data_folder, savepath, fs, figsize, dpi):
    n_t = 4  # number of time points
    n_a = len(alphas)  # number of alphas
    n_p = len(pats)  # number of patients
    conc_simulated = np.zeros((n_t, n_p, n_a)) + np.nan
    conc_experimental = np.zeros((n_t, n_p))  + np.nan

    for alpha_idx, alpha in enumerate(alphas):
        
        experimental_conc_sum = np.zeros(n_t)
        ne_sum = np.zeros(n_t)
        s_sum = np.zeros(n_t) 
        pat_no = 0

        for pat in pats:

            resolution = 32
        
            if "alphatest" in resultfoldername(pat) and not np.isnan(alpha):
                folder = data_folder + str(pat) + "/" + resultfoldername(pat) + "/alpha" + str(alpha) + "/"
            
            elif np.isnan(alpha):
                folder = reaction_resultfolder(pat, best=True)
                assert folder is not None
                assert os.path.isdir(folder)
            else:

                raise ValueError
                sim_folder = 'D%d.0RD0.00e+00RN0.00e+00' % alpha
                folder = data_folder + str(pat) + resultfoldername(pat) + '%s/' % (sim_folder)

                try:
                    dts = os.listdir(folder)
                except:
                    breakpoint()
                dt = sorted([x for x in dts], key=lambda x: int(x[3:]), )[0]
                dt = int(dt[3:])
                folder += "%d_%d/" % (resolution, dt)

                assert os.path.isdir(folder)
            
            try:
                
                params = json.load(open(folder + "params"))

                if not params["concfolder"] == "FIGURES_CONC":
                    breakpoint()

                experimental = pd.read_csv(folder + 'experimental_data.csv')
                concs = pd.read_csv(folder + 'concs.csv')
            except pd.errors.EmptyDataError:
                print("concs.csv Not ready,", folder, "CONTINUE")

            with open(data_folder + str(pat) + '/region_volumes' + str(resolution), 'rb') as f:
                region_volumes = pickle.load(f)
                roi_volume = region_volumes[region] / 1e6

            t0 = concs['t'] / 3600
            c0 = concs[region]

            assert max(np.abs(c0)) < 1
            
            experimental_times = experimental['t'] / 3600
            c1 = experimental[region]

            for i, interval in enumerate(intervals):
                found = False
                for exp_idx, t_exp in enumerate(experimental_times):
                    if not (t_exp >= interval[0] and t_exp <= interval[1]):
                        continue
                    found = True
                    break

                # if not found and pat == 241:
                #     addval = np.nan
                #     found = True

                # else:
                #     ne_sum[i] += 1

                if not found:
                    addval = np.nan
                    # print(pat, interval, "not found, appending nan")
               
                conc_idx = np.argmin(np.abs(t0 - t_exp))

                addval = float(c1[exp_idx] * roi_volume)
                conc_experimental[i, pat_no] = addval
                experimental_conc_sum[i] += addval

                addval_simulated = roi_volume * c0[conc_idx]

                conc_simulated[i, pat_no, alpha_idx + 1 - 1] = addval_simulated   
                s_sum[i] += addval_simulated

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
        # k += j
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
        # return
        colors = ["k" for _ in range(10)]
    else:
        colors = np.load("/home/basti/Pictures/sleep_paper/fig2/colors_" + region + ".npy")
        # breakpoint()

    datacolor = colors[0]

    colors = iter(colors)

    USE_VEGARDS_FORMAT = False

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    if np.nan in alphas:
        label="data"
    else:
        label = None 

    ax.bar(x_e, y_e, yerr=st_e / np.sqrt(n_p), capsize=capsize_, color=datacolor, 
        label=label, 
        # hatch='xx',
        hatch='...'
        )
    print("time  ", " 2h ", "     6h   ", "      24h  ", "  48h   ")
    print("data  ", [format(x, ".2f") for x in y_e]) #, st_e / np.sqrt(n_p))

    for i, alpha in enumerate(alphas):
        # label = r"$\alpha=$"+format(alphas[i], ".0f")
        
        # breakpoint()
        print("alpha", alpha, [format(x, ".2f") for x in y_sim[i::(len(alphas))]],
                # [format(x, ".2f") for x in st_sim[i::(len(alphas))  / np.sqrt(n_p)]]
                ) # , st_sim[i::(len(alphas))] / np.sqrt(n_p))
        
        if len(alphas) == 1:
            label = "simulation"
        
        if len(alphas) == 2 and alpha == 1:
            label = "diffusion"
        
        

        if np.isnan(alpha):
            hatch = "xx"
            # color = "goldenrod"
            label = "diffusion-reaction"
        else:
            hatch = None
        color = next(colors)

        ax.bar(x_sim[i::(len(alphas))], y_sim[i::(len(alphas))], yerr=st_sim[i::(len(alphas))] / np.sqrt(n_p),
                capsize=capsize_, label=label, color=color, hatch=hatch)




    ax.set_xticks([int(len(alphas)/2) + (len(alphas) + 2) * x for x in range(n_t)],)
    ax.set_xticklabels([r"$\sim 2\,$", r"$\sim 6\,$", r"$\sim 24\,$", r"$\sim 48\,$",], fontsize=fs)
    # breakpoint()
    ax.tick_params(axis='x', length=0)

    ax.tick_params(axis='y', labelsize=fs)
    ax.set_ylabel("average tracer in " + region + " (mmol / L)", fontsize=fs)
    ax.set_xlabel("time (hours)", fontsize=fs)
    
    plt.ylim(-0.00, 0.15)
    
    if np.nan in alphas:
        plt.ylim(-0., 0.17)
    
    
    ax.set_yticks([0., 0.03, 0.06, 0.09, 0.12, 0.15])
    # plt.locator_params(axis='y', nbins=6, tight=True)
    if np.nan in alphas:
        plt.legend(fontsize=fs)
    
    # plt.tight_layout()


    if paperformat:
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)