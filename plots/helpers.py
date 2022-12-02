import os
# from definitions import intervals
import numpy as np
from matplotlib.markers import TICKDOWN
from typing import Callable, Union
import json, pathlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate



print_messages = []

def get_data_in_intervals(pat, stored_times, stored_data, intervals):

    data = []
    ts = []

    assert len(stored_data) < 12, "this function should only be called with experimental data"
    assert len(stored_times) < 12, "this function should only be called with experimental data"

    if max(stored_times) < 80:
        message = "Assuming times are in hours"
        if message not in print_messages:
            print(message)

        print_messages.append(message)
        
        stored_times = stored_times * 3600


    for idy, (tmin, tmax) in enumerate(intervals):
        for idx, (ti, datapoint) in enumerate(zip(stored_times, stored_data)):
            if tmin <= ti / 3600 < tmax and len(ts) <= idy:

                # print(pat, format(ti / 3600, ".1f"), "(tmin, tmax) =", tmin, tmax)
                ts.append(ti / 3600)
                data.append(datapoint)
                break
        
        if not len(ts) == (idy + 1):
            ts.append(np.nan)
            data.append(np.nan)

            print_message = str("No data available in interval" + str(intervals[idy]) + " for " + str(pat) + ", appending nan")

            if print_message not in print_messages:
                print(print_message)

                print_messages.append(print_message)

        # print(len(print_messages))

    assert len(data) == len(intervals)
    assert len(ts) == len(intervals)

    return np.array(ts), np.array(data)



def sim_at_mri_times(pat, mritimes: Union[np.ndarray, list], simtimes: Union[np.ndarray, list], simdata: Union[np.ndarray, list]):

    data = []
    ts = []

    # breakpoint()

    if isinstance(simtimes, list):
        simtimes = np.array(simtimes).astype(float)

    assert np.nanmax(mritimes) < 80
    assert np.max(simtimes) < 80, "different time unit for sim and mritime?"

    for mrit in mritimes:

        if np.isnan(mrit):
            print_message = str("Found nan in mritimes for pat " + pat + ", appending nan for simdata")

            if print_message not in print_messages:
                print(print_message)
                print_messages.append(print_message)

            ts.append(np.nan)
            data.append(np.nan)        

            continue

        idx = np.argmin(np.abs(mrit - simtimes))

        # print("mrit=", mrit, "appending", simtimes[idx], ".....", list(simtimes[idx-2:idx+2]))

        ts.append(simtimes[idx])
        data.append(simdata[idx])

    # breakpoint()
    assert len(data) == len(mritimes)

    return np.array(ts), np.array(data)



def extract_data(alphas, pats, data_folder, resultfoldername: Callable, region, average_tracer, intervals):

    n_t = len(intervals)  # number of time points
    n_a = len(alphas)  # number of alphas
    n_p = len(pats)  # number of patients

    assert callable(resultfoldername)
    
    #  TODO FIXME

            # if np.nan not in alphas:
            #     folder = data_folder + str(pat) + "/" + resultfoldername + "/alpha" + str(alpha) + "/"
            
            # else:
            #     # assert np.isnan(alpha)
            #     # folder = reaction_resultfolder(pat, best=True)
            #     # assert folder is not None
            #     # assert os.path.isdir(folder)

            #     folder = reaction_resultfolder(pat, best=True)



    conc_simulated = np.zeros((n_t, n_p, n_a)) + np.inf # -42 # + np.nan
    conc_experimental = np.zeros((n_t, n_p))  + np.inf 

    for alpha_idx, alpha in enumerate(alphas):
        
        for pat_no, pat in enumerate(pats):

            folder = pathlib.Path(data_folder) / resultfoldername(pat, alpha)

            params = json.load(open(folder / "params"))

            assert params["concfolder"] == "FIGURES_CONC"

            if region == "avgds":
                with open(folder / 'region_areas.json') as f:
                    region_volumes = json.load(f)
                    roi_volume = region_volumes[region] / 1e6
            else:                
                with open(folder / 'region_volumes.json') as f:
                    region_volumes = json.load(f)
                    roi_volume = region_volumes[region] / 1e6

            # if alpha == 1:
            #     print(pat, region_volumes[region])

            if average_tracer:
                roi_volume = 1 # 1e-6

            experimental = pd.read_csv(folder / 'experimental_data.csv')
            
            if np.isnan(alpha) or "alphatests" in str(folder):
                concs = pd.read_csv(folder / 'concs.csv')
            else:
            # if np.nan in alphas and (resultfoldername != "alphatests"):
                concs = pd.read_csv(folder / 'concs_plain.csv') 
            
            simulation_times = concs['t'] / 3600
            simulated_tracer = concs[region] * roi_volume

            assert max(np.abs(simulated_tracer)) < 1
            
            experimental_times = experimental['t'] / 3600
            measured_tracer = experimental[region] * roi_volume

            experimental_times , measured_tracer_at_times = get_data_in_intervals(pat, stored_times=experimental_times, stored_data=measured_tracer, intervals=intervals)

            # print(experimental_times, simulation_times)
            simulation_times, simulated_tracer_at_times = sim_at_mri_times(pat, mritimes=experimental_times, simdata=simulated_tracer, simtimes=simulation_times)    

            # simulation_times, simulated_tracer_at_times = get_data_in_intervals(pat, stored_times=simulation_times, stored_data=simulated_tracer, intervals=intervals)

            # del simulated_tracer

            # if pat == "241":
            #     breakpoint()
            conc_experimental[:, pat_no] = measured_tracer_at_times
            conc_simulated[:, pat_no, alpha_idx] = np.where(np.isnan(measured_tracer_at_times), np.nan, simulated_tracer_at_times)

    assert np.nanmax(conc_experimental) < np.inf
    assert np.nanmax(conc_simulated) < np.inf

    return conc_experimental, conc_simulated


def significance_bar(start, end, height, fontsize, displaystring, text_dh, linewidth=1., markersize=5, color='k'):
    assert start != end
    # draw a line with downticks at the ends
    plt.plot([start, end], [height, height], '-', color=color, lw=linewidth,
             marker=TICKDOWN, markeredgewidth=linewidth, markersize=markersize)
    # draw the text with a bounding box covering up the line
    plt.text(0.5 * (start + end), text_dh + height, displaystring, ha='center', va='center',
             # bbox=dict(facecolor='1.', edgecolor='none', # boxstyle='Square,pad=' + str(boxpad)), 
             size=fontsize)



def load_experimental_data(pat, excelpath, roi, intervals):
        
    exceltable = pd.read_csv(os.path.join(excelpath, "experimental_data.csv"))
    
    assert roi in ["avg", "white", "gray", "avgds"]

    loaded_data = exceltable[roi]

    assert max(exceltable["t"]) < 2.5 * 24 * 3600

    times, data = get_data_in_intervals(pat, simtimes=exceltable["t"], simdata=loaded_data, intervals=intervals)

    return times, data




def load_simulation_data(pat, excelpath, roi):

    exceltable = pd.read_csv(excelpath)

    t = np.array(exceltable["t"].tolist()) / 3600
    times = np.linspace(np.min(t), np.max(t), 100)

    assert roi in ["avg", "white", "gray"]

    stored_data = exceltable[roi]

    f = interpolate.interp1d(t, stored_data)

    data = f(times)

    return times.tolist(), data