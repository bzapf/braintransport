import os
# from definitions import intervals
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from matplotlib.markers import TICKDOWN


print_messages = []

def get_data_in_intervals(pat, stored_times, stored_data, intervals):

    data = []
    ts = []

    assert max(stored_times) > 50 

    for idy, (tmin, tmax) in enumerate(intervals):
        for idx, (ti, datapoint) in enumerate(zip(stored_times, stored_data)):
            if tmin <= ti / 3600 < tmax and len(ts) <= idy:
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

    return ts, data





def significance_bar(start, end, height, fontsize, displaystring, text_dh, linewidth=1., markersize=5, color='k'):
    assert start != end
    # draw a line with downticks at the ends
    plt.plot([start, end], [height, height], '-', color=color, lw=linewidth,
             marker=TICKDOWN, markeredgewidth=linewidth, markersize=markersize)
    # draw the text with a bounding box covering up the line
    plt.text(0.5 * (start + end), text_dh + height, displaystring, ha='center', va='center',
             # bbox=dict(facecolor='1.', edgecolor='none', # boxstyle='Square,pad=' + str(boxpad)), 
             size=fontsize)



def load_experimental_data(pat, excelpath, roi):
        
    exceltable = pd.read_csv(os.path.join(excelpath, "experimental_data.csv"))
    
    assert roi in ["avg", "white", "gray", "avgds"]

    loaded_data = exceltable[roi]

    assert max(exceltable["t"]) < 2.5 * 24 * 3600

    times, data = get_data_in_intervals(pat, stored_times=exceltable["t"], stored_data=loaded_data)

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