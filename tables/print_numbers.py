import os, pathlib, json
import numpy as np
import pandas

import sys
import scipy
from pathlib import Path
here = Path(__file__).parent
import sys
sys.path.insert(0, str(here.parent))

from plots.definitions import groups, intervals
from plots.helpers import get_data_in_intervals, sim_at_mri_times
import argparse

pats = groups["all"]

header = ["(" + str(x[0]) + "," + str(x[1]) +")h" for x in intervals]

def read_data(average=False):

    sim_dataframes = {}
    mri_dataframes = {}

    for alpha in range(1, 6):
        alpha = str(alpha)

        for roi in ["avg", "gray", "white", "brainstem", "avgds"]:

            measured_c = {}
            simulated_c = {}

            for pat in pats:

                simulationcsv = pandas.read_csv(datafolder / pat / (resfolder + alpha) / "concs.csv")

                measuredcsv = pandas.read_csv(datafolder / pat / (resfolder + alpha) / "experimental_data.csv")
                

                mritimes, mridata = get_data_in_intervals(pat, stored_data=measuredcsv[roi], stored_times=measuredcsv["t"] / 3600, intervals=intervals)

                simt, simd = sim_at_mri_times(pat=pat, mritimes=mritimes, simtimes=simulationcsv["t"] / 3600, simdata=simulationcsv[roi])

                simtimes, simdata = get_data_in_intervals(pat, stored_times=simt, stored_data=simd, intervals=intervals)

                if roi == "avgds":
                    volumes = json.load(open(datafolder / pat / (resfolder + alpha) / "region_areas.json"))
                else:
                    volumes = json.load(open(datafolder / pat / (resfolder + alpha) / "region_volumes.json"))

                del simulationcsv, measuredcsv, simt, simd
                
                if not average:
                    simdata *= volumes[roi] / 1e6
                    mridata *= volumes[roi] / 1e6

                measured_c[pat] = mridata
                simulated_c[pat] = simdata
            
            measured_df = pandas.DataFrame.from_dict(measured_c, orient="index", columns=header)
            simulated_df = pandas.DataFrame.from_dict(simulated_c, orient="index", columns=header)

            sim_dataframes[(roi, alpha)] = simulated_df
            mri_dataframes[(roi, alpha)] = measured_df

    return sim_dataframes, mri_dataframes


def meanformat(values):
    return format(np.nanmean(values), mmol_digits) + " +- " + format(np.nanstd(values), mmol_digits)


def ratioformat(values1, values2):
    values = values1 / values2
    return format(np.nanmean(values), ratiodigits) + " +- " + format(np.nanstd(values), ratiodigits)


def percentformat(values):
    return format(np.nanmean(100 * values / 0.5), percent_digits) + " +- " + format(np.nanstd(100 * values / 0.5), percent_digits)


def maxformat(values):
    # print(dataframe)
    return format(np.nanmax(values), mmol_digits)

def minformat(values):
    # breakpoint()
    return format(np.nanmin(np.nanmax(values, axis=1)), mmol_digits)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    if "cluster" in os.getcwd():

        default_input_path = "/cluster/projects/nn9279k/bastian/SleepData/"

    else:
        default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"

    parser.add_argument("--inputpath", type=str, default=default_input_path)

    parserargs = parser.parse_args()
    argparse_dict = vars(parserargs)

    datafolder = argparse_dict["inputpath"]

    datafolder = pathlib.Path(datafolder)

    resfolder = "alphatests/alpha"

    sim_dataframes, mri_dataframes = read_data()

    sim_avg_dataframes, mri_avg_dataframes = read_data(average=True)

    #####################


    ratiodigits = ".1f"
    mmol_digits = ".3f"
    percent_digits = ".0f"

    ###############################################################################################################################

    print("-" * 80)
    print("Numbers for 'One-fourth of the tracers enter the brain'")

    regions = {"avg": "Brain-wide", "gray": "Cerebral cortex", "white": "Subcortical white matter"}

    roi = "avg"
    for time_idx in header[1:]:
        print(regions[roi], "tracer in interval", time_idx, meanformat(mri_dataframes[(roi, "1")][time_idx]), "mmol")
        print(regions[roi], "tracer in interval", time_idx, percentformat(mri_dataframes[(roi, "1")][time_idx]), "%")
        print()
        
        # Sanity check: total concentration should be roughly gray + white + brainstem
        total = mri_dataframes[("avg", "1")][time_idx] 
        subregion_sum = mri_dataframes[("gray", "1")][time_idx] + mri_dataframes[("brainstem", "1")][time_idx] + mri_dataframes[("white", "1")][time_idx]
        assert np.mean(np.abs(total- subregion_sum)) < 1e-6

    del time_idx

    roi = "white"
    print(regions[roi], "max tracer in (0, 48) h", maxformat(mri_dataframes[(roi, "1")]), "mmol")

    time_idx = header[-2]
    print()
    print(regions[roi], "group mean avg tracer at", time_idx, meanformat(mri_avg_dataframes[(roi, "1")][time_idx]), "mmol / L")
    print(regions[roi], "max peak avg tracer in (0, 48) h", maxformat(mri_avg_dataframes[(roi, "1")]), "mmol / L")
    print(regions[roi], "min peak avg tracer in (0, 48) h", minformat(mri_avg_dataframes[(roi, "1")]), "mmol / L")
    roi = "gray"
    print()
    print(regions[roi], "group mean avg tracer at", time_idx, meanformat(mri_avg_dataframes[(roi, "1")][time_idx]), "mmol / L")
    print(regions[roi], "max peak avg tracer in (0, 48) h", maxformat(mri_avg_dataframes[(roi, "1")]), "mmol / L")
    print(regions[roi], "min peak avg tracer in (0, 48) h", minformat(mri_avg_dataframes[(roi, "1")]), "mmol / L")

    del time_idx
    ###############################################################################################################################

    # After $\sim$6 hours, more tracer is observed clinically than diffusion simulations predict both in the cerebral cortex 
    # ($0.098 \pm 0.045$ vs $0.057 \pm 0.029$ mmol) and subcortical white matter ($0.016 \pm 0.007$ vs $0.002 \pm 0.002$ mmol)
    #  (Fig.~\ref{fig:fig2}E, F, $\alpha=1$). After 24 hours, observations and simulations of tracer amounts agree in the cerebral cortex
    #   ($0.095 \pm 0.038$ vs $0.105 \pm 0.039$ mmol) and subcortical white matter ($0.029 \pm 0.013$ vs $0.038 \pm 0.014$ mmol).
    #    But after 48 hours, clinically observed tracer amounts were smaller compared to simulations 
    #    ($0.046 \pm 0.025$ vs $0.072 \pm 0.035$ mmol in the cerebral cortex and $0.021 \pm 0.010$ vs $0.050 \pm 0.022$ mmol in the 
    #    subcortical white matter, Fig.~\ref{fig:fig2}D). Thus, extracellular diffusion alone also underestimates the tracer clearance.
    #     

    print("-" * 80)
    print("Numbers for 'Tracer influx and clearance is more rapid than by extracellular diffusion'")

    for time_idx in header[1:]:
        for roi in ["gray", "white"]:
            print(regions[roi], "measured tracer in interval", time_idx, meanformat(mri_dataframes[(roi, "1")][time_idx]), "mmol")
            print(regions[roi], "simulated tracer in interval", time_idx, meanformat(sim_dataframes[(roi, "1")][time_idx]), "mmol")
            print()


    ###############################################################################################################################

    print("-" * 80)
    print("Numbers for 'Enhanced diffusion predicts inaccurate influx and clearance interaction patterns'")

    alpha = "5"

    roi = "gray"
    time_idx = header[1]

    print(regions[roi], "measured tracer in interval", time_idx, meanformat(mri_dataframes[(roi, alpha)][time_idx]), "mmol")
    print(regions[roi], "simulated tracer in interval, alpha =", alpha, ":", time_idx, meanformat(sim_dataframes[(roi, alpha)][time_idx]), "mmol")

    print()

    for alpha in range(1,6):

        alpha = str(alpha)

        roi = "gray"
        time_idx = header[-1]

        print(regions[roi], "simulated tracer in interval, alpha =", alpha, ":", time_idx, meanformat(sim_dataframes[(roi, alpha)][time_idx]), "mmol")
        print(regions[roi], "ratio simulated / measured, alpha =", alpha, ":", time_idx, 
                    ratioformat(values1=sim_dataframes[(roi, alpha)][time_idx], values2=mri_dataframes[(roi, alpha)][time_idx]))


    # $\sim$6 hours ($0.098\pm 0.045$ vs $0.089 \pm 0.045$ mmol, $\alpha=5$, Fig.~\ref{fig:fig2}C). 
    # At 24 hours, the discrepancy increases with increasing $\alpha$. 
    # After 48 hours, the simulated values are essentially independent of $\alpha$ in the cerebral cortex (\mer{$0.07\pm 0.04$ for all $\alpha$}) 
    # and around \mer{two} times the tracer concentrations observed. 
    # 

    print() 

    for alpha in range(1,6):

        alpha = str(alpha)

        roi = "white"
        time_idx = header[1]

        print(regions[roi], "ratio simulated / measured, alpha =", alpha, ":", time_idx, 
                    ratioformat(values1=sim_dataframes[(roi, alpha)][time_idx], values2=mri_dataframes[(roi, alpha)][time_idx]))

    print()

    alpha = str(2)

    roi = "white"
    time_idx = header[2]

    print(regions[roi], "measured tracer in interval, alpha =", alpha, ":", time_idx, meanformat(mri_dataframes[(roi, alpha)][time_idx]), "mmol")
    print(regions[roi], "simulated tracer in interval, alpha =", alpha, ":", time_idx, meanformat(sim_dataframes[(roi, alpha)][time_idx]), "mmol")
    print(regions[roi], "ratio simulated / measured, alpha =", alpha, ":", time_idx, 
                ratioformat(values1=sim_dataframes[(roi, alpha)][time_idx], values2=mri_dataframes[(roi, alpha)][time_idx]))

    time_idx = header[2]
    print()
    print(regions[roi], "ratio simulated / measured, alpha =", alpha, ":", time_idx, 
                ratioformat(values1=sim_dataframes[(roi, alpha)][time_idx], values2=mri_dataframes[(roi, alpha)][time_idx]))




    # In the subcortical white matter, similar observations hold (Fig.~\ref{fig:fig2}D).
    #  After 6 hours, \mer{all} simulations underestimate the amount of tracer, while after 24 hours, simulations 
    #  with $\alpha = 2$ overestimate the data by a factor $2.1 \pm 0.6$ ($0.029\pm 0.013$ vs $0.056\pm 0.022$ mmol) with increasing discrepancy 
    #  for increasing $\alpha$. 
    #  Moreover, tracer concentrations in the subcortical white matter at 48 hours are overestimated by more than $3\times$. 

        
    print("Numbers for 'Tracer clearance by advection is reduced in sleep deprivation group.'")

    roi = "avg"
    time_idx = header[1]
    print(regions[roi], "Sleep group, measured tracer in interval", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]]), "mmol")
    print(regions[roi], "Sleep deprivation group, measured tracer in interval", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]]), "mmol")
    stat, pval = scipy.stats.ttest_ind(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]], 
                mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]], axis=0, equal_var=True, nan_policy="omit")
    print("Students's t-test, p=", pval)
    

    time_idx = header[-1]
    print(regions[roi], "Sleep group, measured tracer in interval", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]]), "mmol")
    print(regions[roi], "Sleep deprivation group, measured tracer in interval", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]]), "mmol")


    stat, pval = scipy.stats.ttest_ind(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]], 
                mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]], axis=0, equal_var=True, nan_policy="omit")
    print("Students's t-test, p=", pval)

    roi = "white"
    time_idx = header[-1]
    print(regions[roi], "Sleep group, measured tracer in interval,", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]]), "mmol")
    print(regions[roi], "Sleep deprivation group, measured tracer in interval", ":", time_idx, 
        meanformat(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]]), "mmol")

    stat, pval = scipy.stats.ttest_ind(mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]], 
                mri_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]], axis=0, equal_var=True, nan_policy="omit")
    print("Students's t-test, p=", pval)

    print("Surface differences")
    for time_idx in header:
        roi = "avgds"
        stat, pval = scipy.stats.ttest_ind(mri_avg_dataframes[(roi, alpha)][time_idx].loc[groups["sleep"]], 
                    mri_avg_dataframes[(roi, alpha)][time_idx].loc[groups["sleepdep"]], axis=0, equal_var=True, nan_policy="omit")
        print("Students's t-test, p=", pval)