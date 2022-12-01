import os
import pylab
import matplotlib
import pandas
import seaborn as sns
import scipy
from collections import OrderedDict
import matplotlib.pyplot as plt
import scipy


#pylab.style.use('dark_background')
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["axes.linewidth"] = 3
matplotlib.rcParams["axes.labelsize"] = "xx-large"
matplotlib.rcParams["axes.titlesize"] = "xx-large"
matplotlib.rcParams["grid.linewidth"] = 1
matplotlib.rcParams["xtick.labelsize"] = "xx-large"
matplotlib.rcParams["ytick.labelsize"] = "xx-large"
matplotlib.rcParams["legend.fontsize"] = "x-large"
matplotlib.rcParams["font.size"] = 12

def format_data(averages_6h, averages_24h):

    speeds = []
    sleep_deprived = set(["199", "227", "230", "235", "236", "241", "249"])

    for subject in averages_6h:
        (brain, gray, white, stem) = averages_6h[subject]

        dep = (subject in sleep_deprived)
        speeds += [(subject, brain, "Brain-wide", "6-24h", dep)]
        speeds += [(subject, gray, "Cerebral cortex", "6-24h", dep)]
        speeds += [(subject, white, "Subcortical white matter", "6-24h", dep)]
        speeds += [(subject, stem, "Brain stem", "6-24h", dep)]

    excludes = ["172", "191", "205"]
    for subject in averages_24h:
        if subject in excludes:
            continue
        (brain, gray, white, stem) = averages_24h[subject]
        dep = (subject in sleep_deprived)
        speeds += [(subject, brain, "Brain-wide", "24-48h", dep)]
        speeds += [(subject, gray, "Cerebral cortex", "24-48h", dep)]
        speeds += [(subject, white, "Subcortical white matter", "24-48h", dep)]
        speeds += [(subject, stem, "Brain stem", "24-48h", dep)]

    df = pandas.DataFrame(speeds)
    df.columns = ["Subject", "Speed", "Region", "Time", "Deprived"]

    # Rescale (from mm/h) to mum/min
    rescale = 1000/60 # mm/h to mum/min
    df["Speed"] = rescale*df["Speed"]

    #print(df[df.Region == "Brain-wide"].describe())
    
    return df



if __name__ == "__main__":

# To generate these tables, run these commands first at Saga in the current directory
# salloc --ntasks=1 --mem-per-cpu=8G --time=00:59:00 --account=NN9279K --qos=devel
# source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

    from results.ocd_averages_6h import phi_averages_6h
    from results.ocd_averages_24h import phi_averages_24h

    df = format_data(phi_averages_6h, phi_averages_24h)

    df.to_csv("./results/dataframe.csv")


    df = pandas.read_csv("./results/dataframe.csv")

    # print(df)

    times = ["6-24h", "24-48h"]

    for t in times:
        print("-"*80)
        for r in ["Brain-wide", "Cerebral cortex", "Subcortical white matter", "Brain stem"]:
            print("%s (reference versus deprived) at %s" % (r, t))
            A = df["Speed"][df.Time == t][df.Region == r][df.Deprived==False]
            B = df["Speed"][df.Time == t][df.Region == r][df.Deprived==True]

            print(len(A), len(B))

    # exit()

    patsa = set(df["Subject"][df.Time == times[0]])
    
    patsb = set(df["Subject"][df.Time == times[1]])

    allpats = set.intersection(patsa, patsb)

    # allpats.remove(228)

    df = df.loc[df['Subject'].isin(allpats)]

    # print(df)
    # exit()


    for r in ["Brain-wide", "Cerebral cortex", "Subcortical white matter", "Brain stem"]:
        A = df["Speed"][df.Time == times[0]][df.Region == r]
        B = df["Speed"][df.Time == times[1]][df.Region == r]
        
        # if len(A) == 0 or len(B) == 0:
        #     print(r)
        #     exit()

        test = scipy.stats.pearsonr(A, B)
        
        print("tid pearson r:", format(test[0], ".2f"))

        plt.figure()
        # label = r + ", r="+ format(test[0], ".2f")
        plt.title(r + ", r = "+ format(test[0], ".2f"))
        plt.plot(A,B, marker="o", linewidth=0, color="k")
        plt.xlabel("$|\overline{v}|$" + " (6 -> 24 h)")
        plt.ylabel("$|\overline{v}|$" + " (24 -> 48 h)")

        plt.tight_layout()

        plt.savefig("/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/Brain influx and clearance during sleep and sleep deprivation/figures/simulations/"
        + r + "-correlation.png", dpi=420)

        # plt.show()
