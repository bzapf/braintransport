import os
import pylab
import matplotlib
import pandas
import seaborn as sns
import scipy
from collections import OrderedDict

#pylab.style.use('dark_background')
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["axes.linewidth"] = 3
matplotlib.rcParams["axes.labelsize"] = "xx-large"
matplotlib.rcParams["axes.titlesize"] = "xx-large"
matplotlib.rcParams["grid.linewidth"] = 1
matplotlib.rcParams["xtick.labelsize"] = "xx-large"
matplotlib.rcParams["ytick.labelsize"] = "xx-large"
matplotlib.rcParams["legend.fontsize"] = "xx-large"
matplotlib.rcParams["font.size"] = 14

def format_data(averages_6h, averages_24h):

    speeds = []
    sleep_deprived = set(["199", "227", "230", "235", "236", "241", "249"])

    for subject in averages_6h:
        (brain, gray, white, stem) = averages_6h[subject]

        dep = (subject in sleep_deprived)
        speeds += [(subject, brain, "Brain-wide", "6-24h", dep)]
        speeds += [(subject, gray, "Gray", "6-24h", dep)]
        speeds += [(subject, white, "White", "6-24h", dep)]
        speeds += [(subject, stem, "Stem", "6-24h", dep)]

    excludes = ["172", "191", "205"]
    for subject in averages_24h:
        if subject in excludes:
            continue
        (brain, gray, white, stem) = averages_24h[subject]
        dep = (subject in sleep_deprived)
        speeds += [(subject, brain, "Brain-wide", "24-48h", dep)]
        speeds += [(subject, gray, "Gray", "24-48h", dep)]
        speeds += [(subject, white, "White", "24-48h", dep)]
        speeds += [(subject, stem, "Stem", "24-48h", dep)]

    df = pandas.DataFrame(speeds)
    df.columns = ["Subject", "Speed", "Region", "Time", "Deprived"]

    # Rescale (from mm/h) to mum/min
    rescale = 1000/60 # mm/h to mum/min
    df["Speed"] = rescale*df["Speed"]

    #print(df[df.Region == "Brain-wide"].describe())
    
    return df

def is_there_a_difference(reg1, reg2, independent=True, equal_var=False):

    print("Mean values in A (%d) vs B (%d) are" % (reg1.count(), reg2.count()))
    print("... %.3g vs %.3g (mum/min)" % (reg1.mean(), reg2.mean()))
    reg1s = reg1.dropna()
    reg2s = reg2.dropna()

    if independent:
        a, pvalue = scipy.stats.ttest_ind(reg1s, reg2s, equal_var=equal_var)
    else:
        a, pvalue = scipy.stats.ttest_rel(reg1s, reg2s)

    if (pvalue < 0.05):
        print("Different (p-value: %.2g, a: %.3g) " % (pvalue, a))
    else:
        print("No significant difference (p-value: %.2g, a: %.3g)" % (pvalue, a))
    print("")

def print_table(phi_averages_6h, phi_averages_24h):
    
    df = format_data(phi_averages_6h, phi_averages_24h)

    for t in ("6-24h", "24-48h"):
        subjects = list(df[df.Region=="Brain-wide"]["Subject"][df.Time==t])
        brain_wide = list(df[df.Region=="Brain-wide"]["Speed"][df.Time==t])
        gray = list(df[df.Region=="Gray"]["Speed"][df.Time==t])
        white = list(df[df.Region=="White"]["Speed"][df.Time==t])
        stem = list(df[df.Region=="Stem"]["Speed"][df.Time==t])
        
        A = pandas.DataFrame({"ID": subjects, "v (mum/min)": brain_wide,
                              "v_g": gray, "v_w": white, "v_s": stem})
        print(A.to_latex(float_format="%.2f"))
    
def main(phi_averages_6h, phi_averages_24h):

    # Convert from whatever format you have to a pandas DataFrame
    # (makes pretty plotting easier)
    df = format_data(phi_averages_6h, phi_averages_24h)

    # Compare reference versus sleep-deprived cohort at the two time slots
    print("\n")
    print("*"*80)
    print("Comparing average regional phi's between the two cohorts at 6-24h, and 24-48h")
    for t in ["6-24h", "24-48h"]:
        print("-"*80)
        for r in ["Brain-wide", "Gray", "White", "Stem"]:
            print("%s (reference versus deprived) at %s" % (r, t))
            A = df["Speed"][df.Time == t][df.Region == r][df.Deprived==False]
            B = df["Speed"][df.Time == t][df.Region == r][df.Deprived==True]
            is_there_a_difference(A, B)

    # Compare regional differences for each cohort and time separately
    print("\n")
    print("*"*80)
    print("Comparing between regional phi's for each cohort at 6-24h and 24-48h")
    pairs = [("Gray", "White"), ("Gray", "Stem"), ("White", "Stem")]
    for t in ["6-24h", "24-48h"]:
        for (r0, r1) in pairs:
            print("-"*80)
            for c in [False, True]:
                label = "Sleep-deprived" if c else "Reference"
                print("%s cohort, %s versus %s at %s" % (label, r0, r1, t))
                A = df["Speed"][df.Time == t][df.Region == r0][df.Deprived==c]
                B = df["Speed"][df.Time == t][df.Region == r1][df.Deprived==c]
                is_there_a_difference(A, B, independent=False)

    # Compare regional differences for both cohorts combined at each tie
    print("\n")
    print("*"*80)
    t = "6-24h"
    print("Comparing between regional phi's for cohort combined at 6-24h and 24-48h")
    print("*"*80)

    print("At 6-24h (combined cohorts)")
    for r in ["Brain-wide", "Gray", "White", "Stem"]:
        A = df["Speed"][df.Time == t][df.Region == r]
        print(r, A.mean())

    for (r0, r1) in pairs:
        print("Combined, %s versus %s at %s" % (r0, r1, t))
        A = df["Speed"][df.Time == t][df.Region == r0]
        B = df["Speed"][df.Time == t][df.Region == r1]
        is_there_a_difference(A, B, independent=False)

    print("*"*80)

    t = "6-24h"
    print("Average flow speed statistics at %s (overall)" % t)
    print("-"*80)
    for r in ["Brain-wide", "Gray", "White", "Stem"]:
        A = df["Speed"][df.Time == t][df.Region == r]
        print(r)
        print(A.describe())

    t = "24-48h"
    print("-"*80)
    print("Average flow speed statistics at %s (overall)" % t)
    print("-"*80)
    for r in ["Brain-wide", "Gray", "White", "Stem"]:
        A = df["Speed"][df.Time == t][df.Region == r]
        print(r)
        print(A.describe())
    print("*"*80)

    pylab.figure(figsize=(8, 8))
    my_pal = {0: "dodgerblue", 1: "powderblue"}
    ax = sns.barplot(data=df[df.Time == "6-24h"], y="Speed", x="Region", hue="Deprived",
                     palette=my_pal, capsize=0.1, errwidth=3) 
    pylab.ylim(0, 6) 
    pylab.ylabel("Flow speed ($\mu$m/min)")
    pylab.xlabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if not os.path.isdir("results/graphics"):
        os.mkdir("results/graphics")
    
    handles, _ = ax.get_legend_handles_labels() 
    ax.legend(handles, ["Sleep", "Sleep-deprived"], loc="best")
    pylab.savefig("results/graphics/figure_ocd_average_6h.pdf")

    pylab.figure(figsize=(8, 8))
    my_pal = {0: "dodgerblue", 1: "powderblue"}
    ax = sns.barplot(data=df[df.Time == "24-48h"], y="Speed", x="Region", hue="Deprived",
                     palette=my_pal, capsize=0.1, errwidth=3) 
    pylab.ylim(0, 6) 
    pylab.ylabel("Flow speed ($\mu$m/min)")
    pylab.xlabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_legend().remove()
    pylab.savefig("results/graphics/figure_ocd_average_24h.pdf")

if __name__ == "__main__":

    from results.ocd_averages_6h import phi_averages_6h
    from results.ocd_averages_24h import phi_averages_24h

    main(phi_averages_6h, phi_averages_24h)

    print("\n"*5)
    print_table(phi_averages_6h, phi_averages_24h)

    pylab.show()
        
