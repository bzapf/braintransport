import os
import pylab
import matplotlib
import pandas
import seaborn as sns
import scipy
from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter

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

def format_data(div_phi_averages_6h, div_phi_averages_24h):

    divs = []
    sleep_deprived = set(["199", "227", "230", "235", "236", "241", "249"])

    for subject in div_phi_averages_6h:
        (brain, gray, white, stem) = div_phi_averages_6h[subject]

        dep = (subject in sleep_deprived)
        divs += [(subject, brain, "Brain-wide", "6-24h", dep)]
        divs += [(subject, gray, "Gray", "6-24h", dep)]
        divs += [(subject, white, "White", "6-24h", dep)]
        divs += [(subject, stem, "Stem", "6-24h", dep)]

    excludes = ["172", "191", "205"]
    for subject in div_phi_averages_24h:
        if subject in excludes:
            continue
        (brain, gray, white, stem) = div_phi_averages_24h[subject]
        dep = (subject in sleep_deprived)
        divs += [(subject, brain, "Brain-wide", "24-48h", dep)]
        divs += [(subject, gray, "Gray", "24-48h", dep)]
        divs += [(subject, white, "White", "24-48h", dep)]
        divs += [(subject, stem, "Stem", "24-48h", dep)]

    df = pandas.DataFrame(divs)
    df.columns = ["Subject", "Div(phi)", "Region", "Time", "Deprived"]

    # Rescale (from 1/h) to 1/min and then by 10000 to make figure more pretty
    rescale = 1.e4*1/60 # mm/h to mum/min

    df["Div(phi)"] = rescale*df["Div(phi)"]

    return df


def main(div_phi_averages_6h, div_phi_averages_24h):

    df = format_data(div_phi_averages_6h, div_phi_averages_24h)

    pylab.figure(figsize=(8, 8))
    my_pal = {0: "hotpink", 1: "pink"}
    ax = sns.boxplot(data=df[df.Region=="Brain-wide"],
                     y="Div(phi)", x="Time", hue="Deprived",
                     palette=my_pal, width=0.9)
    
    pylab.ylim(-1.5, 6) 
    ax.set_yticks([-1, 0, 3, 6])
    pylab.ylabel("Fluid influx rate ($\\times 10^{-4}$ min$^{-1}$)")
    pylab.xlabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pylab.subplots_adjust(left=0.2)
    handles, _ = ax.get_legend_handles_labels() 
    ax.legend(handles, ["Sleep", "Sleep-deprived"], loc="lower right")

    if not os.path.isdir("results/graphics"):
        os.mkdir("results/graphics")
    
    pylab.savefig("results/graphics/figure_div_phi.pdf")
    
    print("Comparing fluid influx rates between cohorts at 6-24h")
    A0 = df[df.Time=="6-24h"][df.Deprived==False][df.Region=="Brain-wide"]["Div(phi)"]
    A0.dropna()
    print(A0.describe())
    A1 = df[df.Time=="6-24h"][df.Deprived==True][df.Region=="Brain-wide"]["Div(phi)"]
    A1.dropna()
    print(A1.describe())

    a, pvalue = scipy.stats.ttest_ind(A0, A1, equal_var=False)
    if (pvalue < 0.05):
        print("Different (p-value: %.2g, a: %.3g) " % (pvalue, a))
    else:
        print("No significant difference (p-value: %.2g, a: %.3g)" % (pvalue, a))
    print("")

    print("Comparing fluid influx rates between cohorts at 24-48h")
    A0 = df[df.Time=="24-48h"][df.Deprived==False][df.Region=="Brain-wide"]["Div(phi)"]
    A0.dropna()
    print(A0.describe())
    A1 = df[df.Time=="24-48h"][df.Deprived==True][df.Region=="Brain-wide"]["Div(phi)"]
    A1.dropna()
    print(A1.describe())

    a, pvalue = scipy.stats.ttest_ind(A0, A1, equal_var=False)
    if (pvalue < 0.05):
        print("Different (p-value: %.2g, a: %.3g) " % (pvalue, a))
    else:
        print("No significant difference (p-value: %.2g, a: %.3g)" % (pvalue, a))
    print("")

    

if __name__ == "__main__":

    from results.div_phi_averages_6h import div_phi_averages_6h
    from results.div_phi_averages_24h import div_phi_averages_24h

    main(div_phi_averages_6h, div_phi_averages_24h)

    pylab.show()
        

