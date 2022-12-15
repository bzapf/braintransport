import argparse
import os
import json
import numpy as np
import itertools

import pandas
import scipy
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from scipy import stats

import pathlib
import warnings

import scipy

from barplots import make_barplot, make_figs
from definitions import pats, groups, datafolder, ylabels, reaction_resultfolder




def r2halflife(r):
    """
    """

    if not np.max(r) < 0.001:
        r = r / (1e4 * 60)


    retval = np.log(2) / (r)
    
    return retval


def stderr(x, axis): 
    return stats.sem(x, axis=axis, nan_policy="omit")

def barplot(qty, qtyname, dpi, figsize, savepath=None, twinax=True):

    if isinstance(qty, list):
        plt.close("all")

    plt.figure(figsize=figsize, dpi=dpi)
 
    global unit
    # unit, scale = "days", 24 * 3600
    unit, scale = "$10^{-4}$ min", 1e4 * 60

    halflifeunit, halflifescale = "h", 1 / (3600) # seconds to hours
    # halflifeunit, halflifescale = "min", 1 / (60) # seconds to min

    if "alpha" in qty:
        unit, scale = "", 1
        halflifeunit, halflifescale = "", 1 # seconds to hours

    my_pal = {0: "hotpink", 1: "lightpink"}

    if "alpha" in qtyname:
        colors = sns.color_palette("Paired")
        purple1 = colors[9]
        purple2 = colors[8]

        my_pal = {0: purple1, 1: purple2}

   
    if max(dfs[qty]) > 0.001:
        scale = 1

    ax = sns.boxplot(data=[dfs[qty] * 1 * scale, dfsdep[qty] * scale], 
            # data=dfs, x=df["sleep"], y=qty,
            palette=my_pal# palette="Set1"
    )

    
    ax.set_ylabel("Reaction rate " + r"$" + qtyname + r"$" + " (" + unit + r"$^{-1}$)")


    if not "alpha" in qty:
        ax.set_ylim(4, 62)

    if "alpha" in qty:
        ax.set_ylabel("best dispersion factor " + r"$\alpha$")
        plt.locator_params(axis='y', nbins=5)

    plt.xticks([0, 1], ["Sleep", "Sleep deprivation"])

    if twinax:
        ax_c = ax.twinx()        
        ax_c.set_ylabel('half life $T_{1/2}=$ln$2/' + qtyname + '$ (' + halflifeunit + ')', #fontsize=fs
        )
        ax_c.tick_params(axis='y', #labelsize=fs
        )

    ax.set_yticks([x for x in ax.get_yticks() if x > 0])

    if twinax:
        y1ticks = ax.get_yticks()
        ax_c.set_yticks(y1ticks, [format(r2halflife(x / scale) * halflifescale, ".1f") for x in y1ticks])
    
    if not twinax:
        sns.despine()
    
    plt.tight_layout()


    if savepath is not None:
        plt.savefig(savepath, dpi=dpi)




def plot_all(region, pats, path_to_files, fs, figsize, savepath=None, legend_outside=False, ylabel="avg. tracer in brain (mmol / L)", average = False,
        savedpi=300, legendfs=12):

    raise NotImplementedError("Deprecated method")

    n = 100

    parameters = ["plain", "best"]

    labels=[r"$\alpha=1, r=0, \phi=0$",      
            r"best $(\alpha, r)$"]

    if not average:
        assert "avg" not in ylabel
    else:
        assert "avg" in ylabel

    days = 2

    x = np.linspace(0, days, n)

    y_sim = np.zeros((n, len(parameters), len(pats))) + np.nan
    y_data = np.zeros((len(intervals), len(pats))) + np.nan
    datatimes = np.zeros((len(intervals), len(pats))) + np.nan

    n_t = len(intervals)
    n_a = len(parameters)  # number of alphas
    n_p = len(pats)  # number of patients


    colors = ["tab:orange", "tab:green"]

    for pat_counter, pat in enumerate(pats):

        # os.chdir("/")

        paths_to_plot, datatable = get_paths(pat)

        if paths_to_plot is None:
            print("paths is None for", pat, "continue")
            continue

        if pat_counter == 0 and "opt_r_tests_const" in paths_to_plot[1]:
            labels=[r"$\alpha=1, r=0, \phi=0$",      
                r"best $(\alpha, r)$"]

        try:
            for f in paths_to_plot:
                assert os.path.isfile(f)
            assert os.path.isfile(datatable)
        except AssertionError:
            breakpoint()
            continue

        expdata = pandas.read_csv(datatable)

        for param_counter, excelpath in enumerate(paths_to_plot):

            if param_counter == 0:

                times, data = get_data_in_intervals(pat, stored_times=expdata["t"], stored_data=expdata[region])

                if not average:
                    with open(os.path.join(datafolder, pat, 'region_volumes' + str(32)), 'rb') as f:
                        region_volumes = pickle.load(f)
                        assert "ds" not in region
                    
                    data = np.array(data)
                    data *= region_volumes[region] * 1e-6



                y_data[:, pat_counter] = data
                datatimes[:, pat_counter] = days * np.array(times) / max(times)
                # datatimes[:, pat_counter] = np.array(times) / 2 * days

            simtimes, simdata, _ = extract_data(pat, excelpath, average=average, load_simulation_results=True, 
                                    path_to_files=path_to_files, resolution="32", 
                                    load_experimental_results=False, regions=[region])
       

            y_sim[:, param_counter, pat_counter] = simdata

            assert len(simdata.shape) == 1
            assert simdata.shape[0] == y_sim.shape[0]
            del simtimes, simdata

    fig = plt.figure(figsize=figsize, dpi=dpi)
    
    xd = np.nanmean(datatimes, axis=1)
    yd = np.nanmean(y_data, axis=1)
    yerr = stderr(y_data, axis=1)
    xerr = stderr(datatimes, axis=1)
    
    plt.errorbar(xd, yd, yerr=yerr, xerr=xerr,
        elinewidth=1, capsize=1,
        marker="o", linewidth=0, markersize=12, color="k")

    # # breakpoint()
    # try:
    #     assert y_sim.min() >= - 1e12
    #     assert np.sum(np.isnan(y_sim)) == 0
    # except AssertionError:
    #     breakpoint()

    y = np.nanmean(y_sim, axis=2)
    yerr = stderr(y_sim, axis=2)

    for idx in range(len(parameters)):

        _label = ""

        lw = 2

        if labels is not None:
            _label = labels[idx]
              
        plt.fill_between(x=x, y1=(y - yerr)[:, idx], y2=(y + yerr)[:, idx], alpha=0.5, color=colors[idx])
        plt.plot(x, y[:, idx], linewidth=lw, label=_label, color=colors[idx])

    plot_night(fs=fs, max_t=50 / 2)

    if legend_outside is True:
        plt.legend(loc="upper left", fontsize=legendfs, bbox_to_anchor=(1.02, 0.8))
    elif legend_outside == "no":
        pass  
    else:
    
        plt.legend(loc="lower right", fontsize=fs, bbox_to_anchor=(0.5, 0.1, 0.5, 0.5))

    # plt.xticks(fontsize=fs)
    # plt.xlabel("time (days)", fontsize=fs)
    
    plt.xticks([6 / 24, 1, 2], labels=[r"$6$", r"$24$", r"$48$"], fontsize=fs)
    plt.xlabel("time (hours)", fontsize=fs)

    plt.yticks(fontsize=fs)
    plt.locator_params(axis='x', nbins=3)
    plt.locator_params(axis='y', nbins=5)
    plt.ylabel(ylabel=ylabel, fontsize=fs)
    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=savedpi)


    




def make_df():

    result_dict = {}

    header = ["sleep", "alpha", "r_d", # "r_n", "rd/rn", 
            "j_plain", # "j_init", 
            "jfinal", "red. in J (%)", "time steps", 
            # "lbfgs-iters", 
            # "compute time / iter"
            ]
    dones = []

    for pat in pats:

        subfolder = reaction_resultfolder(pat, best=True)

        if subfolder is None:
            result_dict[pat] = [np.nan for x in header]
            print("Empty subfolder for,", pat, " CONTINUE")
            continue
        
        try:
            params = json.load(open(subfolder + "/params"))
        except json.decoder.JSONDecodeError:
            print("vim " + subfolder + "/params")
        # jpath = subfolder + "/J_during_optim.txt"
        # j_during_train = np.genfromtxt(jpath, delimiter=',')

        # if np.nan in j_during_train:
        #     breakpoint()
        
        # j = j_during_train[0]

        try:
            params["j_d_plain"]
            params2 = params

        except KeyError:
            subfolder2 = "./" + pat + "/" + "opt_r_tests_const" + "/iter1/k96"
            breakpoint()
            # TODO FIXME
            # this should not be needed when the results for k144 are done for all
            params2 = json.load(open(subfolder2 + "/params"))

        # print(sorted(params.keys()))
        
        red = abs(params2["j_d_plain"] - params["j_d_final"]) / params2["j_d_plain"]

        # if red < 0.1:
        #     print("Improvement < 0.1, CONTINUE")
        #     continue


        sleepq = pat in sleepdep

        row = [not sleepq, params["alpha_final"], 
                1e4* params["r_d_final"] * 60, # 1/s to 10^-4 / min 
        # params["r_n_final"]
        ]
        # row.append(params["r_d_final"] / params["r_n_final"])
        row.append(params2["j_d_plain"])
        # row.append(j)
        row.append(params["j_d_final"])
        row.append(int(red * 100))
        row.append(int(int(pathlib.Path(subfolder).stem[1:])))
        # row.append(int((j_during_train.size-1)/12))
        # row.append(format(params["optimization_time_hours"]/ int((j_during_train.size-1)/12), ".2f"))



        result_dict[pat] = row

        dones.append(pat)

    print("Including", len(dones), "patients:")
    print(dones)

    df = pandas.DataFrame.from_dict(result_dict, orient="index", columns=header)# .transpose()

    return df, dones


iterks = ["96", "144", "288", "576"]

dts = [format(60 * 48 / (float(x)), ".0f") for x in iterks]

dt_lookup = {}
for i, dt in zip(iterks, dts):
    dt_lookup[i] = dt




def resolutiontable2(pats, latexname=None):

    

    def formatfun(x, qty):

        if "r_d_final" in qty and x < 1:
            x = 1e4 * 60 * x
            return format(x, ".0f") 

        elif "j_d" in qty or "j_d_init" in qty:
            return format(x, ".2f")
        else:
            return format(x, ".1f")

    columns = ["alpha_final", "r_d_final", "j_d_plain","j_d_init", "j_d_final"]
    

    # table_header = ["pat", "dt (min)", r"$\alpha$", r"$r$ ($10^{-4}\,$min$^{-1}$) ", r"$J(\alpha=1, r=0)$", r"$J_0(\alpha=3, r=10^{-5}\,$s$^{-1}$)$", r"$J_F$"]
    table_header = ["pat", "dt (min)", r"$\alpha$", r"$r$ ($10^{-4}\,$min$^{-1}$) ", r"$J_p$", r"$J_0$", r"$J_F$"]


    all_rows = []

    for pat in pats:

        lowest = 1e16
        lowest_counter = None
        patrows = []
        
        alphas, rs = [], []

        for counter, iterk in enumerate(iterks):

            subfolder = reaction_resultfolder(pat, k="k"+iterk, best=False)

            # if subfolder is not None:
                
            #     subfolder = os.path.join(datafolder, pat, subfolder, "iter100", "k" + iterk, "")

            #     if not os.path.isdir(subfolder):
            #         subfolder = None
            #     elif not os.path.isfile(os.path.join(subfolder, "params")):
            #         subfolder = None

            # breakpoint()

            if subfolder is not None:
                try:
                    params = json.load(open(subfolder + "/params"))

                    alphas.append(params["alpha_final"])
                    rs.append(params["r_d_final"])

                    computetimes[str(params["iter_k"])].append(params["optimization_time_hours"])

                    if params["j_d_final"] < lowest:
                        lowest = params["j_d_final"]
                        lowest_counter = counter
                except json.decoder.JSONDecodeError:
                    print("vim " + subfolder.replace("Dropbox (UiO)", "Dropbox\ \(UiO\)") + "/params")
                    print("----------------------")
                    exit()

            row = []

            if counter == 0:
                row.append(r"\multirow{4}{*}{\makecell{" + pat + r"}} ")
            else:
                row.append("")
                # row += ""
            # row = pat

            row.append(" & " + dt_lookup[iterk])

            for quantity in columns:                               
                if subfolder is None:
                    row.append(" & ")
                else:

                    if quantity not in params.keys():
                        params[quantity] = np.nan

                    row.append(" & " + formatfun(params[quantity], quantity))
                
            row.append(r" \\ ")    
            patrows.append(row)
        
        mark = False
        
        assert np.nan not in alphas
        
        for a1, a2 in itertools.product(alphas[-2:], alphas[-2:]):
            if np.abs(a1-a2) / a1 > 0.1:
                mark = True
        for a1, a2 in itertools.product(rs[-2:], rs[-2:]):
            if np.abs(a1-a2) / a1 > 0.1:
                mark = True

        if mark:
            for lowest_counter in range(len(iterks)):
                bold_row = [patrows[lowest_counter][0]]

                assert patrows[lowest_counter] is not None
                for element in patrows[lowest_counter][1:-1]:
                    elements = element.split("&")

                    bold_row.append(r" & \textbf{" + elements[1] + "}")    

                bold_row.append(patrows[lowest_counter][-1])
                patrows[lowest_counter] = bold_row
    
        for patrow in patrows:
            row = ""
            for element in patrow:
                row += element
                
            all_rows.append(row)

        all_rows[-1] += r"\midrule"

    header = ""

    for h in table_header:
        header += h + " & "
    
    header = header[:-2] + r" \\ \toprule "

    latexdf = header

    # print(header)

    for row in all_rows:
        # print(row)
        latexdf += row

    if latexname is not None:
        with open(latexname, "w") as text_file:
            text_file.write(latexdf)

    



def methodtable(iterk, savepath=None, latexname=None):

    pats = groups["t1map"]

    digits = ".3f"

    def formatfun(x):
        if x < 1:
            x = 1e4 * 60 * x
            
        return format(x, digits)

    # resultfolders = ["opt_grayalpha_r_tests", "opt_r_tests_const"]
    # method_lookup = {"opt_grayalpha_r_tests": "avg D & avg T1", "opt_r_tests_const": "DTI & T1"}

    resultfolders = ["avgDTIavgT1", "DTI_filteredT1"]
    method_lookup = {"avgDTIavgT1": "avg D & avg T1", "DTI_filteredT1": "DTI & filtered T1"}

    columns = [("alpha_final", resultfolder) for resultfolder in resultfolders]
    columns += [("r_d_final", resultfolder) for resultfolder in resultfolders]

    header = [x + method_lookup[y] for x,y in columns]

    # print(header)

    result_dict = {}

    for pat in pats:

        row = []

        for quantity, resultfolder in columns:
            
            subfolder = os.path.join(datafolder, pat, "diffusion_reaction", resultfolder, "k" + iterk, "")
            # breakpoint()
            if not os.path.isdir(subfolder):
                subfolder = None
            elif not os.path.isfile(os.path.join(subfolder, "params")):
                subfolder = None
            if subfolder is None:
                row.append(np.nan)
            else:
                try:
                    params = json.load(open(subfolder + "/params"))
                except json.decoder.JSONDecodeError:
                    print("vim " + subfolder.replace("Dropbox (UiO)", "Dropbox\ \(UiO\)") + "/params")
                    print("----------------------")
                    exit()
                if resultfolder == "opt_grayalpha_r_tests":
                    # if not ("dti" in params.keys()):
                    #     # print(resultfolder, subfolder)
                    #     print(params["path_to_files"])
                    # assert params["nodti"]
                    assert params["concfolder"] == "FIGURES_CONC"
                    # assert not params["dti"]
                elif resultfolder == "opt_r_tests_const":
                    # if not ("dti" in params.keys()):
                    #     print(resultfolder, subfolder)
                    #     print(params["path_to_files"])
                    # assert params["dti"]
                    assert params["concfolder"] == "FIGURES_CONC_LUT"

                # app_val = formatfun(params[quantity])
                row.append(params[quantity])

                

        result_dict[pat] = ["" if np.isnan(x) else formatfun(x) for x in row]
        
      
    df = pandas.DataFrame.from_dict(result_dict, orient="index", columns=header)# .transpose()
    
    print("Methodtable")
    print(df)
    x = df["alpha_finalavg D & avg T1"].astype(float)
    y = df["alpha_finalDTI & filtered T1"].astype(float)
    print("alpha", np.mean(x), "+-", np.std(x))
    print("alpha", np.mean(y), "+-", np.std(y))

    # breakpoint()
    x2 = df["r_d_finalavg D & avg T1"].astype(float)
    y2 = df["r_d_finalDTI & filtered T1"].astype(float)
    print("r_d", np.mean(x2), "+-", np.std(x2))
    print("r_d", np.mean(y2), "+-", np.std(y2))

    
    result_dict["mean $\pm$ std"]= [
    format(np.mean(x), digits) + " $" + r"\pm " + "$ " + format(np.std(x), digits), 
    format(np.mean(y), digits) + " $" + r"\pm " + "$ " + format(np.std(y), digits),
    format(np.mean(x2), digits) +" $" + r"\pm " + "$ " + format(np.std(x2), digits),
    format(np.mean(y2), digits) +" $" + r"\pm " + "$ " + format(np.std(y2), digits),
    ]

    df = pandas.DataFrame.from_dict(result_dict, orient="index", columns=header)

    print("-----------------------------------------Methodtable--------------------------------------------")
    print(df)
    del x
    del y
    # exit()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        latexdf = df.to_latex(index=True, header=False, escape=False)

    # print(latexdf)

    latexdf = latexdf.split("\n")

    latexdf2 = ""

    for j in range(2, len(latexdf) - 1 - 2):
        entries = latexdf[j].split("&")
        tablerow = latexdf[j]
        entries = [x.split("\\")[0] for x in entries]
        # print(entries)
        
        # try:
        #     alpha96 = float(entries[1])
        #     alpha144 = float(entries[2])

        #     dalpha = abs(alpha96 - alpha144) / alpha144

        #     r96 = float(entries[3])
        #     r144 = float(entries[4])

        #     dr = abs(r96 - r144) / r144
        #     # print(entries[0],dr, dalpha)
        #     if dalpha > 0.1 or dr > 0.1:
        #         # print(entries[0],dr, dalpha)
        #         # tablerow = r" \textbf{ " + tablerow.replace(r"\\", "") + r"}" + r"\\"
        #         tablerow = r"  \textbf{ " + entries[0] + r"}" 
        #         for x in entries[1:]:
        #             tablerow +=  "&"
        #             tablerow += r" \textbf{ " + x + r"}"
        #         tablerow += r"\\"
            
        # except ValueError:
        #     print(j, entries[0], "pass")
        #     pass

        latexdf2 += tablerow
        # breakpoint()

    latexdf = latexdf2

    latexdf = latexdf.replace("\\begin{tabular}{lllllll}\n\\toprule\n", "")
    latexdf = latexdf.replace("\\\\\n\\bottomrule\n\\end{tabular}\n", "")
    latexdf = latexdf.replace("\\begin{tabular}{llllllllll}\n\\toprule\n", "")

    # print(latexdf)  
    assert "tabular" not in latexdf
    assert "begin" not in latexdf

    if latexname is not None:
        with open(latexname, "w") as text_file:
            text_file.write(latexdf)

    if savepath is not None:
        df.to_csv(savepath, index=True)

    return df







if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    default_input_path = "/home/basti/Dropbox (UiO)/Sleep/"
    default_output_path = "/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/"
    default_output_path += "Brain influx and clearance during sleep and sleep deprivation/figures/simulations/"

    parser.add_argument("--inputpath", type=str, default=default_input_path)
    parser.add_argument("--outputpath", type=str, default=default_output_path)

    parserargs = parser.parse_args()
    argparse_dict = vars(parserargs)


    path_to_files = argparse_dict["inputpath"]
    plotpath = argparse_dict["outputpath"]

    if path_to_files != default_input_path:
        assert plotpath != default_output_path

    latexpath = plotpath

    os.chdir(path_to_files)

    FS = 40

    dpi = 400
    figsize = (10, 10)

    matplotlib.rcParams["lines.linewidth"] = 2
    matplotlib.rcParams["axes.linewidth"] = 2
    matplotlib.rcParams["axes.labelsize"] = FS # "xx-large"
    matplotlib.rcParams["grid.linewidth"] = 1
    matplotlib.rcParams["xtick.labelsize"] = FS # "xx-large"
    matplotlib.rcParams["ytick.labelsize"] = FS # "xx-large"
    matplotlib.rcParams["legend.fontsize"] = FS # "xx-large"
    matplotlib.rcParams["font.size"] = FS

    sns.set_style(rc=dict(matplotlib.rcParams))
    sns.set_context(context=sns.plotting_context(), rc=dict(matplotlib.rcParams))

    fs = None

    dpi = 400
    figsize = (10, 10)


    assert "Sleep" in os.getcwd()

    sleepdep = groups["sleepdep"]

    pats = sorted(pats)

    def alphaformat(x):
        if np.isnan(x):
            return x
        else:
            return format(float(x), ".2f")

    def rformat(x):
        if np.isnan(x):
            return x
        else:
            return format(float(x)*1e5, ".2f")

    computetimes = {"96": [], "144": [],"288": [],"576": [],}
    dr = resolutiontable2(pats=[x for x in pats if float(x) < 200], latexname=latexpath + "alpha-r-convergence-table<200.tex")
    dr = resolutiontable2(pats=[x for x in pats if float(x) > 200], latexname=latexpath + "alpha-r-convergence-table>200.tex")
    
    print("Average compute times:")

    for iterk, timelist in computetimes.items():
        print(iterk, np.mean(timelist), "+-", np.std(timelist))

    df, dones = make_df()

    print(df)

    print("average red.", df["red. in J (%)"].mean(), "+-", np.std(df["red. in J (%)"]))

    print("average red.", df["j_plain"].mean(), "+-", np.std(df["j_plain"]))
    print("average red.", df["jfinal"].mean(), "+-", np.std(df["jfinal"]))

    # exit()

    iterk = "288"
    methodtable(iterk, savepath=plotpath + "alpha-r-method-tablek" + iterk +".csv", 
                latexname=latexpath + "alpha-r-method-tablek" + iterk + ".tex")

    sleep = [x for x in groups["sleep"] if x in dones]
    sleepdep = [x for x in groups["sleepdep"] if x in dones]

    dfs = df.loc[sleep]
    dfsdep = df.loc[sleepdep]

    for idx, qty in enumerate(["alpha", "r_d"]):

        print(qty, ",Sleep vs non sleep")
        fac = 1
        if qty not in ["alpha"]:
            fac = 1
        print(qty, "(Sleep group)            ", format(np.mean(dfs[qty] * fac), ".2f"), "+-", format(np.std(dfs[qty] * fac), ".2f"))
        if qty not in ["alpha", "rd/rn"]:
            print("half life", "(Sleep group)", format(np.median(r2halflife(dfs[qty])/3600), ".2f"), "+-", 
                format(np.std(r2halflife(dfs[qty])/3600), ".2f"), "(hours)")

        print(qty, "(Sleep deprivation group)", format(np.mean(dfsdep[qty] * fac), ".2f"), "+-", 
            format(np.std(dfsdep[qty] * fac), ".2f"))
        # r2halflife(r)
        if qty not in ["alpha", "rd/rn"]:
            print("half life", "(Sleep deprivation group)", format(np.median(r2halflife(dfsdep[qty])/3600), ".2f"), "+-", 
                format(np.std(r2halflife(dfsdep[qty])/3600), ".2f"), "(hours)")
        _, pval = scipy.stats.ttest_ind(dfs[qty], dfsdep[qty], axis=0, equal_var=True, nan_policy="omit")
        print("t-test, p =", format(pval, ".3f"))


    print("mean alpha", format(np.mean(df["alpha"]), ".1f"), "pm", format(np.std(df["alpha"]), ".1f"))
    print("mean r", format(np.mean(df["r_d"]), ".0f"), "pm", format(np.std(df["r_d"]), ".0f"))

    df.to_csv(plotpath + "dframe.csv", index=True)

    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("Values excluding pat 091")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    
    # breakpoint()
    df2 = df.drop('091')
    
    del df

    dfs = df2.loc[[x for x in sleep if x != "091"]]
    dfsdep = df2.loc[sleepdep]

    for idx, qty in enumerate(["alpha", "r_d"]):

        print(qty, ",Sleep vs non sleep")
        fac = 1
        if qty not in ["alpha"]:
            fac = 1
        print(qty, "(Sleep group)            ", format(np.mean(dfs[qty] * fac), ".2f"), "+-", format(np.std(dfs[qty] * fac), ".2f"))
        if qty not in ["alpha", "rd/rn"]:
            print("half life", "(Sleep group)", format(np.median(r2halflife(dfs[qty])/3600), ".2f"), "+-", 
                format(np.std(r2halflife(dfs[qty])/3600), ".2f"), "(hours)")

        print(qty, "(Sleep deprivation group)", format(np.mean(dfsdep[qty] * fac), ".2f"), "+-", 
            format(np.std(dfsdep[qty] * fac), ".2f"))
        # r2halflife(r)
        if qty not in ["alpha", "rd/rn"]:
            print("half life", "(Sleep deprivation group)", format(np.median(r2halflife(dfsdep[qty])/3600), ".2f"), "+-", 
                format(np.std(r2halflife(dfsdep[qty])/3600), ".2f"), "(hours)")
        _, pval = scipy.stats.ttest_ind(dfs[qty], dfsdep[qty], axis=0, equal_var=True, nan_policy="omit")
        print("t-test, p =", format(pval, ".3f"))


    print("mean alpha", format(np.mean(df2["alpha"]), ".1f"), "pm", format(np.std(df2["alpha"]), ".1f"))
    print("mean r", format(np.mean(df2["r_d"]), ".0f"), "pm", format(np.std(df2["r_d"]), ".0f"))

    print("median alpha", format(np.median(df2["alpha"]), ".1f"), "pm", format(np.std(df2["alpha"]), ".1f"))
    print("median r", format(np.median(df2["r_d"]), ".0f"), "pm", format(np.std(df2["r_d"]), ".0f"))
    print("median half life", format(np.median(df2["r_d"]), ".0f"), "pm", format(np.std(df2["r_d"]), ".0f"))


    print("min alpha", format(np.min(df2["alpha"]), ".1f"), " max ", format(np.max(df2["alpha"]), ".1f"))
    print("min r", format(np.min(df2["r_d"]), ".0f"), " max ", format(np.max(df2["r_d"]), ".0f"))

    del dfs   

    qtyname = "r"

    
    # barplot(qty="r_d", qtyname="r", savepath=figurepath +"r.png", twinax=False, dpi=dpi, figsize=figsize)
    # barplot(qty="alpha", qtyname=r"\alpha", savepath=figurepath + "alpha.png", twinax=False, dpi=dpi, figsize=figsize)


    paperformat = True

    pats.remove("091")

    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("Creating plots excluding pat 091")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------")
    
    # make_figs(region="white", pats=pats, alphas=alphas, data_folder=datafolder, average_tracer=False)

    width = 0.8 + 0.1


    def resultfoldername(pat, alpha):
        return reaction_resultfolder(pat, best=True, k=None, subfoldername="avgDTIavgT1")

    pers2permin = 60 * 1e4

    fig2, ax2 = plt.subplots(figsize=figsize, dpi=dpi)

    # fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    bestalpha, bestr = [], []

    for pat in pats:
        foldername = resultfoldername(pat, alpha=None)

        data = json.load(open(foldername + "params"))

        bestalpha.append(data["alpha_final"])
        bestr.append(data["r_d_final"])

        # ax.plot([0,1], [data["alpha_final"], 1e5 * data["r_d_final"]], marker="o")


        ax2.plot(data["alpha_final"], pers2permin * data["r_d_final"], marker="o", markersize=20, color="darkred")

    # ax.set_xlim(-0.1, 1.1)
    # ax.tick_params(axis='x', width=0)
    # ax.tick_params(axis='y', labelsize=fs)
    # plt.xticks([0, 1], [r"$\alpha$", "$r \, (10^{-5}$s$^{-1}$)" ], fontsize=FS)
    # plt.yticks([1, 3, 5, 7, 9])
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')


    
    # plt.tight_layout()
    plt.savefig(plotpath + "bests.png", dpi=400)
    

    test = scipy.stats.pearsonr(bestalpha, bestr)
    
    print("tid best alpha, r pearson r:", format(test[0], ".2f"))

    
    plt.sca(ax2)
    #ax2.tick_params(axis='x', width=0)
    # ax2.tick_params(axis='y', labelsize=fs)

    x = np.linspace(min(bestalpha), max(bestalpha), 100)
    y = np.mean(bestr) + test[0] * np.std(bestalpha) * np.std(bestr) * (x - np.mean(bestalpha))
    plt.plot(x, pers2permin * y, color="k", label='$r=$' + format(test[0], ".2f") + ", $p=$" + format(test[1], ".2f"))
    
    plt.xticks(range(1,8))
    
    plt.yticks(np.array(range(1, 8)) * 10)
    
    plt.xlabel(r"best $\alpha$", fontsize=FS)
    plt.ylabel("best $r \, (10^{-4}$min$^{-1}$)", fontsize=FS)

    # plt.title('Correlation r=' + format(test[0], ".2f"), fontsize=FS)
    
    plt.legend(loc="upper left", 
            bbox_to_anchor=(0., 1.0, 1., .07),
            fontsize=matplotlib.rcParams["legend.fontsize"]-2)
    dw = 0.04

    plt.tight_layout(rect=[-dw, -dw, 1 + dw, 1 + dw])

    plt.savefig(plotpath + "bestscorrelation.png", dpi=400)
    
    # plt.show()

    # exit()

    alphas = [1, np.nan]

    if np.nan in alphas:
        suffix = "best"


    GREY_WHITE = False

    if not GREY_WHITE:
        plotname = "colors"

    for region in ["white", "gray"]:

        ylabel = ylabels[region]

        # # display best vs. plain for every patient in a single plot
        # make_figs(region, pats, alphas, paperformat, 
        #             resultfoldername=resultfoldername, data_folder=datafolder, fs=fs,
        #             savepath=plotpath + "barplot" + region + ".png", 
        #             figsize=figsize,
        #             dpi=dpi)


        make_barplot(region, pats, alphas, paperformat, width=width, ylabel=ylabel,
                    resultfoldername=resultfoldername, data_folder=datafolder, FS=FS,
                    savepath=plotpath + plotname + "barplot" + region + ".png", GREY_WHITE=GREY_WHITE,
                    figsize=figsize,
                    dpi=dpi)