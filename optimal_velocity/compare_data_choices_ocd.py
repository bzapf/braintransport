import os
import pandas
import numpy as np
import pandas
import pickle
import copy 

# conda activate /cluster/projects/nn9279k/bastian/mri_inv_env

print("Starting to run script")
try:
    import parse
except ModuleNotFoundError:
    pass

def get_params_from_foldername(resfolder):
    parse.parse("results_red_ocd_pat{}_n{}_beta{}_{}", resfolder)
    pat = parse.parse("results_red_ocd_pat{}_n{}_beta{}_{}", resfolder)[0]
    n = parse.parse("results_red_ocd_pat{}_n{}_beta{}_{}", resfolder)[1]
    beta = parse.parse("results_red_ocd_pat{}_n{}_beta{}_{}", resfolder)[2]
    timekey = parse.parse("results_red_ocd_pat{}_n{}_beta{}_{}", resfolder)[3]

    return pat, n, beta, timekey


path = "/cluster/projects/nn9279k/meg/modelling-sleep-deprivation-2021/"
targetpath = "/cluster/home/bazapf/nobackup/braintransport/optimal_velocity/results/"

fa1 = "results-red-ocd-dti"
fb1 = "results-red-ocd-avg-dti"

labels = ["T1&DTI", "avgT1&avgD"]

qtys = [# "phi_avg (mm/h)", 
    "phi_mag_avgs",
    "div_phi_avgs"
    # "|div phi|_L2"
    ]

rois = {}

for q in qtys:
    if q == "phi_mag_avgs":
        rois[q] = ["avg (mm/h)", "avg_gray" ,"avg_white", "avg_stem"]
    elif q == "div_phi_avgs":
        rois[q] = ["avg (1/h)", "avg_gray" ,"avg_white", "avg_stem"]

roinames = ["avg", "avg_gray" ,"avg_white", "avg_stem"]



header = []

for q in qtys:
    for roi in roinames:
        for m in labels:
            header.append(q + " (" + m + ")" + roi)

replaced_header = False

t1pats= ["105", "175", "178", "183", "190", "191", "205", "215", "228", "240", "199", "227", "235", "236", "241"]

compare_n = "32"
compare_beta = "1.0e-04"


vols = {}
dfs = {}

mean_df = {}

timekeys = ["6h", "24h"]

for compare_t in timekeys:

    result_dict = {}


    for pat in t1pats:

        if compare_t == "24h":
            if pat == "191" or pat == "205":
                print("Exclude", pat)
                continue

        volpath = os.path.join("/cluster/projects/nn9279k/bastian/SleepData/", pat, "region_volumes32")
        
        vols[pat] = pickle.load(open(volpath, "rb"))["avg"]

        result_dict[pat] = []

        fa = os.path.join(path, fa1, pat, "simulations", "")
        fb = os.path.join(path, fb1, pat, "simulations", "")

        for key in qtys:

            for readkey in rois[key]:
            
                for resfolder in [x for x in os.listdir(fa) if "log" not in x]:

                    pat, n, beta, timekey = get_params_from_foldername(resfolder)
                    
                    if not (n == compare_n and beta == compare_beta and timekey == compare_t):
                        continue

                    if resfolder in os.listdir(fb):
                                
                        try:

                            if "div_phi_avgs" in key or "phi_mag_avgs" in key:
                                resa = os.path.join(fa, resfolder, key + ".csv")
                                resb = os.path.join(fb, resfolder, key + ".csv")
                            
                                valsa = pandas.read_csv(resa)
                                valsb = pandas.read_csv(resb)

                                valsa = valsa[readkey][0]
                                valsb = valsb[readkey][0]

                                # assert valsa is not None
                                # print(valsa, valsb)
                                # exit()

                                if "phi_mag_avgs" in key:
                                    valsa = 1e3 * valsa / (60)
                                    valsb = 1e3 * valsb / (60)
                                elif "div_phi_avgs" in key:
                                    valsa = 1e4 * valsa / 60
                                    valsb = 1e4 * valsb / 60


                                    # if not replaced_header:
                                    #     header = [x.replace("div_phi_avgs", "influx (1e-4 / min)") for x in header]
                                    #     replaced_header = True

                            else:
                                raise ValueError
                                resa = os.path.join(fa, resfolder, "values.csv")
                                resb = os.path.join(fb, resfolder, "values.csv")
                            
                                valsa = pandas.read_csv(resa)
                                valsb = pandas.read_csv(resb)

                                valsa = valsa[key][0]
                                valsb = valsb[key][0]

                                assert "div_phi_avgs" not in key

                                if "phi_avg" in key:

                                    valsa = 1e3 * valsa / 60
                                    valsb = 1e3 * valsb / 60

                                    if not replaced_header:
                                        header = [x.replace("mm/h", "mum/min") for x in header]
                                        replaced_header = True


                            
                            # result_dict[pat].append(format(valsa, ".0f"))
                            # result_dict[pat].append(format(valsb, ".0f"))
                            result_dict[pat].append(valsa)
                            result_dict[pat].append(valsb)

                        except FileNotFoundError:
                            # print("vim " + resa)
                            # exit()
                            if n == compare_n and beta == compare_beta and timekey == compare_t:
                                result_dict[pat].append(np.nan)
                                result_dict[pat].append(np.nan)
                            continue

                    else:
                        # pass
                        print("--", pat, resfolder, "not in results-b")
                        print(pat, n, beta, timekey)


    # print(result_dict)
    df = pandas.DataFrame.from_dict(copy.deepcopy(result_dict), orient="index", columns=header)
    dfs[compare_t] = df
    print(df)

# exit()

header = []
methods = [" (T1&DTI)", " (avgT1&avgD)"]
qtys = ["phi_mag_avgs", "div_phi_avgs"]
subheader = []

mean_df = {}

for timekey in timekeys:
    # for method in methods:
    
    header.append(timekey)
    header.append(timekey)

    subheader.append("(T1&DTI)")
    subheader.append("(avgT1&avgD)")


def myformat(x):
    return format(x, ".1f")


for q in qtys:
    for roi in roinames:
        row = []
        for timekey in timekeys:
            for method in methods:
                
                key = q + method + roi

                # print(dfs[timekey])

                vals = dfs[timekey][key]
                # print(vals)
                # exit()

                cellval = format(np.nanmean(vals), ".2f")

                row.append(cellval)

        # print(row, header)
        mean_df[q + roi] = row
# 
# print(result_dict["105"])
mean_df = pandas.DataFrame.from_dict(mean_df, orient="index", columns=header)
print(mean_df)

mean_df.to_csv(targetpath + "ocd_dti_vs_nodti.csv", index=True)