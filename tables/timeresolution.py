
# salloc --ntasks=1 --mem-per-cpu=1G --time=00:59:00 --account=NN9279K --qos=devel

import os
import json
import pandas
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--quantity", choices=["relL2err", "c48"])

parserargs = vars(parser.parse_args())

quantity = parserargs["quantity"]

intervals = [
    # (0, 0.001),
    (1.2, 2.6), 
    (4.5, 7.4), (20, 30), (40, 56)
]


def get_data_in_intervals(pat, stored_times, stored_data):

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
            print("No data available in interval", intervals[idy], "for", pat, ", appending nan")

    return ts, data


datapath = "/cluster/projects/nn9279k/bastian/SleepData/"

iterks = ["96", "144", # "192", 
            "288", "576"]

dts = [format(60 * 48 / (float(x)), ".0f") for x in iterks]

dt_lookup = {}
for i, dt in zip(iterks, dts):
    dt_lookup[i] = dt

def make_table(columns, resolutions):
    
    nc = ""
    for _ in resolutions:
        nc += "c"
    print(r"\begin{tabular}[t]{c|c|" + nc + r"}")

    header = r"$(\alpha,r)$ & $dt$ "
    for res in resolutions:
        header += " & $n=" + str(res) + "$"
    header += r"\\"
    print(header)
    print(r"\toprule")

    printed = False

    for alpha, r in columns:

        for counter, iterk in enumerate(iterks):
            
            exp_conc = r"\multicolumn{2}{c}{data}"


            if counter == 0:
                row = r"\multirow{" + str(len(iterks)) + "}{*}{\makecell{(" + str(alpha) + "," +str(r) + r")}} & "
            else:
                row = " & "

            row += dt_lookup[iterk] + " & "

            for resolution in resolutions:

                path = datapath + pat + "/alpha_r_tests_" + resolution + "/"

                assert os.path.isdir(path)

                folders = os.listdir(path)

                f = "alpha" + str(alpha)

                f = path + f + "/r" + str(r) + "/"
                
                if not os.path.isdir(f):
                    row += " x & "
                    continue

                for iterpath in os.listdir(f):
                    
                    if not str(iterk) in iterpath:
                        continue

                    iterpath = f + iterpath

                    try:
                        if not os.path.isfile(iterpath + "/params"):
                            # print(iterpath + "/params", "not found")
                            row += "XXX & "
                        
                        else:

                            if quantity == "relL2err":
                                d1 = json.load(open(iterpath + "/params"))
                                rowval = format(d1["j_d_final"] * 100, ".1f")
                            elif quantity == "c48":


                            
                                sim = pandas.read_csv(iterpath + "/concs.csv")
                                exp = pandas.read_csv(iterpath + "/experimental_data.csv")
                                ts, simd = get_data_in_intervals(pat, sim["t"], sim["avg"])
                                ts, expd = get_data_in_intervals(pat, exp["t"], exp["avg"])
                                rowval = format(simd[-1], ".3f")
                                exp_conc += " & " + format(expd[-1], ".3f")

                            row += rowval + " & "

                            

                    except json.decoder.JSONDecodeError:
                        print("vim " + iterpath + "/params")

            if not printed and quantity == "c48":
                print(exp_conc + r"\\ \midrule")
                printed = True
            row = row[:-2] + r"\\"
            print(row)


        print(r"\midrule")

resolutions = ["16", "32", "64"]

iterks = ["96", "144", "288"]

r = 0
rr = "5e-5"
columns = [(1, r), (3, r), (5, r),
            (1, rr), (3, rr), (5, rr)
            ]
del r



print(r"\begin{table}")

if quantity == "c48":
    print(r"\caption{Concentration (mmol / L) at 48 hours, data and simulation}")

elif quantity == "relL2err":
    print(r"\caption{Rel. L2 error in $\%$ between data and simulations for different parameters $(\alpha, r)$ and different time/mesh resolutions.}")

for pat in ["105", "205"]:

    print(r"\begin{subtable}[t]{0.5\textwidth}")
    print(r"\centering")

    make_table(columns=columns, resolutions=resolutions)


    print(r"\end{tabular}")
    print(r"\caption{Pat ID " + pat + "}")
    print(r"\end{subtable}")


print(r"\end{table}")
# print(r"")
# print(r"")
