import sys
import os
import csv
import datetime as dt
from collections import OrderedDict

def read_all_div_phi_values(pat, pat_dir, n, t, beta, debug=False):
    # Read and collect average div(phi) brain-wide, gray, matter, white matter and stem
    name = "results_red_ocd_pat%s_n%s_beta%s_%s" % (pat, n, beta, t)
    filename = os.path.join(pat_dir, name, "div_phi_avgs.csv")
    try:
        with open(filename) as csvfile:
            reader = csv.reader(csvfile)
            rows = [row for row in reader]
            avg_phis = [float(r) for r in rows[1]]
    except:
        if debug:
            print("Failed to read %s, continuing" % filename)
        avg_phis = None
            
    return avg_phis

def read_all_avg_phi_values(pat, pat_dir, n, t, beta, debug=False):
    # Read and collect average phi values in total, gray, matter, white matter and stem
    name = "results_red_ocd_pat%s_n%s_beta%s_%s" % (pat, n, beta, t)
    filename = os.path.join(pat_dir, name, "phi_mag_avgs.csv")
    try:
        with open(filename) as csvfile:
            reader = csv.reader(csvfile)
            rows = [row for row in reader]
            avg_phis = [float(r) for r in rows[1]]
    except:
        if debug:
            print("Failed to read %s, continuing" % filename)
        avg_phis = None
            
    return avg_phis

def read_avg_phi_values(pat, pat_dir, n, t, betas, debug=False):
    # Read and collect average phi values printed during simulation
    # runs for given patient, given n, given time interval t, but all
    # betas
    avg_phis = []
    for beta in betas:
        name = "results_red_ocd_pat%s_n%s_beta%s_%s" % (pat, n, beta, t)
        filename = os.path.join(pat_dir, name, "values.csv")
        try:
            with open(filename) as csvfile:
                reader = csv.reader(csvfile)
                rows = [row for row in reader]
                avg_phis += [float(rows[1][1])]
                
        except:
            if debug:
                print("Failed to read %s, continuing" % filename)
                avg_phis += [0.0]
            
    return avg_phis

# def extract_avg_phis(xh, datadir, patients):

#     # xh is "0h" or "6h" or "24h"
#     n = "32"  # "16" or "32", latter is probably more relevant
#     betas = ["1.0e-03", "1.0e-04", "1.0e-05"]
    
#     print("# Extracting from %d patients: " % len(patients), patients)
#     print("# Reading from %s" % datadir)
#     now = dt.datetime.now(dt.timezone.utc)
    
#     results = {}
#     for pat in patients:
#         pat_dir = os.path.join(datadir, pat, "simulations")
#         avg_phis = read_avg_phi_values(pat, pat_dir, n, xh, betas)
#         if avg_phis:
#             results[pat] = avg_phis
#     print("# n = %s, t (xh) = %s" % (n, xh))
#     print("# Successfully collected values from %s subjects" % len(results.keys()))
#     print("# Extracted at %s (UTC)" % now)
#     print("results_%s = " % n, results)

def extract_regional_avg_phis(xh, datadir, patients):

    # xh is "0h" or "6h" or "24h"
    n = "32"  # "16" or "32", latter is probably more relevant
    beta = "1.0e-04"
    
    print("# Extracting from %d patients: " % len(patients), patients)
    print("# Reading from %s" % datadir)
    now = dt.datetime.now(dt.timezone.utc)
    
    results = OrderedDict({})
    for pat in patients:
        pat_dir = os.path.join(datadir, pat, "simulations")
        avg_phis = read_all_avg_phi_values(pat, pat_dir, n, xh, beta)
        if avg_phis:
            results[pat] = avg_phis
    print("# n = %s, t (xh) = %s" % (n, xh))
    print("# Successfully collected values from %s subjects" % len(results.keys()))
    print("# Extracted at %s (UTC)" % now)
    print("# avg |v|, avg |v|_gray, avg |v|_white, avg |v|_stem (mm/h)")
    print("from collections import OrderedDict)
    print("phi_averages_%s = " % xh, results)

def extract_div_phis(xh, datadir, patients):

    n = "32"  # "16" or "32", latter is probably more relevant
    beta = "1.0e-04"
    
    print("# Extracting from %d patients: " % len(patients), patients)
    print("# Reading from %s" % datadir)
    now = dt.datetime.now(dt.timezone.utc)
    
    results = OrderedDict({})
    for pat in patients:
        pat_dir = os.path.join(datadir, pat, "simulations")
        div_phis = read_all_div_phi_values(pat, pat_dir, n, xh, beta)
        if div_phis:
            results[pat] = div_phis
    print("# n = %s, t (xh) = %s" % (n, xh))
    print("# Successfully collected values from %s subjects" % len(results.keys()))
    print("# Extracted at %s (UTC)" % now)
    print("# avg div(v), avg div(v)_gray, avg div(v)_white, avg div(v)_stem (1/h)")
    print("from collections import OrderedDict)
    print("div_phi_averages_%s = " % xh, results)

if __name__ == "__main__":

    ## Run these commands first at Saga in the current directory
    # salloc --ntasks=1 --mem-per-cpu=8G --time=00:59:00 --account=NN9279K --qos=devel
    # source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

    # Where do we want to read patient data from
    datadir = "/cluster/projects/nn9279k/meg/modelling-sleep-deprivation-2021/"
    datadir = os.path.join(datadir, "results-red-ocd-avg-dti")

    # Which patient directories do we want
    patients = ["002", "078", "091", "105", "127", "172", "175", "176",
                "178", "183", "190", "191", "199", "205", "215", "218",
                "227", "228", "230", "235", "236", "240", "241", "249"]
    patients.sort()

    # Set False -> True below, and run with:
    # python3 read_ocd_values.py 6h > results/ocd_averages_6h.py
    # python3 read_ocd_values.py 24h > results/ocd_averages_24h.py
    xh = sys.argv[1]
    if False:
        extract_regional_avg_phis(xh, datadir, patients)

    # OR set False -> True below, and run with:
    # python3 read_ocd_values.py 6h > results/div_phi_averages_6h.py
    # python3 read_ocd_values.py 24h > results/div_phi_averages_24h.py
    if True:
        extract_div_phis(sys.argv[1], datadir, patients)
