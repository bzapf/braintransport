import sys
import os
import glob
import datetime
import csv
import itertools
from dolfin import *

def extract_relative_times(patient_dir):
    image_dir = os.path.join(patient_dir, "FIGURES_CONC_LUT")
    if not os.path.isdir(image_dir):
        print("Missing image directory %s" % image_dir)
        return []

    images = glob.glob(os.path.join(image_dir, "*.mgz"))

    times = []
    for image in images:
        date, time, _ = os.path.basename(image).split("_")
        d = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
        t = datetime.time(int(time[0:2]), int(time[2:4]), int(time[4:6]))
        times += [datetime.datetime.combine(d, t)]
    times.sort()
    t0 = times[0]
    relatives = [(t - t0) for t in times]
    
    return relatives

def list_times_for_all_patients():
    print("Extracting relative examination times for all patients.")
    workdir = "/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer"
    patient_dirs = glob.glob(os.path.join(workdir, "[0-9][0-9][0-9]"))
    patient_dirs.sort()
    print(patient_dirs)
    sec2h = 60*60

    df = []
    for d in patient_dirs:
        patient_num = os.path.basename(d)
        dts = extract_relative_times(d)

        df += [[int(patient_num)] + [dt.total_seconds()/sec2h for dt in dts]]

    with open("results/times_per_patient.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(df)     

def extract_mesh_info(patient_dir, n):

    mesh = Mesh()
    meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d_with_DTI.h5" % n)
    meshfile2 = os.path.join(patient_dir, "mesh", "parenchyma%d.h5" % n)
    if os.path.isfile(meshfile):
        print("Reading mesh from %s" % meshfile)
    elif os.path.isfile(meshfile2):
        print("Reading mesh from %s" % meshfile2)
        meshfile = meshfile2
    else:
        return (None, 0, 0, 0, 0, 0)
    
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)
    hdf.close()
    
    return (mesh, mesh.num_vertices(), mesh.num_edges(), mesh.num_cells(),
            mesh.hmin(), mesh.hmax())

def list_dims_for_all_patients():

    print("Extracting mesh dimensions for all patients.")
    workdir = "/cluster/projects/nn9279k/Vegard/SleepDeprivation/Freesurfer"
    patient_dirs = glob.glob(os.path.join(workdir, "[0-9][0-9][0-9]"))
    patient_dirs.sort()
    print(patient_dirs)
    
    dims = []
    for patient_dir in patient_dirs:
        patient_num = os.path.basename(patient_dir)
        
        dim = [patient_num]
        for n in [16, 32, 64]:
            meshinfo = extract_mesh_info(patient_dir, n)
            if not meshinfo:
                print("Unable to extract info for %s with n %d" % (patient_num, n))
                dim += [0,]*4
            else:
                (mesh, v, e, c, hmin, hmax) = meshinfo
                print("... with %d vertices, and %d cells" % (v, c))
                dim += [v, c, "%.2f" % hmin, "%.2f" % hmax]
                
        dims += [dim]
        for d in dims:
            print(d)

    with open("results/dims_per_patient.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(dims)     

def analyze_optimization(csvname):
    """
    csvname: File with optimization values output
    """
    # Read Js and phi_max'es from the output file
    js = []
    phi_maxs = []
    with open(csvname) as csvfile:
        reader = csv.reader(csvfile)
        for row in list(reader)[1::2]:
            (j, phi_max) = row
            js += [float(j)]
            phi_maxs += [float(phi_max)]

    return js, phi_maxs

def list_avg_phi_for_all_patients(timekey, outfilename, resultfoldername):
    #print("\% Extracting \bar{\phi} for all patients and simulations for from %s." % xh)
    workdir = "/cluster/projects/nn9279k/meg/modelling-sleep-deprivation-2021/" + resultfoldername
    patient_dirs = glob.glob(os.path.join(workdir, "[0-9][0-9][0-9]"))
    patient_dirs.sort()
    
    ns = [16, 32]
    betas = ["1.0e-03", "1.0e-04", "1.0e-05"]

    # These patients are excluded due to missing MR scans - or
    # non-converged optimization ("172", "191", "205")
    excluded = {"0h": None, "6h": None, "24h": None}
    excluded["0h"] = set(("091", "241"))
    excluded["6h"] = set(("091", "241"))
    excluded["24h"] = set(("002", "078", "127", "172", "191", "205"))

    # Create custom table with relevant information
    table = {}
    for patient_dir in patient_dirs:

        # Skip excluded patients
        ID = os.path.basename(patient_dir)
        if ID in excluded[timekey]:
            continue

        table[ID] = {}
        for n in [16, 32]:
            table[ID][n] = []
            for beta in betas:
                key = "results_red_ocd_pat%s_n%d_beta%s_%s" % (ID, n, beta, timekey)
                results_dir = os.path.join(workdir, ID, "simulations", key)
                if os.path.isdir(results_dir):
                    values_csv = os.path.join(results_dir, "values.csv")
                    with open(values_csv) as csvfile:
                        reader = csv.reader(csvfile)
                        rows = [row for row in reader]
                        avg_phi = float(rows[1][1]) # Average phi (mm/h)
                        c_L2 = float(rows[1][3])    # L2-norm of computed c
                        c1_L2 = float(rows[1][4])   # L2-norm of target c1
                        delta_c = abs(c_L2 - c1_L2)/c1_L2*100 # Relative diff in percent
                        table[ID][n] += [avg_phi, c_L2, c1_L2, delta_c]
                else:
                    table[ID][n] += [0.]*4

    outfile = open(outfilename, "w")

    # Custom print the table
    IDs = sorted(table.keys())
    for ID in IDs:
        row = ["%s" % ID, "16"]
        row += ["%.2f" % i for i in table[ID][16]]
        outfile.write(" & ".join(row) + " \\\\")
        outfile.write("\n")

        row = ["", "32"]
        row += ["%.2f" % i for i in table[ID][32]]
        outfile.write(" & ".join(row) + " \\\\")
        outfile.write("\n")

        #print("\\midrule")
        
    outfile.close()


def list_Js_for_all_patients(timekey, outfilename, resultfoldername):

    #print("Extracting Js for all patients and simulations.")
    workdir = "/cluster/projects/nn9279k/meg/modelling-sleep-deprivation-2021/" + resultfoldername
    patient_dirs = glob.glob(os.path.join(workdir, "[0-9][0-9][0-9]"))
    patient_dirs.sort()
    
    outfile = open(outfilename, "w")

    dims = []
    times = [timekey,]#, "0h", "6h", "24h"]
    ns = [16, 32]
    betas = ["1.0e-03", "1.0e-04", "1.0e-05"]

    outfile.write("ID, t0 & $n$ & $J_0$ & $J_N$ & $\max \phi_{N}$ ($\Delta \phi$\%) & $J_0$ & $J_N$ & $\max \phi_{N}$ ($\Delta \phi$ \%) & $J_0$ & $J_N$ & $\max \phi_{N}$ ($\Delta \phi$ \%) \\\\")
    outfile.write("\n")

    max_delta_phi = 0
    min_delta_j = 100000
    max_delta_j = 0
    
    for t0 in times:
        for patient_dir in patient_dirs:
            outfile.write("\midrule")
            outfile.write("\n")
            patient_num = os.path.basename(patient_dir)
            sim_dir = os.path.join(workdir, patient_num, "simulations")
            corner = patient_num + ", " + t0
            for n in [16, 32]:
                row = [corner + " & " + str(n)]
                for beta in betas:
                    key = "results_red_ocd_pat%s_n%d_beta%s_%s" % \
                        (patient_num, n, beta, t0)
                    res_dir = os.path.join(sim_dir, key)
                    if os.path.isdir(res_dir):
                        opt_csv = os.path.join(res_dir, "opts",
                                               "optimization_values.csv")
                        
                        js, max_phis = analyze_optimization(opt_csv)
                        j_0 = js[0]
                        j_N = js[-1]
                        delta_j = j_0/j_N
                        phi_N_1 = max_phis[-2]
                        phi_N = max_phis[-1]
                        delta_phi = (abs(phi_N - phi_N_1))/phi_N*100
                        row += ["%.1f" % j_0, "%.2f" % j_N,
                                "%.2f (%.2f)" % (phi_N, delta_phi)]
                        if delta_phi > max_delta_phi:
                            max_delta_phi = delta_phi
                        if delta_j > max_delta_j:
                            max_delta_j = delta_j
                        if delta_j < min_delta_j and delta_j > 6:
                            min_delta_j = delta_j
                    else:
                        row += ["-", "-", "- (-)"]

                row_str = " & ".join(row) + " \\" + "\\"
                outfile.write(row_str)
                outfile.write("\n")
                corner = ""

        if True:
            outfile.write("%% max_delta_phi = " + str(max_delta_phi))
            outfile.write("\n")
            outfile.write("%% max_delta_j = " + str(max_delta_j))
            outfile.write("\n")
            outfile.write("%% min_delta_j > 6 = " + str(min_delta_j))
            outfile.write("\n")
                
if __name__ == "__main__":

    ## Run these commands first at Saga in the current directory
    # salloc --ntasks=1 --mem-per-cpu=8G --time=00:59:00 --account=NN9279K --qos=devel
    # source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

    ## Then run the script at will
    # python3 misc.py

    os.makedirs("results/tables", exist_ok=True)

    # # To generate table of times for all patients (mainly used pre-analysis).
    list_times_for_all_patients()

    # # To generate table of mesh dimensions and sizes:
    list_dims_for_all_patients()


    for timekey in ["0h", "6h", "24h"]:

        list_Js_for_all_patients(timekey, outfilename="results/tables/table_S4_" + timekey + ".tex", 
            resultfoldername="results-red-ocd-avg-dti")
        list_avg_phi_for_all_patients(timekey, outfilename="results/tables/table_S5_" + timekey + ".tex", 
            resultfoldername="results-red-ocd-avg-dti")


    # To generate table of J, max \phi etc for all patients:
    # run python3 misc.py 0h > results/tables/table_S4_0h.tex

    if False:
        xh = str(sys.argv[1])
        list_Js_for_all_patients(xh)

    # To generate table of average velocity for all patients etc.
    # run python3 misc.py 0h > results/tables/table_S5_0h.tex

    
    if False:
        xh = str(sys.argv[1])
        list_avg_phi_for_all_patients(xh)

