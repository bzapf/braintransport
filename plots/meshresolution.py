from fenics import *
import os
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd

from definitions import datafolder, intervals

from helpers import extract_data, get_data_in_intervals, sim_at_mri_times

pat = "205"
pat = "105"

for pat in ["105", "205"]:

    c16, c32, c64 = None, None, None

    alpha=1
    r="0"

    resfolders = {}

    resultfolder = datafolder + str(pat) + "/convergence_tests/"


    for meshn in ["16", "32", "64"]:
        for timesteps in ["96", "144", "288"]:
            resfolders[(meshn, timesteps)] = resultfolder + "alpha_r_tests_" + meshn + "/alpha" + str(alpha) + "/r" + r + "/k" + timesteps + "/"

    dtlookup = {"96": "30","144": "20", "288": "10"}

    linestyles = {"96": "-","144": "-.", "288": ":"}

    volumes_avg = {}
    if pat == "105":

        vertices_cells_hmax = {16: (76920, 392357, 10.991521987123088), 32: (229201, 1085103, 5.511652460497721), 64: (927444, 423353,2.817558781723474)}
    elif pat == "205":
        volumes_avg = {"16": 1214932, "32": 1207988, "64": 1209862.}
        volumes_white = {"16": 565397.7643429264, "32": 563620.135340711, "64": 562008.1748970981}
        volumes_gray = {"16": 629790.6226778422, "32": 624661.7064872711, "64": 628180.7425056563}
        vertices_cells_hmax = {16: (82439, 421944, 10.92), 32: (247904, 1179462, 5.47), 64: (950203, 4453739, 2.75)}


    # volumes = {"avg": volumes_avg, "white": volumes_white, "gray": volumes_gray}

    assert len(resfolders) > 0

    # for meshfile in ["parenchyma16.h5", "parenchyma32.h5", "parenchyma64.h5"]:
    #     meshpath = "/home/basti/Dropbox (UiO)/Sleep/" + str(pat) + "/mesh/" + meshfile

    #     if not os.path.isfile(meshpath):
    #         continue
    #     else:
    #         print(meshpath)
    #     mesh = Mesh()
    #     hdf = HDF5File(mesh.mpi_comm(), meshpath, "r")
    #     hdf.read(mesh, "/mesh", False)
    #     SD = MeshFunction("size_t", mesh, mesh.topology().dim())
    #     hdf.read(SD, "/subdomains")

    #     print(mesh.num_vertices(), mesh.num_cells(), mesh.hmax())

    #     dx_SD = Measure('dx')(domain=mesh, subdomain_data=SD)
    #     vol = assemble(Constant(1)*dx_SD(domain=mesh))
    #     gray_volume = assemble(Constant(1)*dx_SD(1))
    #     white_volume = assemble(Constant(1)*dx_SD(2))

    #     print(meshfile, vol)
    #     print("white", white_volume)
    #     print("gray", gray_volume)

    # colors_ = list(cm.jet(np.linspace(0, 1, len(dts))))
    # colors = {}
    # for idx, dt in enumerate(dts):
    #     colors[dt] = colors_[idx]

    # colors_ = list(cm.brg(np.linspace(0, 1, 3)))
    # # colors_ = seaborn.color_palette("crest", 3)

    # colors = {}
    # for idx, key in enumerate(["16", "32", "64"]):
    #     colors[int(key)] = colors_[idx]

    colors = {16: "orangered", 32: "forestgreen", 64: "mediumblue"}

    # colors = {200: "red", 600: "green", 1800: "blue", 3200: "yellow", 5400: "k"}

    plt.figure(dpi=222)
    fs = 15
    lw = 1.5

    c48_list32 = []
    c48_list64 = []

    # for idr, roi in enumerate(["avg", "white", "gray"]):
    for idr, roi in enumerate(["avg"]):
        for idx, key in enumerate(sorted(resfolders.keys(), key=lambda x: x[1])):

            resfolder = resfolders[key]

            exceltable = pd.read_csv(os.path.join(resfolder, "concs.csv"))
            experimentaltable = pd.read_csv(os.path.join(resfolder, "experimental_data.csv"))
            with open(os.path.join(resfolder, "params")) as f:
                data = json.load(f)

            c = "k"
            v = 1

            if int(key[0]) == 16 and int(key[1]) == 144:
                assert c16 is None
                c16 = exceltable["avg"]


                simulation_times, simulated_tracer_at_times = sim_at_mri_times(pat, mritimes=experimentaltable["t"] / 3600, simdata=exceltable["avg"], simtimes=exceltable["t"] / 3600)    

                
                _, c16intervals = get_data_in_intervals(pat="205", stored_times=simulation_times, stored_data=simulated_tracer_at_times, intervals=intervals)
                # breakpoint()

            if int(key[0]) == 32 and int(key[1]) == 144:
                assert c32 is None
                c32 = exceltable["avg"]

                simulation_times, simulated_tracer_at_times = sim_at_mri_times(pat, mritimes=experimentaltable["t"] / 3600, simdata=exceltable["avg"], simtimes=exceltable["t"] / 3600)    

                
                _, c32intervals = get_data_in_intervals(pat="205", stored_times=simulation_times, stored_data=simulated_tracer_at_times, intervals=intervals)
            if int(key[0]) == 64 and int(key[1]) == 144:
                assert c64 is None
                c64 = exceltable["avg"]

                simulation_times, simulated_tracer_at_times = sim_at_mri_times(pat, mritimes=experimentaltable["t"] / 3600, simdata=exceltable["avg"], simtimes=exceltable["t"] / 3600)    

                
                _, c64intervals = get_data_in_intervals(pat="205", stored_times=simulation_times, stored_data=simulated_tracer_at_times, intervals=intervals)

            if int(key[0]) == 64:
                c48_list64.append(c64intervals[1])
            if int(key[0]) == 32:
                c48_list32.append(c32intervals[1])

            if int(key[1]) == 96:

                label = "$h_{\mathrm{max}}=" + format(vertices_cells_hmax[int(key[0])][2], ".0f") + "$ mm" # + " " + str(dtlookup[key[1]])
                plt.plot(experimentaltable["t"] / 3600, v * experimentaltable[roi], marker="x", 
                        label="data", 
                        linewidth=0,
                        color=colors[int(key[0])],
                        )
            else:
                label = None

            # if not omit_time:
            #     label += ", dt = " + format(int(data["dt"]) / 60, ".0f") + " min"

            # if not idr == 0:
            #     label = None

            


            plt.plot(exceltable["t"] / 3600, v * exceltable[roi], 
                color=colors[int(key[0])], 
                label=label,
                linestyle=linestyles[key[1]],
                )

    assert c16.shape == c64.shape
    assert c32.shape == c64.shape

    p = 1

    print("rel. difference between 16 and 64", np.mean(np.abs(c16-c64)**p) / np.mean(np.abs(c64)**p))            
    print("rel. difference between 32 and 64", np.mean(np.abs(c32-c64)*p) / np.mean(np.abs(c64)**p))
    print("Ration c32/c64 at data points", c32intervals / c64intervals)

    print("Ration c16/c64 at data points", c16intervals / c64intervals)

    l2diff48 = ((np.array(c48_list32)-np.array(c48_list64))**2) / (np.array(c48_list64)**2)
    print("Maximum rel. difference at 48", np.max(l2diff48))
    # plt.title("Simulated gadobutrol concentration", fontsize=fs)
    plt.xlabel("Time (hours)", fontsize=fs)
    plt.ylabel("Brain-wide (mmol / L)", fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(fontsize=fs-4)
    plt.tight_layout()

    plt.savefig("/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/Brain influx and clearance during sleep and sleep deprivation/figures/mesh_resolution" + pat + ".png", dpi=444)

    # plt.show()