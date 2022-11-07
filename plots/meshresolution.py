from fenics import *
import os
import matplotlib.pyplot as plt
import numpy as np
import json
# from IPython import embed
import pandas as pd
from matplotlib.pyplot import cm
import seaborn
import scipy.integrate as integrate
from scipy.interpolate import interp1d


pat = 205
pat = 191
paths = ["/home/basti/Dropbox (UiO)/Sleep/" + str(pat) + "/results_alpha_dt/"]

pat = 235
paths = ["/home/basti/Dropbox (UiO)/Sleep/" + str(pat) + "/results_alpha_dt/",
        "/home/basti/Dropbox (UiO)/Sleep/" + str(pat) + "/results_alpha/",
]

omit_time = False
mesh32_only = True

alpha = 2
rd = 5e-5
rn = 5e-5


alpha = 4
rd = 1e-4
rn = 1e-4
#rd = 0
#rn = 0


resfolders = {}
dts = set()

settings = set()

for path in paths:
    for wp in os.listdir(path):

        wp = os.path.join(path, wp)
        # breakpoint()

        for paramfolder in os.listdir(wp):
            
            resfolder = os.path.join(wp, paramfolder)
            assert os.path.isdir(resfolder)

            # breakpoint()
            
            resfolder = os.path.join(paramfolder, resfolder)

            with open(os.path.join(resfolder, "parserarguments")) as f:
                data = json.load(f)

            if mesh32_only and str(data["resolution"]) != "32":
                continue

            settings.add((data["dt"], data["resolution"] == "32", data["alpha"], data["gm_day_clearance"], data["gm_night_clearance"]))

            if not data["alpha"] == alpha:
                continue
            if not data["gm_night_clearance"] == rn:
                continue
            if not data["gm_day_clearance"] == rd:
                continue
            
            if int(data["dt"]) > 3200 and omit_time:
                continue

            if int(data["dt"]) < 600 and omit_time:
                continue
            if omit_time:
                if int(data["dt"]) != 1800:
                    continue

            resfolders[(data["dt"], int(data["resolution"]))] = resfolder
            dts.add(data["dt"])

        for se in settings:
            print(se)
        
if pat == 191:
    volumes_avg = {"16": 1258992.328196706, "32": 1250652.68007822, "64": 1254042.2486467457}
    volumes_white = {"16": 516430, "32": 514565, "64": 512466}
    volumes_gray = {"16": 724921, "32": 718369, "64": 723875}
    vertices_cells_hmax = {16: (86046, 443576, 99999999), 32: (257001, 1216305, 99999999), 64: (973581, 4554988, 99999999)}
elif pat == 205:
    volumes_avg = {"16": 1214932, "32": 1207988, "64": 1209862.}
    volumes_white = {"16": 565397.7643429264, "32": 563620.135340711, "64": 562008.1748970981}
    volumes_gray = {"16": 629790.6226778422, "32": 624661.7064872711, "64": 628180.7425056563}
    vertices_cells_hmax = {16: (82439, 421944, 10.92), 32: (247904, 1179462, 5.47), 64: (950203, 4453739, 2.75)}

elif pat == 235:
    volumes_avg = {"32": 1, "64": 1209862.}
    volumes_white = {"16": 565397.7643429264, "32": 1, "64": 562008.1748970981}
    volumes_gray = {"16": 629790.6226778422, "32": 1, "64": 628180.7425056563}
    vertices_cells_hmax = {16: (82439, 421944, 10.92), 32: (1, 1, 42), 64: (950203, 4453739, 2.75)}

volumes = {"avg": volumes_avg, "white": volumes_white, "gray": volumes_gray}

assert len(resfolders) > 0

# for meshfile in ["parenchyma16.h5", "parenchyma32.h5", "parenchyma64_with_DTI.h5", "parenchyma64.h5"]:
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

colors_ = list(cm.brg(np.linspace(0, 1, 3)))
# colors_ = seaborn.color_palette("crest", 3)

colors = {}
for idx, key in enumerate(["16", "32", "64"]):
    colors[int(key)] = colors_[idx]

colors = {16: "orangered", 32: "forestgreen", 64: "mediumblue"}

plt.figure(dpi=222)
fs = 15
lw = 1.5

# for idr, roi in enumerate(["avg", "white", "gray"]):
for idr, roi in enumerate(["avg"]):
    for idx, key in enumerate(sorted(resfolders.keys(), key=lambda x: x[1])):

        # breakpoint()
        resfolder = resfolders[key]

        # Plot concentration data as well
        # raise NotImplementedError("check which mesh resolution this is and scale accordingly")
        # if idx == 0:
        #     exceltable = pd.read_csv(os.path.join(resfolder, "experimental_data.csv"))
        #     plt.plot(exceltable["t"] / 3600, exceltable["avg"], marker="x", color="k", label="data", linewidth=0)

        # roi = "avg"

        exceltable = pd.read_csv(os.path.join(resfolder, "concs.csv"))
        experimentaltable = pd.read_csv(os.path.join(resfolder, "experimental_data.csv"))
        with open(os.path.join(resfolder, "parserarguments")) as f:
            data = json.load(f)

        c = "k"
        v = -1

        if "64" in resfolder:
            linestyle = "-"
            v = volumes[roi]["64"] * 1e-6
        elif "32" in resfolder:
            v = volumes[roi]["32"] * 1e-6
            linestyle = "--"
        elif "16" in resfolder:
            linestyle = ":"
            v = volumes[roi]["16"] * 1e-6

        linestyle = "-"
        assert str(data["resolution"]) in resfolder

        # f = interp1d(experimentaltable["t"], experimentaltable["avg_ds"])
        # data_avg = integrate.quad(f, [experimentaltable["t"][0], experimentaltable["t"][-1]])

        # v = v / data_avg

        # label = format(vertices_cells_hmax[int(data["resolution"])][0], ".0e")[0] + r"$\times 10^" + format(vertices_cells_hmax[int(data["resolution"])][0], ".0e")[-1] + "$"
        # label += " cells"

        label = "$h_{\mathrm{max}}=" + format(vertices_cells_hmax[int(data["resolution"])][2], ".0f") + "$ mm"

        if not omit_time:
            label += ", dt = " + format(int(data["dt"]) / 60, ".0f") + " min"

        if not idr == 0:
            label = None

        colors = {200: "red", 600: "green", 1800: "blue", 3200: "yellow", 5400: "k"}


        plt.plot(exceltable["t"] / 3600, v * exceltable[roi], 
            # color=colors[int(data["resolution"])], 
            color=colors[int(data["dt"])], 
            # linewidth=lw,
            # linewidth=lw * float(data["dt"]) / 1800    ,
            label=str(data["resolution"]) + " " + str(data["dt"]), 
            # label=label,
            linestyle=linestyle
            )
        
        # plt.plot(experimentaltable["t"] / 3600, v * experimentaltable[roi], marker="x", color=colors[int(data["resolution"])], linewidth=0.,)
            # "label=str(data["resolution"]) + " " + str(data["dt"]), linestyle=linestyle)

# bg = (1, 1, 1, 0.5)
#transparency = 1
plt.text(25, 0.0205, "white", fontsize=fs - 2, 
        #backgroundcolor=bg, alpha=transparency)
        )
plt.text(30, 0.04, "gray", fontsize=fs - 2, rotation=-20,
    # backgroundcolor=bg, alpha=transparency,
    )

if fs == 16:
    yavg = 0.045
else:
    yavg = 0.047
plt.text(37, yavg, "brain average", fontsize=fs - 2, rotation=-22,
    # backgroundcolor=bg, alpha=transparency,
    )

plt.xlabel("Time (hours)", fontsize=fs)
plt.ylabel("Tracer concentration (mmol)", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
# for dt in dts:
#     plt.plot([], [], linestyle="-", label="64" + " " + str(dt), color=colors[dt])
#     plt.plot([], [], linestyle="--", label="32" + " " + str(dt), color=colors[dt])
plt.legend(fontsize=fs)
plt.tight_layout()

# plt.savefig("/home/basti/Dropbox (UiO)/Apps/Overleaf (1)/Brain influx and clearance during sleep and sleep deprivation/figures/mesh_resolution.pdf")

plt.show()

# embed()