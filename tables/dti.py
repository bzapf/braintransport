import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
import json
import os

"""
Print average numbers for DTI
"""





path_to_files = "/home/basti/Dropbox (UiO)/Sleep/"


if not os.path.isfile("/home/basti/Dropbox (UiO)/Sleep/dti_values.json"):
    from dolfin import *
    from definitions import groups

    group = "dti"
    pats = groups[group]

    data_dict = {}

    for idx, pat in enumerate(sorted(pats)):

        dt = 1800
        res = 32
        if pat == 241:
            res = 64

        mesh = Mesh()

        hdf = HDF5File(mesh.mpi_comm(), "%s/%.3d/mesh/parenchyma%d_with_DTI.h5" %
                (path_to_files, int(pat), res), "r")

        hdf.read(mesh, "/mesh", False)
        SD = MeshFunction("size_t", mesh, mesh.topology().dim())
        hdf.read(SD, "/subdomains")
        
        bnd = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
        hdf.read(bnd, "/boundaries")
        lookup_table = MeshFunction("size_t", mesh, mesh.topology().dim())
        hdf.read(lookup_table, '/lookup_table')
        # lookup_boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
        # hdf.read(lookup_boundary, "/lookup_boundary")

        TensorSpace = TensorFunctionSpace(mesh, 'DG', 0)
        MDSpace = FunctionSpace(mesh, 'DG', 0)
        MD = Function(MDSpace)
        Kt = Function(TensorSpace)
        hdf.read(MD, '/MD')
        hdf.read(Kt, '/DTI')

        V = FunctionSpace(mesh, "CG", 1)

        write_regions = np.unique(lookup_table.array()[:])

        # GRAY = 1. WHITE = 2. BRAIN STEM = 3.
        dx_SD = Measure('dx')(domain=mesh, subdomain_data=SD)
        dx_lookup = Measure('dx')(domain=mesh, subdomain_data=lookup_table)

        mean_d_gray = assemble(MD * dx_SD(1)) / assemble(1 * dx_SD(1))

        mean_d_stem = assemble(MD * dx_SD(3)) / assemble(1 * dx_SD(3))

        mean_d_white = assemble(tr(Kt) * dx_SD(2)) / (3 * assemble(1 * dx_SD(2)))


        data_dict[pat] = [mean_d_gray, mean_d_white, mean_d_stem]

        print(idx, "/", len(pats), "pat:", pat, [mean_d_gray, mean_d_white, mean_d_stem])

        with open("/home/basti/Dropbox (UiO)/Sleep/dti_values.json", 'w') as outfile:
            json.dump(data_dict, outfile, sort_keys=True, indent=4)


with open("/home/basti/Dropbox (UiO)/Sleep/dti_values.json") as data_file:    
    data = json.load(data_file)

d1, d2, d3 = 0,0,0

print("Mean diffusion coefficient from DTI (for water!)")
for key, item in data.items():
    print("pat:", key, item)

    d1 += item[0]
    d2 += item[1]
    d3 += item[2]

MD_water = 3e-3
MD_gadovist = 3.8e-4 # (taken from Gd-DPTA)
scale_diffusion_gad = MD_gadovist / MD_water
D_scale = scale_diffusion_gad

print("Diffusivities for gadobutrol:")
print("group average mean D in gray", d1 * D_scale / len(data))
print("group average mean D in white", d2 * D_scale/ len(data))
print("group average mean D in brainstem", d3* D_scale/ len(data))
