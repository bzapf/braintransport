import resource
import sys
import os
import os.path
import numpy
import csv
from IPython import embed

from dolfin import *
#from mpi4py import MPI

def read_mesh(meshfile):
    mesh = Mesh()
    if not os.path.isfile(meshfile):
        raise Exception("Missing meshfile %s" % meshfile)
    info("Reading mesh, subdomains and look-up table from %s" % meshfile)
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)
    subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(subdomains, "/subdomains") # Gray, white, brainstem
    lookup_table = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(lookup_table, "/lookup_table") # All FreeSurfer regions
    hdf.close()
    info("... with %d vertices and %d cells" % (mesh.num_vertices(), mesh.num_cells()))
    return (mesh, subdomains, lookup_table)

def read_field(mesh, phifile, V=None, phi=None):
    if V is None:
        V = VectorFunctionSpace(mesh, "CG", 1)
    if phi is None:
        phi = Function(V)
    hdf = HDF5File(mesh.mpi_comm(), phifile, "r")
    hdf.read(phi, "/function")
    hdf.close()
    return phi

def compute_averages(v, mesh, subdomains, unit=None):

    # Compute total average
    dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
    one = Constant(1.0)
    V = assemble(one*dx(domain=mesh))
    avg_b = assemble(v*dx())/V
    
    # Compute volumes of the gray matter (V_g), white matter (V_w) and
    # brain stem (V_s) regions
    V_g = assemble(one*dx(1, domain=mesh))
    V_w = assemble(one*dx(2, domain=mesh))
    V_s = assemble(one*dx(3, domain=mesh))

    # Compute average of v over given regions
    avg_g = assemble(v*dx(1))/V_g
    avg_w = assemble(v*dx(2))/V_w
    avg_s = assemble(v*dx(3))/V_s

    averages = (avg_b, avg_g, avg_w, avg_s)
    if not unit:
        unit = "mm/h"
    
    headers = ("avg (%s)" % unit, "avg_gray", "avg_white", "avg_stem") 
    return (averages, headers)

def compute_magnitude_field(phi):

    # Split phi into its components (really split it)
    (phi0, phi1, phi2) = phi.split(deepcopy=True)

    # Get the map from vertex index to dof
    v2d = [vertex_to_dof_map(phii.function_space()) for phii in (phi0, phi1, phi2)]

    # Create array of component values at vertex indices so p0[0] is
    # value at vertex 0 e.g. 
    p0 = phi0.vector().get_local(v2d[0])
    p1 = phi1.vector().get_local(v2d[1])
    p2 = phi2.vector().get_local(v2d[2])

    # Take element-wise magnitude
    m = numpy.sqrt((numpy.square(p0) + numpy.square(p1) + numpy.square(p2)))

    # Map array of values at vertices back to values at dofs and
    # insert into magnitude field.
    M = phi0.function_space()
    magnitude = Function(M)
    d2v = dof_to_vertex_map(M)
    magnitude.vector()[:] = m[d2v]
    return magnitude

def compute_vorticity(phi):
    V = phi.function_space()
    curl_phi = Function(V)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u, v)*dx()
    L = inner(curl(phi), v)*dx()
    A = assemble(a)
    b = assemble(L)
    solve(A, curl_phi.vector(), b, "mumps")
    return curl_phi

def compute_divergence(phi):
    Q = phi.function_space().sub(0).collapse()
    div_phi = Function(Q)
    u = TrialFunction(Q)
    v = TestFunction(Q)
    a = u*v*dx()
    L = div(phi)*v*dx()
    A = assemble(a)
    b = assemble(L)
    solve(A, div_phi.vector(), b, "mumps")
    return div_phi

def extract_max_mins(u, eps=1.0e-12, unit=None):
    values = u.vector().get_local()
    non_zeros = values[(abs(values) > eps)]
    u_max = values.max()
    u_min = values.min()
    nz_min = non_zeros.min()

    maxmins = (u_max, u_min, nz_min)
    if not unit:
        unit = "mm/h"
    headers = ("max (%s)" % unit, "min", "min > eps, eps=%.3e" % eps )
    
    return (maxmins, headers)

def hdf5write(name, field):
    mesh = field.function_space().mesh()
    file = HDF5File(mesh.mpi_comm(), name, "w")
    file.write(field, "/function", 0)
    file.close()

def hdf5write_mf(name, mf, mesh):
    hfile = HDF5File(mesh.mpi_comm(), name, "w")
    hfile.write(mf, "/meshfunction")
    hfile.close()

def csvwrite(name, values, header, mode=None, debug=True): 
   # Default to write, set mode to "a" for append
   if not mode:
       mode = "w"
   # If mode is append, but the file does not exist (yet), write instead
   elif mode == "a" and not os.path.isfile(name):
       mode = "w"
   if debug:
       info(name)
       info(" ".join(str(s) for s in header))
       info(" ".join("%.3e" % g for g in values))
   with open(name, mode) as f:
       writer = csv.writer(f)
       writer.writerow(header)     
       writer.writerow(values)     

def compute_regional_averages(u, mesh, lookup_table, eps=1.0e-12):
       
    # Compute average velocity magnitude for each FreeSurfer region
    dy = Measure("dx", domain=mesh, subdomain_data=lookup_table)
    regions = numpy.unique(lookup_table.array())
    regions.sort()
    avgs = {}
    one = Constant(1.0)
    for r in regions:
        vol = assemble(one*dy(int(r), domain=mesh))  # Compute volume of region r
        if (abs(vol) < eps):
            avgs[r] = 0.0
        else:
            avgs[r] = assemble(u*dy(int(r)))/vol   # Compute average |phi| over region r

    # Map the regional averages into a mesh function
    u_avg = MeshFunction("double", mesh, mesh.topology().dim()) 
    cell2region = lookup_table.array()
    for c in range(u_avg.array().size):
        u_avg[c] = avgs[cell2region[c]]

    return (regions, avgs, u_avg)
       
def postprocess_fields(datadir, n):

    # Utility lambda function for writing datadir/f:
    output = lambda f: os.path.join(datadir, f)

    # Read the mesh
    timer = Timer("XXX: Read mesh")
    meshfile = output("hdf5/parenchyma%d_with_DTI.h5" % n) 
    if not os.path.isfile(meshfile):
        meshfile = output("hdf5/parenchyma%d.h5" % n)
    (mesh, subdomains, lookup_table) = read_mesh(meshfile)
    timer.stop()

    # Read phi
    timer = Timer("XXX: Read phi.h5")
    phifile = output("hdf5/phi.h5")
    phi = read_field(mesh, phifile)
    timer.stop()

    # Compute phi magnitude field and store as .h5 and .pvd
    timer = Timer("XXX: Compute magnitude of phi")
    info("Computing the magnitude (abs) of phi")
    phi_mag = compute_magnitude_field(phi)
    # hdf5write(output("hdf5/phi_mag.h5"), phi_mag)
    # phi_mag_file = File(output("pvd/phi_mag.pvd"))
    # phi_mag_file << phi_mag
    timer.stop()

    # Compute average phi magnitudes:
    averages, headers = compute_averages(phi_mag, mesh, subdomains)
    avgsfile = output("phi_mag_avgs.csv")
    csvwrite(avgsfile, averages, headers)

    # Compute regional averages and store
    # (regions, avgs, phi_mag_avg) = compute_regional_averages(phi_mag, mesh, lookup_table)
    # hdf5write_mf(output("hdf5/phi_mag_reg_avgs.h5"), phi_mag_avg, mesh)
    # regfile = File(output("pvd/phi_mag_reg_avgs.pvd"))
    # regfile << phi_mag_avg
    
    # Compute max/min/min > 0 phis :
    maxmins, headers = extract_max_mins(phi_mag)
    maxfile = output("phi_mag_maxmins.csv")
    csvwrite(maxfile, maxmins, headers)
    
    # Compute and store div(phi)
    timer = Timer("XXX: Compute div(phi)")
    info("Computing the divergence (div) of phi")
    div_phi = compute_divergence(phi)
    
    # hdf5write(output("hdf5/div_phi.h5"), div_phi)
    # div_phi_file = File(output("pvd/div_phi.pvd"))
    # div_phi_file << div_phi
    timer.stop()

    # Compute average div(phi)
    averages, headers = compute_averages(div_phi, mesh, subdomains,
                                         unit="1/h")
    avgsfile = output("div_phi_avgs.csv")
    csvwrite(avgsfile, averages, headers)

    # Compute curl(phi) and store it
    # timer = Timer("XXX: Compute curl of phi")
    # info("Computing the vorticity (curl) of phi")
    # curl_phi = compute_vorticity(phi)
    # hdf5write(output("hdf5/curl_phi.h5"), curl_phi)
    # curl_phi_file = File(output("pvd/curl_phi.pvd"))
    # curl_phi_file << curl_phi
    # timer.stop()

    # Compute magnitude of curl(phi) and store it
    # timer = Timer("XXX: Compute magnitude of curl_phi")
    # info("Computing the magnitude of curl(phi)")
    # curl_phi_mag = compute_magnitude_field(curl_phi)
    # hdf5write(output("hdf5/curl_phi_magnitude.h5"), curl_phi_mag)
    # curl_phi_mag_file = File(output("pvd/curl_phi_mag.pvd"))
    # curl_phi_mag_file << curl_phi_mag
    # timer.stop()

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    info("\nMax memory usage (MB) = %g\n" % (mem/1024))

if __name__ == "__main__":

    ## For single runs on Saga,
    ## Run these commands first at Saga in the current directory
    # salloc --ntasks=1 --mem-per-cpu=10G --time=00:59:00 --account=NN9279K --qos=devel
    # source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

    # datadir e.g.: /cluster/projects/nn9279k/meg/modelling-sleep-deprivation-2021/results-red-ocd/002/simulations/results_red_ocd_pat002_n32_beta1.0e-04_6h
    datadir = sys.argv[1]
    n = int(sys.argv[2])
    info("Running from %s with n=%d" % (datadir, n))

    postprocess_fields(datadir, n)

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    info("\nMax memory usage (MB) = %g\n" % (mem/1024))
        

    info("Success!")
