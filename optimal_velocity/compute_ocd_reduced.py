import sys
import os
import glob
import shutil
import resource
import numpy
import math
import argparse
import time
import datetime
import csv

from dolfin import *

from compute_omt import map_image_onto_mesh
from compute_ocdr import non_negativize, diffusion
# from postprocess_omt import csvwrite, compute_magnitude_field
from postprocess_phi import csvwrite, compute_magnitude_field


from dolfin_adjoint import *

def bar(phi, volume=None):
    phi_mag = compute_magnitude_field(phi)
    if not volume:
        mesh = phi.function_space().mesh()
        volume = assemble(1.0*dx(domain=mesh))
    avg_phi_mag = assemble(phi_mag*dx())/volume
    return avg_phi_mag
    
def compute_ocd_reduced(c0, c1, tau, D, alpha, results_dir, space="CG", reg="H1"):
    # c0:
    # c1:
    # tau:   time step
    # D:     diffusion coefficient
    # alpha: regularization parameter
    # space: finite element space for the velocity field ("CG" | "RT" | "BDM")
    # reg:   if "H1", use H1-regularization, else use H(div)
    
    info("Computing OCD via reduced approach")

    # Define mesh and function space for the concentration
    mesh = c0.function_space().mesh()
    C = FunctionSpace(mesh, "CG", 1)
    
    # Space for the convective velocity field phi
    if space == "CG":
        Q = VectorFunctionSpace(mesh, "CG", 1)
    else:
        Q = FunctionSpace(mesh, space, 1)
    phi = Function(Q, name="Control")

    # Regularization term
    def R(phi, alpha, mesh):
        if reg == "H1":
            form = 0.5*alpha*(inner(phi, phi) + inner(grad(phi), grad(phi)))*dx(domain=mesh)
        else:
            form = 0.5*alpha*(inner(phi, phi) + inner(div(phi), div(phi)))*dx(domain=mesh)
        return form

    # Define previous solution
    c_ = Function(C)
    c_.assign(c0) # Hack to make dolfin-adjoint happy, maybe just start tape here?

    c2 = Function(C)
    c2.assign(c1)
    
    # Define variational problem
    c = TrialFunction(C)
    d = TestFunction(C)
    F = (1.0/tau*(c - c_)*d + div(c*phi)*d)*dx() + diffusion(D, c, d)
    a, L = system(F)
    bc = DirichletBC(C, c2, "on_boundary")

    # ... and solve it once
    c = Function(C, name="State")
    solve(a == L, c, bc, solver_parameters={"linear_solver": "mumps"})

    # Output max values of target and current solution for progress
    # purposes
    info("\max c_1 = %f" % c2.vector().max())
    info("\max c = %f" % c.vector().max())

    # Define the objective functional
    j = 0.5*(c - c2)**2*dx(domain=mesh) + R(phi, alpha, mesh)
    J = assemble(j)
    info("J (initial) = %f" % J)
    
    # Define control field
    m = Control(phi)

    # Define call-back for output at each iteration of the optimization algorithm
    name = lambda s: os.path.join(results_dir, "opts", s)
    dirname = os.path.join(results_dir, "opts")
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    header = ("j", "\max \phi")

    def eval_cb(j, phi):
        values = (j, phi.vector().max())
        mem = resource.getrusage(resource.RUSAGE_SELF)[2]
        info("Current memory usage: %g (MB)" % (mem/1024))
        info("\tj = %f, \max phi = %f (mm/h)" % values)
        csvwrite(name("optimization_values.csv"), values, header, mode="a")

        # Read the optimization counter file, update counter, and
        # write it back, geez. 
        counter = 0
        with open(os.path.join(results_dir, "counter.csv"), "r") as f:
            reader = csv.reader(f)
            for row in reader:
                counter = int(row[0])
        counter += 1

        with open(os.path.join(results_dir, "counter.csv"), "w") as f:
            info("Updating counter file, counter is now %d " % counter)
            writer = csv.writer(f)
            writer.writerow((counter,))     

        # Write current control variable to file in HDF5 and PVD formats
        file = HDF5File(mesh.mpi_comm(), name("opt_phi_%d.h5" % counter), "w")
        file.write(phi, "/function", 0)
        file.close()

    # Define reduced functional in terms of J and m
    Jhat = ReducedFunctional(J, m, eval_cb_post=eval_cb)

    # Minimize functional
    tol = 1.0e-8
    phi_opt = minimize(Jhat,
                       tol=tol, 
                       options={"gtol": tol, "maxiter": 80, "disp": True})
    pause_annotation()

    # Update phi, and do a final solve to compute c
    phi.assign(phi_opt)
    solve(a == L, c, bc, solver_parameters={"linear_solver": "mumps"})

    J = assemble(j)
    j0 = 0.5*(c - c2)**2*dx(domain=mesh)
    jr = R(phi, alpha, mesh)
    J0 = assemble(j0)
    Jr = assemble(jr)
    info("J  = %f" % J)
    info("J0 = %f" % J0)
    info("Jr = %f" % Jr)
    
    return (c, phi) 

def find_image_pair(images, key):

    # Read dates/times of all images
    times = []
    for image in images:
        date, time, _ = os.path.basename(image).split("_")
        d = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
        t = datetime.time(int(time[0:2]), int(time[2:4]), int(time[4:6]))
        times += [datetime.datetime.combine(d, t)]

    # Set all times to the default (the smallest time point, should be 0.0)
    t0 = min(times)
    t1 = min(times)
    t2 = min(times)

    # "0h" refers to pair of first image at 0h and image after 4.9h - 6.7h ('6h')
    if key == "0h":
        dt0 = datetime.timedelta(hours=0)
        dt1 = datetime.timedelta(hours=4.9)
        dt2 = datetime.timedelta(hours=6.7)

        # Set t1 to be the first time and index
        t1 = min(times)
        i1 = 0

        # Find the second time and index
        for (i, t) in enumerate(times):
            dt = t - t0 
            if (dt1 < dt < dt2) and t > t2: 
                t2 = t
                i2 = i

    # "6h" refers to pair of image after 4.9h - 6.7h ('6h') and image 8-30 hours
    elif key == "6h":
        # For day after to day day after 24h to 48h    
        dt0 = datetime.timedelta(hours=4.9)
        dt1 = datetime.timedelta(hours=6.7)
        dt2 = datetime.timedelta(hours=26)
        
        # Extract pair of images between dt0 and dt1, and dt1 and dt1
        for (i, t) in enumerate(times):
            dt = t - t0
            if dt0 < dt < dt1:
                t1 = t
                i1 = i
            elif dt1 < dt < dt2:
                t2 = t
                i2 = i
            else:
                pass

    elif key == "24h":
        # For day after to day day after 24h to 48h    
        dt0 = datetime.timedelta(hours=7)
        dt1 = datetime.timedelta(hours=30)
        dt2 = datetime.timedelta(hours=60)
        
        # Extract pair of images between dt0 and dt1, and dt1 and dt1
        for (i, t) in enumerate(times):
            dt = t - t0
            if dt0 < dt < dt1:
                t1 = t
                i1 = i
            elif dt1 < dt < dt2:
                t2 = t
                i2 = i
            else:
                pass

    else:
        info("Unrecognized image pair starting at %s, exiting" % key)
        exit()

    if (i2 == i1):
        info("Unable to identify pair of images, indices are: %d, %d" % (i1, i2))
        exit()

    info("Images at %s and %s with dt = tau (h) %s\n" % (str(t1), str(t2), str(t2-t1)))

    # Return image filenames
    return (images[i1], images[i2], t2-t1)


def read_image_pair(patient_dir, mesh, key, imfolder):

    # "0h" refers to pair of first image at 0h and image after 4.9h - 6.7h ('6h')
    # "6h" refers to pair of image after '6h' and image after 24h
    # "24h" refers to pair of image after 24h and image after 48h
    
    # Read set of images for this patient
    image_dir = os.path.join(patient_dir, imfolder) #"FIGURES_CONC")
    if not os.path.isdir(image_dir):
        info("Missing image directory %s, exiting" % image_dir)
        exit()
    images = glob.glob(os.path.join(image_dir, "*.mgz"))
    images = sorted(images)
    # Find the contrast images we want
    i0, i1, dt = find_image_pair(images, key)
    tau = dt.seconds/(60*60) # Compute rate in hours

    # Map images onto the mesh
    c0 = map_image_onto_mesh(i0, mesh)
    c1 = map_image_onto_mesh(i1, mesh)

    # Adjust negative values if any
    c0 = non_negativize(c0)
    c1 = non_negativize(c1)

    info("||c_0||_0 = %g" % norm(c0, "L2"))
    info("||c_1||_0 = %g" % norm(c1, "L2"))

    return (c0, c1, tau)
        
def map_patient(patient_dir, n, alpha, key, results_dir, use_dti, imfolder):
    """
    patient_dir: Absolute patient directory e.g. "../../241"
    n:           Mesh resolution e.g. 16 (depending on available meshes)
    alpha:       Regularization parameter
    key:         0h, 6h or 24h (starting time)
    results_dir: Name for results directory (subdirectories will be created)
    """

    # Read FEniCS mesh 
    mesh = Mesh()
    

    meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d.h5" % n)
    
    if (not os.path.isfile(meshfile)) or use_dti:
        meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d_with_DTI.h5" % n)
    
    if not os.path.isfile(meshfile):
        raise Exception("Missing meshfile %s" % meshfile)
    
    info("Reading mesh from %s" % meshfile)
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)

    info("Reading subdomains (gray 1, white 2, brainstem 3) from %s" % meshfile)
    subdomains = MeshFunction("size_t", mesh,mesh.topology().dim())
    hdf.read(subdomains, "/subdomains")

    # Also read DTI information (mean diffusivity MD, and DTI K) from meshfile
    # info("Reading MD and DTI from %s" % meshfile)
    Q = FunctionSpace(mesh, "DG", 0)

    if use_dti:
        MD = Function(Q)
        hdf.read(MD, "/MD")

        W = TensorFunctionSpace(mesh, 'DG', 0)
        K = Function(W)
        hdf.read(K, "/DTI")

        print("Succesfully read DTI from mesh", meshfile)
    
    hdf.close()
    info("... with %d vertices and %d cells" % (mesh.num_vertices(), mesh.num_cells()))

    # Construct the region-wise diffusion tensor
    D_H20 = 0.001049 # mm^2
    D_gad = 0.000130 # mm^2
    Dstar = D_gad/D_H20

    if use_dti:
        D1 = Dstar * MD    # Diffusion coefficient in gray matter (subdomain 1)
        D3 = D1            # Diffusion coefficient in brain stem (subdomain 3)
        D2 = Dstar * K     # Diffusion tensor in white matter (subdomain 2)
    else:

        print("Not using DTI, using average D instead")

        D1 = Function(Q)
        D2 = Function(Q)
        D3 = Function(Q)

        D1.vector()[:] = Dstar * 0.00092582074369026 # Diffusion coefficient in gray matter (subdomain 1)
        D2.vector()[:] = Dstar * 0.000867643624860488 # Diffusion tensor in white matter (subdomain 2)
        D3.vector()[:] = Dstar *  0.0010602528519216016 # Diffusion coefficient in brain stem (subdomain 3)



    D = ((D1, D2, D3), subdomains)

    # Read set of images for this patient
    (c0, c1, tau) = read_image_pair(patient_dir, mesh, key, imfolder)
    
    # Store images and solution to .h5 format for postprocessing 
    name = lambda s: os.path.join(results_dir, "hdf5", s)
    info("Storing images and data to %s" % name("..."))
    file = HDF5File(mesh.mpi_comm(), name("c0.h5"), "w")
    file.write(c0, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c1.h5"), "w")
    file.write(c1, "/function", 0)
    file.close()

    # Make an optimization counter file
    with open(os.path.join(results_dir, "counter.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow((0,))     
            
    # Compute OCD approximation via reduced method
    space = "CG"
    reg = "H1"
    (c, phi) = compute_ocd_reduced(c0, c1, tau, D, alpha, results_dir,
                                   space=space, reg=reg)

    # -----------------------------------------------------------
    # The following is boiler-plate post-processing code
    # -----------------------------------------------------------
    file = HDF5File(mesh.mpi_comm(), name("phi.h5"), "w")
    file.write(phi, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c.h5"), "w")
    file.write(c, "/function", 0)
    file.close()

    # Copy the mesh there as well
    if "with_DTI" in meshfile:
        shutil.copy(meshfile, name("parenchyma%d_with_DTI.h5" % n))
    else:
        shutil.copy(meshfile, name("parenchyma%d.h5" % n))
    

    # Compute and output some key numbers
    phi_L2 = norm(phi, "L2")
    phi_Hdiv0 = norm(phi, "Hdiv0")
    omega = assemble(1*dx(domain=mesh))
    avg_phi = bar(phi, omega)
    c_L2 = norm(c, "L2")
    c1_L2 = norm(c1, "L2")

    values = (phi_L2, avg_phi, phi_Hdiv0, c_L2, c1_L2)
    header = ("|phi|_L2", "phi_avg (mm/h)", "|div phi|_L2" , "|c|_L2", "|c1|_L2")
    name = lambda s: os.path.join(results_dir, s)
    csvwrite(name("values.csv"), values, header)

    # Store images and phi to .pvd format as well
    name = lambda s: os.path.join(results_dir, "pvd", s)
    info("Storing images and data to %s" % name("..."))
    file = File(name("c_e.pvd"))
    file << c0
    file << c1

    file = File(name("c.pvd"))
    file << c

    file = File(name("phi.pvd"))
    if space != "CG":
        V = VectorFunctionSpace(mesh, "CG", 1)
        phi_CG = project(phi, V)
        file << phi_CG
    else:
        file << phi

    info("Reduced OCD computed successfully.\n")
    return (c, phi)

def run_ocd_reduced(num=None, n=None, beta=None, key=None):
    # Take patient number, mesh resolution parameter n (16, 32, 64),
    # beta > 0 and night/day as input

    # Handle arguments
    t0 = time.time()
    set_working_tape(Tape())

    if not num:
        num = sys.argv[1]
        n = int(sys.argv[2])
        beta = float(sys.argv[3])
        key = sys.argv[4]
        use_dti = sys.argv[5]
        assert use_dti in ["true", "false"]
        use_dti = use_dti == "true"
        assert isinstance(use_dti, bool)
        imfolder = sys.argv[6]

    info("Handling patient %s with mesh %d, beta=%.1e at %s" % (num, n, beta, key))

    # Use ./num as input directory and set output ./key as output
    # directory
    path = os.getcwd()
    patient = os.path.join(path, num)
    dirname = "results_red_ocd_pat%s_n%d_beta%.1e_%s" % (num, n, beta, key)
    output_dir = os.path.join(patient, dirname)
    info("Setting %s as output directory" % output_dir)

    # Create OCD(R) map
    (c, phi) = map_patient(patient, n, beta, key, output_dir, use_dti, imfolder)

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t1 = time.time()

    info("\nMax memory usage (MB) = %g\n" % (mem/1024))
    info("\nTotal time (min) = %g\n" % ((t1 - t0)/60))

if __name__ == "__main__":

    # python3 compute_ocd_reduced.py 240 32 1.0e-04 night
    run_ocd_reduced()

