import sys
import os
import glob
import shutil
import resource
import numpy
import math

from dolfin import *

from compute_omt import map_image_onto_mesh, find_image_pair
# from postprocess_omt import csvwrite
from postprocess_phi import csvwrite

def non_negativize(I):
    """Given a scalar Function c, return the field with negative values
    replaced by zero.

    """
    c = I.copy(deepcopy=True) 
    vec = c.vector().get_local()
    # Extract dof indices where the vector values are negative
    negs = numpy.where(vec < 0)[0]
    # Set these values to zero
    c.vector()[negs] = 0.0
    return c

def diffusion(D, d, y):
    if not hasattr(D, "__len__"):
        F = inner(D*grad(d), grad(y))*dx()
    else:
        (Ds, subdomains) = D
        mesh = subdomains.mesh()
        dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
        F = sum([inner(Di*grad(d), grad(y))*dx(i+1) for (i, Di) in enumerate(Ds)])
    return F
        
def compute_ocdr(c0, c1, tau, D, alpha):
    info("Computing OCDR")

    cell = triangle # tetrahedron
    C = FiniteElement("CG", cell, 1) # H1
    Y = FiniteElement("CG", cell, 1) # H10
    Q = VectorElement("CG", cell, 1) # H(div), so RT?
    R = FiniteElement("R", cell, 0)
    G = FiniteElement("R", cell, 0)
    Mx = MixedElement([C, Y, Q, R, G])
    
    mesh = c0.function_space().mesh()
    M = FunctionSpace(mesh, Mx)

    m = Function(M)
    (c, y, phi, r, g) = split(m)
    (d, z, psi, s, k) = TestFunctions(M)

    # Handle diffusion coefficient either as single value, or given in
    # terms of subdomains labelled starting at 1
    F1 = (c*d - c1*d + (1./tau*d + div(d*phi) + r*d)*y)*dx() + diffusion(D, d, y)
    F2 = (div(c*psi)*y + alpha*(inner(phi, psi) + div(phi)*div(psi)))*dx()
    F3 = ((1./tau*c + div(c*phi) + r*c - g)*z - 1./tau*c0*z)*dx() + diffusion(D, c, z)
    F4 = c*y*s*dx()
    F5 = y*k*dx()
    F = F1 + F2 + F3 + F4 + F5

    assign(m.sub(0), c1)

    bcs = [DirichletBC(M.sub(0), c1, "on_boundary"),
           DirichletBC(M.sub(1), 0.0, "on_boundary")]

    ps = {"newton_solver": {"linear_solver": "mumps", "maximum_iterations": 50}}
    solve(F == 0, m, bcs, solver_parameters=ps)
    (c, y, phi, r, g) = m.split(deepcopy=True)

    return (c, y, phi, r, g) 

def compute_ocd(c0, c1, tau, D, alpha, space="CG", auto=True, reg="H1"):

    # c0:
    # c1:
    # tau:   time step
    # D:     diffusion coefficient
    # alpha: regularization parameter
    # space: finite element space for the velocity field
    # auto:  if true, use auto-differentiation to define PDE
    # reg:   if "H1", use H1-regularization, else use H(div)
    
    info("Computing OCD")
    
    mesh = c0.function_space().mesh()
    cell = triangle if mesh.topology().dim() == 2 else tetrahedron

    C = FiniteElement("CG", cell, 1) # H1
    Y = FiniteElement("CG", cell, 1) # H10
    if space == "CG":
        Q = VectorElement("CG", cell, 1)
    else:
        Q = FiniteElement("RT", cell, 1)
    Mx = MixedElement([C, Y, Q])
    M = FunctionSpace(mesh, Mx)

    m = Function(M)
    (c, y, phi) = split(m)
    dm = TestFunction(M)
    (d, z, psi) = split(dm)

    assign(m.sub(0), c1)
    bcs = [DirichletBC(M.sub(0), c1, "on_boundary"),
           DirichletBC(M.sub(1), 0.0, "on_boundary")]

    def R(c, phi, alpha):
        if reg == "H1":
            form = 0.5*alpha*(inner(phi, phi) + inner(grad(phi), grad(phi)))*dx()
        else:
            form = 0.5*alpha*(inner(phi, phi) + inner(div(phi), div(phi)))*dx()
        return form
        
    if auto:
        L = 0.5*inner(c - c1, c - c1)*dx() \
            + R(c, phi, alpha) \
            + (1.0/tau*(c - c0) + div(c*phi))*y*dx() \
            + diffusion(D, c, y)
        F = derivative(L, m, dm)
        
    else:
        F1 = (c*d + (1./tau*d + div(d*phi))*y - c1*d)*dx() + diffusion(D, d, y)
        F2 = (div(c*psi)*y + alpha*(inner(phi, psi) + div(phi)*div(psi)))*dx()
        F3 = ((1./tau*c + div(c*phi))*z - 1./tau*c0*z)*dx() + diffusion(D, c, z)
        F = F1 + F2 + F3
        
    #ps = {"nonlinear_solver": "snes",
    #      "snes_solver": {"linear_solver": "mumps",
    #                      "maximum_iterations": 50,
    #                      "relative_tolerance": 1.e-6,
    #                      "line_search": "l2"}}
    ps = {"newton_solver": {"linear_solver": "mumps",
                            "maximum_iterations": 50,
                            "relaxation_parameter": 1.0}}
    solve(F == 0, m, bcs, solver_parameters=ps)
    (c, y, phi) = m.split(deepcopy=True)

    return (c, y, phi) 

def map_patient(patient_dir, n, alpha, night, results_dir):
    """
    patient_dir: Absolute patient directory e.g. "../../241"
    n:           Mesh resolution e.g. 16 (depending on available meshes)
    alpha:       Regularization parameter
    night:       True if comparison should be during night, false if during day 1.
    results_dir: Name for results directory (subdirectories will be created)
    """

    # Read FEniCS mesh 
    mesh = Mesh()
    meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d_with_DTI.h5" % n)
    if not os.path.isfile(meshfile):
        meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d.h5" % n)
        assert os.path.isfile(meshfile)
        # raise Exception("Missing meshfile %s" % meshfile)
    info("Reading mesh from %s" % meshfile)
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)

    info("Reading subdomains (gray 1, white 2, brainstem 3) from %s" % meshfile)
    subdomains = MeshFunction("size_t", mesh,mesh.topology().dim())
    hdf.read(subdomains, "/subdomains")

    # All FreeSurfer regions available in the stored mesh:
    #hdf.read(lookup_table, "/lookup_table", False)

    # Also read DTI information (mean diffusivity MD, and DTI K) from meshfile
    info("Reading MD and DTI from %s" % meshfile)
    Q = FunctionSpace(mesh, "DG", 0)
    MD = Function(Q)
    hdf.read(MD, "/MD")

    W = TensorFunctionSpace(mesh, 'DG', 0)
    K = Function(W)
    hdf.read(K, "/DTI")
    hdf.close()
    info("... with %d vertices and %d cells" % (mesh.num_vertices(), mesh.num_cells()))

    # Construct the region-wise diffusion tensor
    D_H20 = 0.001049 # mm^2
    D_gad = 0.000130 # mm^2
    Dstar = D_gad/D_H20
    D1 = Dstar*MD    # Diffusion coefficient in gray matter (subdomain 1)
    D3 = D1          # Diffusion coefficient in brain stem (subdomain 3)
    D2 = Dstar*K     # Diffusion tensor in white matter (subdomain 2)
    D = ((D1, D2, D3), subdomains)

    # Read set of images for this patient
    image_dir = os.path.join(patient_dir, "FIGURES_CONC_LUT")
    if not os.path.isdir(image_dir):
        info("Missing image directory %s, exiting" % image_dir)
        exit()
    images = glob.glob(os.path.join(image_dir, "*.mgz"))

    # Find the contrast images we want
    i0, i1, dt = find_image_pair(images, night)
    tau = dt.seconds/(60*60) # Compute rate in hours

    # Map images onto the mesh
    c0 = map_image_onto_mesh(i0, mesh)
    c1 = map_image_onto_mesh(i1, mesh)

    # Adjust negative values if any
    c0 = non_negativize(c0)
    c1 = non_negativize(c1)
    info("||c_0||_0 = %g" % norm(c0, "L2"))
    info("||c_1||_0 = %g" % norm(c1, "L2"))

    space = "CG"
    reg = "H1"
    auto = True
    (c, y, phi) = compute_ocd(c0, c1, tau, D, alpha,
                              space=space, auto=auto, reg=reg)

    # -----------------------------------------------------------
    # The following is boiler-plate post-processing code
    # -----------------------------------------------------------
    # Store images and solution to .h5 format for postprocessing
    name = lambda s: os.path.join(results_dir, "hdf5", s)
    info("Storing images and data to %s" % name("..."))
    file = HDF5File(mesh.mpi_comm(), name("c0.h5"), "w")
    file.write(c0, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c1.h5"), "w")
    file.write(c1, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("phi.h5"), "w")
    file.write(phi, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c.h5"), "w")
    file.write(c, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("y.h5"), "w")
    file.write(y, "/function", 0)
    file.close()

    # Copy the mesh there as well
    if "with_DTI" in meshfile:
        shutil.copy(meshfile, name("parenchyma%d_with_DTI.h5" % n))
    else:
        shutil.copy(meshfile, name("parenchyma%d.h5" % n))

    # Compute and output some key numbers
    phi_L2 = norm(phi, "L2")
    phi_Hdiv0 = norm(phi, "Hdiv0")
    omega = math.sqrt(assemble(1*dx(domain=mesh)))
    avg_phi = phi_L2/omega # mm/h
    c_L2 = norm(c, "L2")
    y_L2 = norm(y, "L2")
    c1_L2 = norm(c1, "L2")

    values = (phi_L2, avg_phi, phi_Hdiv0, c_L2, y_L2, c1_L2)
    header = ("|phi|_L2", "phi_avg (mm/h)", "|div phi|_L2" , "|c|_L2", "|y|_L2","|c1|_L2")
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
    if space == "RT":
        V = VectorFunctionSpace(mesh, "CG", 1)
        phi_CG = project(phi, V)
        file << phi_CG
    else:
        file << phi

    info("OCDR computed successfully.\n")
    return (c, phi, y)

def run_ocd():
    # Take patient number, mesh resolution parameter n (16, 32, 64),
    # beta > 0 and night/day as input
    # FIXME: argparse
    num = sys.argv[1]
    n = int(sys.argv[2])
    beta = float(sys.argv[3])
    is_night = sys.argv[4] == "night"
    night = "night" if is_night else "day"
    info("Handling patient %s with mesh %d, beta=%.1e at %s" % (num, n, beta, night))

    # Use ./num as input directory and set output ./key as output
    # directory
    path = os.getcwd()
    patient = os.path.join(path, num)
    key = "results_ocd_pat%s_n%d_beta%.1e_%s" % (num, n, beta, night)
    output_dir = os.path.join(patient, key)
    info("Setting %s as output directory" % output_dir)

    # Create OCD(R) map
    (c0, c1, phi) = map_patient(patient, n, beta, is_night, output_dir)

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    info("\nMax memory usage (MB) = %g\n" % (mem/1024))

if __name__ == "__main__":

    run_ocd()


    
