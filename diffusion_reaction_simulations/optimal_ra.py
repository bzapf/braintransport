
from dolfin import *
from dolfin_adjoint import *
import numpy as np
import argparse
import os
import numpy as np
import json
import pandas
import time
from copy import deepcopy
import json
from measurements import MRI_Measurements, interpolate_bc_guess
from parameters import PositiveParam

set_log_level(LogLevel.ERROR)


def store_stuff(context):
    global pat
    global params
    global data

    region_volumes = {}
    region_areas = {"avgds": Surface}
    region_volumes["avg"] = Volume
    region_volumes["white"] = white_volume
    region_volumes["gray"] = gray_volume
    region_volumes["brainstem"] = brain_stem_volume

    with open(context.outfolder + 'region_volumes.json', 'w') as f:
        json.dump(region_volumes, f)

    with open(context.outfolder + 'region_areas.json', 'w') as f:
        json.dump(region_areas, f)
                
    expfile = open('%s/experimental_data.txt' %
                context.outfolder, 'w')  # Data read from experiments
    expfile.write('t avg avgds gray white brainstem ')

    expfile.write('\n')

    for tkey, image in data.measurements.items():

        exptime = float(tkey)

        print("Writing exptime", exptime / 3600)

        expfile.write('%g ' % exptime)
        average = assemble(image*dx)/Volume
        average_ds = assemble(image*ds)/Surface
        expfile.write('%g ' % average)
        expfile.write('%g ' % average_ds)
        gray_avg = assemble(image*context.dx_SD(1))/gray_volume
        white_avg = assemble(image*context.dx_SD(2))/white_volume
        brain_stem_avg = assemble(image*context.dx_SD(3))/brain_stem_volume
        
        expfile.write('%g ' % gray_avg)
        expfile.write('%g ' % white_avg)
        expfile.write('%g ' % brain_stem_avg)

        expfile.write('\n')

    expfile.close()






class Context(object):
    def __init__(self, # mesh_config, 
                    params, V, delta_alpha, delta_r_d, g_list, data, outfolder, MD, Kt, dx_SD):
        """
        :param mesh_config: dictionary which contains ds, dx
        :param V: FunctionSpace for state variable
        :param delta: diffusion coefficient
        :param params: dictionary which contains time step size dt, regularization parameter alpha
        :param g_list: boundary conditions
        :param data: class
        """

        self.added_times = set()

        self.MD = MD
        self.dx_SD = dx_SD
        self.Kt = Kt

        self.params = params
        self.ds = ds # mesh_config["ds"]
        self.dx = dx # mesh_config["dx"]
        self.V = V
        self.write_solution = False
        self.delta_alpha, self.delta_r_d = delta_alpha, delta_r_d

        self.g_list = g_list
        self.current_g_index = 0
        self.g = self.g_list[self.current_g_index]
        self.t = 0
        self.dt = params["dt"]
        self.T = params["T"]

        self.data = data
        self.tau = data.measurement_points()
        self.next_tau = 0
        self.d = Function(V)
        self.J = 0.0
        self.datanorm = 0.0

        # self.linear_solver_args = ("gmres", "amg")
        self.state_at_measurement_points = []
        self.outfolder = outfolder

        

    def check_termination(self):
        return (not self.t + self.dt / 2 <= self.T)

    def save_state(self, u):
        self.state_at_measurement_points.append(u)


    def advance_time(self):
        """
        update of time variables and boundary conditions
        """
        self.t += self.dt     # update time-step
        self.current_g_index += 1  # update index for g
        self.g = self.g_list[self.current_g_index]  # update Dirichlet-BC


    def update_objective(self, Uu):
        """
        update objective
        """

        if not self.next_tau == len(self.tau) and abs(self.t - self.tau[self.next_tau]) < abs(
                self.t + self.dt - self.tau[self.next_tau]):
            
            # global mesh
            # comm = mesh.mpi_comm()
            # rank = comm.Get_rank()
            # if rank == 0:

            print("Adding to objective,","self.next_tau", self.next_tau)

            self.d.assign(self.data.get_measurement(self.tau[self.next_tau]))  # Read observation
            
            J_d = assemble((Uu + self.g_list[self.current_g_index] - self.d) ** 2 * self.dx) 

            datanorm = assemble((self.d) ** 2 * self.dx)
            
            self.datanorm += datanorm
            
            # comm = mesh.mpi_comm()
            # rank = comm.Get_rank()
            # if rank == 0:
            print("self.J", self.J, "J_d", J_d)
            self.J += J_d

            self.added_times.add((self.t, self.tau[self.next_tau]))
            
            # save state at measurement points for testing reasons
            self.save_state(project(Uu + self.g_list[self.current_g_index], Uu.function_space()))
            
            # Move on to next observation
            self.next_tau += 1

    def boundary_condition(self):
        return DirichletBC(self.V, Constant(0.0), 'on_boundary')

    def return_value(self):
        return self.J / self.datanorm


def forward(context, write_solution=False, filename=None):
    V = context.V
    # Define trial and test-functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Solution at current and previous time
    U_prev = Function(V)
    U = Function(V)

    dt = context.dt
    dx = context.dx

    fac = 1 / 2 # for Crank-Nicolson method
    iter_k = 0

    if write_solution and store:

        if filename == "plainstate":
            concfile = open('%s/concs_plain.txt' % context.outfolder, 'w')
        else:
            concfile = open('%s/concs.txt' % context.outfolder, 'w')
        concfile.write('t avg avgds gray white brainstem ')
        concfile.write('\n')

    # term 1
    a = inner(u, v) * dx

    # term 2

    if params["global"]:
        scale = alpha(context.delta_alpha)
    else:
        assert params["gray"]
        scale = 1

    a += fac * dt * alpha(context.delta_alpha) * D_scale * context.MD * inner(grad(u), grad(v)) * context.dx_SD(1)
    
    if context.params["dti"]:
        a += fac * dt * D_scale * inner(dot(context.Kt, grad(u)), grad(v)) * context.dx_SD(2)
        a += fac * dt * D_scale * context.MD * inner(grad(u), grad(v)) * context.dx_SD(3)
    else:
        assert context.params["nodti"]
        a += fac * dt * scale * D_scale * MD_white * inner(grad(u), grad(v)) * context.dx_SD(2)
        a += fac * dt * scale * D_scale * MD_brainstem * inner(grad(u), grad(v)) * context.dx_SD(3)


    a += fac * dt * reaction_rate(context.delta_r_d) * inner(u, v) * dx
    bc = context.boundary_condition()
    A = assemble(a)
    bc.apply(A)

    solver = LUSolver()
    solver.set_operator(A)
    # solver = PETScKrylovSolver('gmres', 'amg')
    # solver.set_operators(A, A)

    while not context.check_termination():

        print("iter k", iter_k)

        iter_k += 1

        sol = project(U + context.g_list[context.current_g_index], U.function_space(), solver_type="lu")

        sol.rename("sol", "sol")

        U_prev.assign(U)

        context.advance_time()
        
        # current g index
        g_ind = context.current_g_index
        glist = context.g_list

        # term 6
        L = (U_prev + glist[g_ind - 1]) * v * dx
        
        # term 1
        L -= inner(context.g_list[context.current_g_index], v) * dx

        # term 2
        L -= fac * dt * alpha(context.delta_alpha) * D_scale * context.MD * inner(grad(context.g_list[context.current_g_index]), grad(v)) * context.dx_SD(1)
        
        if context.params["dti"]:
            L -= fac * dt * D_scale * inner(dot(context.Kt, grad(context.g_list[context.current_g_index])), grad(v)) * context.dx_SD(2)
            L -= fac * dt * D_scale * context.MD * inner(grad(context.g_list[context.current_g_index]), grad(v)) * context.dx_SD(3)        
        else:
            assert context.params["nodti"]
            L -= fac * dt * scale * D_scale * MD_white * inner(grad(context.g_list[context.current_g_index]), grad(v)) * context.dx_SD(2)
            L -= fac * dt * scale * D_scale * MD_brainstem * inner(grad(context.g_list[context.current_g_index]), grad(v)) * context.dx_SD(3)

        if write_solution and store:

            print(iter_k, context.t)

            concfile.write('%g ' % context.t)
            average = assemble(sol*dx)/Volume
            average_ds = assemble(sol*ds)/Surface
            gray_avg = assemble(sol*dx_SD(1))/gray_volume
            white_avg = assemble(sol*dx_SD(2))/white_volume
            brain_stem_avg = assemble(sol*dx_SD(3))/brain_stem_volume

            concfile.write('%g ' % average)
            concfile.write('%g ' % average_ds)
            concfile.write('%g ' % gray_avg)
            concfile.write('%g ' % white_avg)
            concfile.write('%g ' % brain_stem_avg)
            concfile.write('\n')


        # term 4
        L -= fac * dt * alpha(context.delta_alpha) * D_scale * context.MD * inner(grad(U_prev + glist[g_ind - 1]), grad(v)) * context.dx_SD(1)
        if context.params["nodti"]:
            L -= fac * dt * D_scale * scale * MD_white * inner(grad(U_prev + glist[g_ind - 1]), grad(v)) * context.dx_SD(2)
            L -= fac * dt * D_scale * scale * MD_brainstem * inner(grad(U_prev + glist[g_ind - 1]), grad(v)) * context.dx_SD(3)
        else:
            assert context.params["dti"]
            L -= fac * dt * D_scale * inner(dot(context.Kt, grad(U_prev + glist[g_ind - 1])), grad(v)) * context.dx_SD(2)
            L -= fac * dt * D_scale * context.MD * inner(grad(U_prev + glist[g_ind - 1]), grad(v)) * context.dx_SD(3)

        L -= fac * dt * reaction_rate(context.delta_r_d) * inner(context.g_list[context.current_g_index], v) * dx
        L -= fac * dt * reaction_rate(context.delta_r_d) * inner(U_prev + glist[g_ind - 1], v) * dx
    
        # Assemble RHS and apply DirichletBC
        b = assemble(L)
        bc.apply(b)

        solver.solve(U.vector(), b)

        # solve(A, U.vector(), b, 'lu')
        
        context.update_objective(U)

    context.t = 0
    context.current_g_index = 0
    context.next_tau = 0


    if write_solution and store:

        concfile.close()

        store_stuff(context=context)

        assert filename is not None
        states = context.state_at_measurement_points
        assert len(states) > 1

        assert context.outfolder is not None
        file = File(context.outfolder + filename + '.pvd')

        Hdf = HDF5File(mesh.mpi_comm(), context.outfolder + filename + ".hdf", "w")
        Hdf.write(mesh, "mesh")
        
            
        for idx, state in enumerate(states):
            state.rename(filename, filename)
            file << state

            Hdf.write(state, "state" + str(idx))

            obs = context.data.get_measurement(context.data.measurement_points()[idx])
            
            err = assemble(((state-obs)**2)*dx)
            # print(err)
            try:
                params[filename + "err"].append(err)
            except:
                params[filename + "err"] = [err]

        Hdf.close()


    return context.return_value()


def make_exportfolder(exportpath):

    subfolder = exportpath
    print("Creating subfolder", subfolder)
    try:
        os.makedirs(subfolder, exist_ok=True)
    except FileExistsError:
        pass

    
    return subfolder


def iter_cb(m):
    global alpha_list
    global iter_cnt
    global store
    global context
    global params



    if iter_cnt == 4:
        dt = time.time() - t0
        time_per_iter = dt / (iter_cnt + 1)
        
        print("First", iter_cnt + 1, "iterations took", format(dt / 60, ".0f"), "minutes")

        expected_time_in_hours = time_per_iter * params["iters"] / 3600
        params["expected_time_in_hours"] = expected_time_in_hours
        
        print("Optimization for", params["iters"], " will probably take", format(expected_time_in_hours, ".2f"), "hours")
    
    vals = m[0]
    alpha_list.append(alpha(vals))
    
    rd = m[1]
    r_d_list.append(reaction_rate(rd))
    

    print("Iter", iter_cnt, " Alpha =", format(alpha(vals), ".2e"),
            " R_d =", format(reaction_rate(rd), ".2e"),
            )
        

    if store:
        global outfolder
        with open(outfolder + "alpha_during_optim.txt", "a") as myfile:
            myfile.write(str(alpha_list[-1]) + ",")

        with open(outfolder + "r_d_during_optim.txt", "a") as myfile:
            myfile.write(str(r_d_list[-1]) + ",")


    iter_cnt += 1               


if __name__ == "__main__":


    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--store", action="store_true", default=False)

    parser.add_argument("--pat", default="205", type=str)

    parser.add_argument("--exportpath", type=str, help="global path to store results to", default=None)

    parser.add_argument("--mesh", type=str, help="meshname", default="parenchyma4_with_DTI.h5")
    parser.add_argument("--concfolder", type=str, choices=["FIGURES_CONC", "FIGURES_CONC_LUT", "FIGURES_CONC_WGB", "FIGURES_CONC_T1"],
                default="FIGURES_CONC_LUT", help="concentration folder")
                        
    parser.add_argument("--path_to_files", type=str, help="patfolders", default="/home/basti/Dropbox (UiO)/Sleep/")

    parser.add_argument("--iter_k", type=int, default=48, help="time steps")
    parser.add_argument("--outfolder", type=str, default="outs", help="name of folder under exportpath/ to store pvd files to")
    parser.add_argument("--overwrite", action="store_true", default=False)
    parser.add_argument("--iters", default=20, type=float, help="max iterations of lbfgs")

    parser.add_argument("--global", action="store_true", default=False)
    parser.add_argument("--gray", action="store_true", default=False)

    parser.add_argument("--dti", action="store_true", default=False)
    parser.add_argument("--nodti", action="store_true", default=True)

    parser.add_argument("--alpha", help="perform forward simulatino with given alpha, needs --reaction_rate", default=0., type=float)
    parser.add_argument("--reaction_rate", help="perform forward simulatino with given reaction rate, needs --alpha", default=0., type=float)

    parser.add_argument("--alpha_init", default=3., type=float)
    parser.add_argument("--r_d_init", default=1e-5, type=float)
        
    params = vars(parser.parse_args())

    pat = params["pat"]

    if params["gray"]:
        assert not params["global"]
    
    if params["global"]:
        assert not params["gray"]

    exportpath = params["exportpath"]
    if exportpath is not None:
        params["store"] = True

    store = params["store"]

    if exportpath is not None and not exportpath.endswith("/"):
        exportpath += "/"

    if params["overwrite"] and params["store"] and os.path.isdir(exportpath):
        os.system("rm -r " + exportpath)

    if params["store"]:
        exportfolder = make_exportfolder(exportpath)
        params["outfolder"] = os.path.join(exportfolder, params["outfolder"])

    params["meshpath"] = os.path.join(params["path_to_files"], params["pat"], "mesh", params["mesh"])
    print(params["meshpath"])
    assert os.path.isfile(params["meshpath"])

    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), params["meshpath"], "r")
    hdf.read(mesh, "/mesh", False)
    SD = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(SD, "/subdomains")

    try:
        lookup_table = MeshFunction("size_t", mesh, mesh.topology().dim())
        hdf.read(lookup_table, '/lookup_table')
        print("Loaded lookup table")
        lookup_available = True
    except:
        print("No lookup table found in meshfile", params["meshpath"])
        lookup_available = False

    if lookup_available:
        write_regions = np.unique(lookup_table.array()[:])
        dx_lookup = Measure('dx')(domain=mesh, subdomain_data=lookup_table)

    else:
        write_regions = []

    # GRAY = 1. WHITE = 2. BRAIN STEM = 3.
    dx_SD = Measure('dx')(domain=mesh, subdomain_data=SD)

    V = FunctionSpace(mesh, "CG", 1)
    TensorSpace = TensorFunctionSpace(mesh, 'DG', 0)
    MDSpace = FunctionSpace(mesh, 'DG', 0)
    MD = Function(MDSpace)

    Volume = assemble(Constant(1)*dx_SD(domain=mesh))
    Surface = assemble(Constant(1)*ds(domain=mesh))
    gray_volume = assemble(Constant(1)*dx_SD(1))
    white_volume = assemble(Constant(1)*dx_SD(2))
    brain_stem_volume = assemble(Constant(1)*dx_SD(3))

    if params["dti"]:
        params["nodti"] = False
        Kt = Function(TensorSpace)
        params["nodti"] = False
        hdf.read(MD, '/MD')
        hdf.read(Kt, '/DTI')
        print("Using DTI values")
    else:
        assert params["nodti"]
        MD.vector()[:] += 0.00092582074369026 # mean from all pats with DTI in gray
        MD_white = Function(MDSpace)
        MD_white.vector()[:] += 0.000867643624860488 # mean from all pats with DTI in white
        MD_brainstem = Function(MDSpace)
        MD_brainstem.vector()[:] += 0.0010602528519216016 # mean from all pats with DTI in brainstem
        print("Using patient cohort mean values for scalar D in white, gray, brainstem")        
        Kt = None

    assert params["dti"] is not params["nodti"]

    params["path_to_data"] = os.path.join(params["path_to_files"], params["pat"], params["concfolder"], "")

    assert os.path.isdir(params["path_to_data"])
    
    data = MRI_Measurements(params, function_space=V, data_filter=None)

    measurement_points = data.measurement_points()

    Tmax = max(measurement_points)
    T = Tmax
    k = params["iter_k"]

    print("Loaded data")

    if not T > 3600:
        raise ValueError("Time unit is seconds, so T should usually be > 1 h")
    
    outfolder = params["outfolder"]
    
    if not store:
        outfolder = None

    if outfolder is not None and not outfolder.endswith("/"):
        outfolder += "/"
    
    if outfolder is not None:
        try:
            os.makedirs(outfolder, exist_ok=True)
        except FileExistsError:
            pass

    if outfolder is not None:
        file = File(outfolder + 'measurements.pvd')
        for i in range(len(measurement_points)):
            func = data.get_measurement(measurement_points[i])
            func.rename("measure", "measure")
            file << func

    assert min(measurement_points) > 0, "Initial condition should not be part of measurements"
    
    if T == 0:
        raise Exception('check if there are measurements in the interval (0, Tmax]')

    params["T"] = T
    params["dt"] = T / k

    iter_cnt = 0
    
    alpha = PositiveParam(min=1)
    alpha_init = params["alpha_init"]
    delta_alpha_0 = Constant(alpha.f_inverse(alpha_init))
    delta_alpha = deepcopy(delta_alpha_0)
    alpha_list = [alpha_init]

    reaction_rate = PositiveParam(min=1e-7)
    r_d_init = params["r_d_init"]
    delta_r_d_0 = Constant(reaction_rate.f_inverse(r_d_init))
    delta_r_d = deepcopy(delta_r_d_0)
    r_d_list = [r_d_init]

    MD_water = 3e-3
    MD_gadovist = 3.8e-4 # (taken from Gd-DPTA)
    scale_diffusion_gad = MD_gadovist / MD_water
    D_scale = Constant(scale_diffusion_gad) 

    # ("Using data to interpolate starting guess forBC")
    g_listb = interpolate_bc_guess(data, k, params=params)

    g_list = [interpolate(Constant(0.0), V)] + g_listb

    # evaluate J
    context = Context(params, V, delta_alpha, delta_r_d, g_list, data, outfolder, MD=MD, Kt=Kt, dx_SD=dx_SD)

    if outfolder is not None and params["alpha"] > 0:
        """
        Only perform forward simulation for given --alpha and --reaction_rate
        """

        j_d_context = Context(params, V, delta_alpha=Constant(params["alpha"]), delta_r_d=Constant(params["reaction_rate"]), 
                        g_list=g_list, data=data, outfolder=outfolder, MD=MD, Kt=Kt, dx_SD=dx_SD)

        print("Evaluating forward pass")
        J_d = forward(j_d_context, write_solution=True, filename="finalstate")
        params["j_d_final"] = float(J_d)

        csv_e = pandas.read_csv('%s/experimental_data.txt' % context.outfolder, sep=' ', header=0)

        csv_c = pandas.read_csv(context.outfolder + '/concs' + '.txt', sep=' ', header=0)
        csv_c.to_csv(context.outfolder + '/concs' + '.csv')

        
        csv_e.to_csv('%s/experimental_data.csv' % context.outfolder)

        if store:
            with open(outfolder + 'params', 'w') as outfile:
                json.dump(params, outfile, sort_keys=True, indent=4)
                
        exit()

    else:


        print("Evaluating forward pass")
        J = forward(context)
        
        print('J', J)

        ctrls = [Control(delta_alpha), Control(delta_r_d)]

        print("Creating ReducedFunctional")
        Jhat = ReducedFunctional(J, ctrls)
        
        print("Evaluating reduced functional")
        t0_rf = time.time()
        
        jhat0 = Jhat([delta_alpha, delta_r_d])
        
        opt_time_rf = time.time() - t0_rf

        print("Evaluating reduced functional took", format(opt_time_rf / 60, ".0f"), "min", "(", format(opt_time_rf / 3600, ".0f"), "h ).")
        params["rf_eval_time_hours"] = opt_time_rf / 3600
        
        
        print('Jhat', jhat0,)

        params["j_d_init"] = jhat0


        j_d_plain_context = Context(params, V, delta_alpha=Constant(1), delta_r_d=Constant(0),
                        g_list=g_list, data=data, outfolder=outfolder, MD=MD, Kt=Kt, dx_SD=dx_SD)

        J_d_plain = forward(j_d_plain_context, write_solution=True, filename="plainstate")
        params["j_d_plain"] = float(J_d_plain)

        


        if store:
            with open(outfolder + "alpha_during_optim.txt", "a") as myfile:
                myfile.write(str(alpha_list[-1]) + ",")

            with open(outfolder + "r_d_during_optim.txt", "a") as myfile:
                myfile.write(str(r_d_list[-1]) + ",")


        t0 = time.time()

        if params["iters"] > 0:
            
            bounds = [[1, 1e-7], [10, 1e-3]]

            params["bounds"] = bounds


            opt_ctrls = minimize(Jhat, method="L-BFGS-B", callback=iter_cb, bounds=bounds,
                                    options={"disp": True, "maxiter": int(params["iters"]), "ftol": 1e-12, "gtol": 1e-6})

            opt_time = time.time() - t0
            print("Optimization took", format(opt_time / 60, ".0f"), "min", "(", format(opt_time / 3600, ".0f"), "h ).")
            params["optimization_time_hours"] = opt_time / 3600
            params["total_time_hours"] = (opt_time_rf + opt_time) / 3600

            vals = np.zeros(1)
            ctrls[0].data().eval(vals, np.zeros(1))
            params["alpha_final"] = alpha(vals[0])

            vals = np.zeros(1)
            ctrls[1].data().eval(vals, np.zeros(1))
            params["r_d_final"] = reaction_rate(vals[0])

            params["added_times"] = list(context.added_times)

            if outfolder is not None:

                print(delta_alpha, delta_r_d)
                print(params["alpha_final"], params["r_d_final"])
            
                j_d_context = Context(params, V, delta_alpha=opt_ctrls[0], delta_r_d=opt_ctrls[1], 
                                g_list=g_list, data=data, outfolder=outfolder, MD=MD, Kt=Kt, dx_SD=dx_SD)

                print("Evaluating forward pass")
                J_d = forward(j_d_context, write_solution=True, filename="finalstate")
                params["j_d_final"] = float(J_d)

            if store:
                csv_e = pandas.read_csv('%s/experimental_data.txt' % context.outfolder, sep=' ', header=0)
                csv_c = pandas.read_csv(context.outfolder + '/concs.txt', sep=' ', header=0)
                csv_c.to_csv(context.outfolder + '/concs.csv')

                csv_c = pandas.read_csv(context.outfolder + '/concs_plain.txt', sep=' ', header=0)
                csv_c.to_csv(context.outfolder + '/concs_plain.csv')
                
                csv_e.to_csv('%s/experimental_data.csv' % context.outfolder)


    if store:
        with open(outfolder + 'params', 'w') as outfile:
            json.dump(params, outfile, sort_keys=True, indent=4)
