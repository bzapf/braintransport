import nibabel
from dolfin import *
from dolfin_adjoint import *
import pathlib
import numpy
from nibabel.affines import apply_affine
import abc
import os
from scipy.ndimage import gaussian_filter, median_filter
try:

    from scripts.simulations.run_nightly import find_files, get_delta_t
    from scripts.plots.definitions import firsts
    from scripts.simulations.run_nightly import Tmax

except ModuleNotFoundError:
    try:
        from run_nightly import find_files, get_delta_t
        from definitions import firsts
        from run_nightly import Tmax
    except ModuleNotFoundError:
        import sys
        sys.path.insert(0, "/cluster/home/bazapf/nobackup/sleepCode/")
        from scripts.simulations.run_nightly import find_files, get_delta_t
        from scripts.plots.definitions import firsts
        from scripts.simulations.run_nightly import Tmax


parameters['allow_extrapolation'] = True



def float_to_filestring(t):
    """format t to have two digits before and two digits after comma

    Args:
        t ([type]): [description]

    Returns:
        string: formated time
    """
    if type(t) == str:
        t = float(t)
    t = format(t, ".2f")
    if len(t) < 5:
        t = "0" + t
    return t



def interpolate_bc_guess(data, k, params):
    
    t_ = 0
    g_list = []
    idx = 0
    g0 = Function(data.function_space)

    # extended_data = data
    data.measurements["0.00"] = g0
    data.time_filename_mapping["0.00"] = None

    
    for i in range(k):
        # print(idx)
        t_ += params["dt"]

        if t_ > data.measurement_points()[idx + 1] and not (t_ > params["T"]):
            idx += 1

        # if t_ > extended_t[idx + 1] and not (t_ > params["T"]):
        #     idx += 1

        g = Function(data.function_space)
        
        # g_i = extended_data[extended_t[idx]].vector()[:]
        # g_ii = extended_data[extended_t[idx + 1]].vector()[:]
        # t_i = extended_t[idx]
        # t_ii = extended_t[idx + 1]

        g_i = data.get_measurement(data.measurement_points()[idx]).vector()[:]
        g_ii = data.get_measurement(data.measurement_points()[idx + 1]).vector()[:]
        t_i = data.measurement_points()[idx]
        t_ii = data.measurement_points()[idx + 1]
        inter_dt = t_ii - t_i
        # print("t_i", t_i, "t,", t_, "t_ii", t_ii)
        
        g.vector()[:] = g_i + (t_ - t_i) * (g_ii - g_i) / inter_dt
        # print(t_, t_i, t_ii, (t_ - t_i))

        g_list.append(g)


    data.measurements.pop("0.00")
    data.time_filename_mapping.pop("0.00")


    return g_list


def load_fenics_files(functionspace, hdfs, mesh2path, keys_to_load=None):
    """_summary_

    Args:
        functionspace (FunctionSpace): Space to project loaded functions onto
        hdfs (_type_): folder to hdf files storing the solution
        mesh2path (_type_): path to mesh 
        keys_to_load (_type_, optional): If not None, load only specific timepoints

    Returns:
        _type_: _description_
    """

    simulations = {}

    mesh2 = Mesh(mesh2path)
    V = FunctionSpace(mesh2, "CG", 2)

    files = sorted(os.listdir(hdfs), key=lambda x: float(x[:-7]))

    if keys_to_load is not None:
        files_ = []
        for file in files:
            t = float_to_filestring(float(file[:-7]) / 3600)
            if t in keys_to_load:
                files_.append(file)
        files = files_

    assert len(files) > 1

    print("Will load files:")
    for file in files:
        print(file)

    for file in files:
    
        t = float_to_filestring(float(file[:-7]) / 3600)

        # assert t in mapping.keys()
        
        # if t == float_to_filestring(float(file[:-7]) /  3600):

        fun = Function(V)
        hdf5 = HDF5File(mesh2.mpi_comm(), hdfs + file, "r")
        hdf5_name = "u"
        if not hdf5_name.startswith("/"):
            hdf5_name = "/" + hdf5_name
        
        hdf5.read(fun, hdf5_name)
        hdf5.close()

        # breakpoint()

        simulations[t] = project(fun, V=functionspace)
        
    return simulations








class Nan_Filter():
    def __init__(self, maskpath):

        # raise ValueError

        assert os.path.isfile(maskpath)

        ajdacent_idx = []
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:
                    ajdacent_idx.append([x, y, z])
        ajdacent_idx.remove([0, 0, 0])
        self.ajdacent_idx = numpy.array(ajdacent_idx)

        self.mask = numpy.load(maskpath)


    def __call__(self, data, ijk, i, j, k):
                
        data = numpy.where(self.mask, data, numpy.nan)

        nan_idx = numpy.argwhere(numpy.isnan(data[i, j, k]))

        nan_ijk = ijk[:, nan_idx[:, 0]]
        ni, nj, nk = numpy.rint(nan_ijk).astype("int")

        data2 = numpy.copy(data)
        idx = numpy.zeros_like(self.ajdacent_idx)

        for x, y, z in zip(ni, nj, nk):
            
            idx[:, 0] = self.ajdacent_idx[:, 0] + x
            idx[:, 1] = self.ajdacent_idx[:, 1] + y
            idx[:, 2] = self.ajdacent_idx[:, 2] + z
            
            sr = data[idx[:, 0], idx[:, 1], idx[:, 2]]
            
            sr = sr[~numpy.isnan(sr)]

            assert sr.size > 1

            data2[x, y, z] = numpy.median(sr)


        return data2







def read_image(filename, space, data_filter=None):
    
    print("Loading", filename)
    
    image2 = nibabel.load(filename)
    data = image2.get_fdata()

    
    # breakpoint()

    # if apply_filter:
    #     print("filtering data")
    #     data = gaussian_filter(data, sigma=sigma)


    # print(filename, format(data[mask].sum(), ".2e"))

    u_data = Function(space)
    ras2vox = image2.header.get_ras2vox()

    ras2vox_tkr_inv = numpy.linalg.inv(image2.header.get_vox2ras_tkr())
    # if tkr is True:
    #     ras2vox = ras2vox_tkr_inv
    ras2vox = ras2vox_tkr_inv

    xyz = space.tabulate_dof_coordinates()
    ijk = apply_affine(ras2vox, xyz).T
    i, j, k = numpy.rint(ijk).astype("int")
    
    if data_filter is not None:
        data = data_filter(data, ijk, i, j, k)
        u_data.vector()[:] = data[i, j, k]
    else:

        print("No filter used, setting", numpy.where(numpy.isnan(data[i, j, k]), 1, 0).sum(), "/", i.size, " nan voxels to 0")
        data[i, j, k] = numpy.where(numpy.isnan(data[i, j, k]), 0, data[i, j, k])
        u_data.vector()[:] = data[i, j, k]
        #data = my_filter(data, ijk, i, j, k)
    
    # breakpoint()
    #assert numpy.isnan(data[i, j, k]).sum() == 0

    

    
    # breakpoint()

    # nnan = numpy.sum(numpy.where(numpy.isnan(u_data.vector()[:]), 1, 0))
    # assert nnan == 0
    
    # print("nan in u_data", nnan)
    # print("in percent", nnan / u_data.vector()[:].size)

    # # IF NAN I SET VALUES EQUAL TO LEFT NEIGHBOR
    # if numpy.isnan(u_data.vector()[0]):
    #     nans = numpy.isnan(u_data.vector())
    #     real = ~nans
    #     real_idx = numpy.where(real)[0]
    #     u_data.vector()[0] = u_data.vector()[real_idx[0]]
    
    # for i in range(1, len(u_data.vector()[:])):
    #     if numpy.isnan(u_data.vector()[i]):
    #         u_data.vector()[i] = 10
    #         breakpoint()
    #         # breakpoint()
    #         print("processing nan")
    #         # u_data.vector()[i] = u_data.vector()[i-1]


    return u_data



class Measurements(abc.ABC):
    def __init__(self, function_space, expected_D=None):
        
        self.function_space = function_space
        # mesh = function_space.mesh()
        self.expected_D = expected_D

    def get_measurement(self, t):
        raise NotImplementedError
        # u = Function(self.function_space)
        # return u
    
    def measurement_points(self):
        raise NotImplementedError
        # return []

    def dump_pvd(self, vtkpath):
        """Dump all data snapshots to a pvd file

        Args:
            vtkpath (str): Path to export to, without file extension. Example: /home/tux/test_measurements
        """

        assert not os.path.isdir(vtkpath)
        assert not vtkpath.endswith("/")
        
        u_ = Function(self.function_space)

        vtkfile = File(vtkpath + '.pvd')
        
        for t in self.measurement_points():
            u = self.get_measurement(t)
            # breakpoint()
            u_.assign(u)
            vtkfile << u_ # float(t))




class MRI_Measurements(Measurements):

    def __init__(self, params, 
                data_filter, print_message=True,
                *args, **kwargs):
        
        super(MRI_Measurements, self).__init__(*args, **kwargs)
        
        self.params = params

        self.print_message = print_message

        files = sorted(find_files(params["pat"], params["concfolder"], data_folder=params["path_to_files"]))

        if params["concfolder"] == "FIGURES":
            assert False not in ["_masked_increase" in x for x in files]
            files = [firsts[params["pat"]] + "_masked_increase"] + [pathlib.Path(x).stem for x in files]

        else:
            assert False not in ["concentration" in x for x in files]

            files = [firsts[params["pat"]] + "_concentration"] + [pathlib.Path(x).stem for x in files]

        # NOTE "00.00" is the IC, should not be included here
        self.time_filename_mapping = {}

        for _, f in enumerate(files[1:]):
            # "dt = t2 - t1"
            # idx = idx + 1
            # dt = get_delta_t(files[idx - 1], files[idx])

            dt = get_delta_t(files[0], f)

            # breakpoint()

            if dt > Tmax:
                print("Omit image", f)
                continue

            # key = format(float(dt) / 3600, "0.2f")
            
            key = self.timeformat(dt)

            # while len(key) < 5:
            #    key = "0" + key

            self.time_filename_mapping[key] = f + ".mgz"

            if self.print_message:

                print("Added ", key, "=", self.time_filename_mapping[key], "to self.time_filename_mapping")

            assert dt > 0

        self.data_filter = data_filter

        self.read_fun = read_image

        self.measurements = {}


    def timeformat(self, t):
        return format(t, ".2f")

    def get_measurement(self, t):

        try:
            t = self.timeformat(t)
            return self.measurements[t]

        except KeyError:
        
            if isinstance(t, float):
                t = self.timeformat(t)

            assert t in self.time_filename_mapping.keys()
            
            filename = self.params["path_to_data"] + self.time_filename_mapping[t]

            

            self.measurements[t] = self.read_fun(filename, space=self.function_space, # tkr=True, 
                                                data_filter=self.data_filter,
                                                )
            
        return self.measurements[t]
    
    def measurement_points(self):
        
        # Use 
        # list(map(lambda x: float(x), list(self.time_filename_mapping.keys())))
        # to get a list of floats
        
        return sorted(list(map(lambda x: float(x), list(self.time_filename_mapping.keys()))))

