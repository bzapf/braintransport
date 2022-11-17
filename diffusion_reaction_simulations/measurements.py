import nibabel
from dolfin import *
from dolfin_adjoint import *
import pathlib
import numpy
from nibabel.affines import apply_affine
import abc
import os
import datetime

firsts = {
    "068": '20160912_080147',
    "002": '20151007_075055',
    "078": '20161031_074500',
    "091": '20170116_074933',
    "127": '20170925_075956',
    "172": '20180212_081154',
    "249": '20191111_075145',
    "105": '20170508_081611',
    "175": '20180416_080959',
    "176": '20180514_080612',
    "178": '20180430_075417',
    "183": '20181029_074933',
    "190": '20180528_074812',
    "191": '20180806_080756',
    "205": '20181008_074949',
    "215": '20190128_073807',
    "218": '20190114_074903',
    "228": '20190304_074149',
    "240": '20191028_074716',
    "199": '20190916_075735',
    "227": '20190401_074904',
    "230": '20190527_075031',
    "235": '20190603_075620',
    "236": '20190624_081053',
    "241": '20191216_075927'
 }


Tmax = 2.5 * 24 * 3600
parameters['allow_extrapolation'] = True

def find_files(pat, concfolder, data_folder):
    folder = data_folder + '/%s/' % pat + concfolder + '/'
    files = os.listdir(folder)
    new_files = [folder + files[i]
                 for i in range(len(files)) if files[i].startswith('2')]
    new_files = sorted(new_files)
    return new_files

def get_delta_t(f1, f2):
    frmt = '%Y%m%d_%H%M%S'
    A1 = f1.split("/")[-1]
    A2 = f2.split("/")[-1]

    try:
        date_time1 = A1.split('_concentration')[0]
        date_time2 = A2.split('_concentration')[0]
        date_obj1 = datetime.strptime(date_time1, frmt)
        date_obj2 = datetime.strptime(date_time2, frmt)      
    except ValueError:
        date_time1 = A1.split('_masked_increase')[0]
        date_time2 = A2.split('_masked_increase')[0]
        date_obj1 = datetime.strptime(date_time1, frmt)
        date_obj2 = datetime.strptime(date_time2, frmt)


    
    difference = date_obj2 - date_obj1
    time = difference.days * 3600 * 24 + difference.seconds
    return time


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

        g = Function(data.function_space)
        
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

        fun = Function(V)
        hdf5 = HDF5File(mesh2.mpi_comm(), hdfs + file, "r")
        hdf5_name = "u"
        if not hdf5_name.startswith("/"):
            hdf5_name = "/" + hdf5_name
        
        hdf5.read(fun, hdf5_name)
        hdf5.close()

        simulations[t] = project(fun, V=functionspace)
        
    return simulations





def read_image(filename, space, data_filter=None):
    
    print("Loading", filename)
    
    image2 = nibabel.load(filename)
    data = image2.get_fdata()


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
            u_.assign(u)
            vtkfile << u_




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

        self.time_filename_mapping = {}

        for _, f in enumerate(files[1:]):

            dt = get_delta_t(files[0], f)

            if dt > Tmax:
                print("Omit image", f)
                continue
            
            key = self.timeformat(dt)

            self.time_filename_mapping[key] = f + ".mgz"

            if self.print_message:

                print("Added ", key, "=", self.time_filename_mapping[key], "to self.time_filename_mapping")

            # NOTE "00.00" is the IC, should not be included here
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

            self.measurements[t] = self.read_fun(filename, space=self.function_space, data_filter=self.data_filter)
            
        return self.measurements[t]
    
    def measurement_points(self):
        
        return sorted(list(map(lambda x: float(x), list(self.time_filename_mapping.keys()))))

