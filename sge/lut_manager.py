import multiprocessing as mp
import os
import socket

import h5py
import jobmanager as jm
import numpy as np

import ggf


USERNAME = os.environ["USER"]
#SERVER = "127.0.0.1"
SERVER = "guck-paulm-pc"
AUTHKEY = "d10fj31"
PORT = 42521
NCPUS = mp.cpu_count()
MPIRUN = "mpirun"

RUNMOD = RUNLIB = RUNEXP = ""

# try to put data there:
SERVERDIR = os.environ["HOME"]+"/sge_out"

try:
    MYIP = [(s.connect(('8.8.8.8', 80)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]
except OSError:
    MYIP = "127.0.0.1"

try:
    SERVERIP = socket.gethostbyname(SERVER)
except:
    SERVERIP = "127.0.0.1"


def map_lut2geom(kwargs):
    kw = kwargs.copy()
    relative_object_index = kw.pop("relative_object_index")
    kw["object_index"] = relative_object_index * kw["medium_index"]

    stretch_ratio = kw.pop("stretch_ratio")
    kw["semi_major"] = stretch_ratio * kw["semi_minor"] + kw["semi_minor"]
    return kw


class PM_Client(jm.JobManager_Client):
    def __init__(self):
        super(PM_Client, self).__init__(server=SERVER, 
                                        authkey=AUTHKEY, 
                                        port=PORT, 
                                        nproc=1,
                                        no_warnings=True, 
                                        verbose=0)

    @staticmethod
    def func(args, const_arg):
        try:
            result = ggf.get_ggf(*args[2:])
        except:
            result = np.nan
        return (MYIP, USERNAME, result)


class PM_Server(jm.JobManager_Server):
    def __init__(self, server_args):
        
        # setup init parameters for the ancestor class
        fname_dump = None
        verbose = 1
        msg_interval = 1
        
        # init ancestor class
        super(PM_Server, self).__init__(authkey=AUTHKEY,
                                        port=PORT,
                                        const_arg={},
                                        verbose=verbose,
                                        msg_interval=msg_interval,
                                        fname_dump=fname_dump)

        for lut_name in server_args:
            lut_set = server_args[lut_name]
            # default keyword arguments
            default_kw = {}
            for key in lut_set:
                if not isinstance(lut_set[key], list):
                    default_kw[key] = lut_set[key]
            # initiate h5 file and get meshgrid with label axes
            out, labels, mesh = self.get_output_h5(lut_name, lut_set)
            for ii in range(mesh[0].size):
                if np.isnan(out.flat[ii]):
                    # build kwargs
                    kw = default_kw.copy()
                    for jj, ll in enumerate(labels):
                        kw[ll] = mesh[jj].flat[ii]
                    putargs = []
                    # convert kwargs to ggf kwargs (stretch_ratio, relative_index)
                    kw_ggf = map_lut2geom(kw)
                    # convert kwargs to args
                    for kk in ggf.get_ggf.__code__.co_varnames:
                        if kk in kw_ggf:
                            putargs.append(kw_ggf[kk])
                    self.put_arg(tuple([lut_name, ii] + putargs))

    def get_h5_filename(self, lut_name):
        return "lut_{}.h5".format(lut_name)

    def get_output_h5(self, lut_name, lut_set=None):
        fn = self.get_h5_filename(lut_name)
        if os.path.exists(fn) and lut_set is None:
            with h5py.File(fn, mode="r") as h5:
                data = h5["lut"].value
        else:
            assert lut_set is not None, "First call must include 'lut_set'"
            labels = [ll for ll in lut_set if isinstance(lut_set[ll], list)]
            labels = sorted(labels)
            # setup meshgrid for client args
            meshgrid_args = []
            dims = []
            for ll in labels:
                vals = np.linspace(*lut_set[ll])
                meshgrid_args.append(vals)
                dims.append(len(vals))
            mesh = np.meshgrid(*meshgrid_args)
            # initialize lut
            with h5py.File(fn, mode="a") as h5:
                if "lut" not in h5:
                    lut = np.nan * np.zeros(dims, dtype=float)
                    h5lut = h5.create_dataset("lut", data=lut)
                    h5lut.attrs["dimension_order"] = ",".join(labels)
                    for ll in labels:
                        h5lut.attrs["{} min".format(ll)] = lut_set[ll][0]
                        h5lut.attrs["{} max".format(ll)] = lut_set[ll][1]
                        h5lut.attrs["{} num".format(ll)] = lut_set[ll][2]
                else:
                    lut = h5["lut"].value
            data = lut, labels, mesh
        return data
    
    def save_to_output_h5(self, lut_name, coord, value):
        fn = self.get_h5_filename(lut_name)
        with h5py.File(fn, mode="a") as h5:
            lut = h5["lut"].value
            lut.flat[coord] = value
            h5["lut"][:] = lut

    def process_new_result(self, arg, result):
        _client_ip, _client_user, res = result
        idset = arg[0]
        coord = arg[1]
        self.save_to_output_h5(idset, coord, res)
