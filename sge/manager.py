import multiprocessing as mp
import os
import socket

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


class PM_Client(jm.JobManager_Client):
    def __init__(self):
        super(PM_Client, self).__init__(server=SERVER, 
                                        authkey=AUTHKEY, 
                                        port=PORT, 
                                        nproc=1,
                                        nice=19, 
                                        no_warnings=True, 
                                        verbose=1)

    @staticmethod
    def func(args, const_arg):
        try:
            result = ggf.get_ggf(*args[2:])
        except:
            result = np.nan
        return (MYIP, USERNAME, result)


class PM_Server(jm.JobManager_Server):
    def __init__(self, server_args, defaults):
        
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

        # output file identifier
        for arg in server_args:
            ids = []
            dims = []
            for aset in arg:
                ids.append("{}-{:.3e}-{:.3e}-{:05d}".format(aset[0][:3],
                                                            aset[1],
                                                            aset[2],
                                                            aset[3]))
                dims.append(aset[3])
            idset = "_".join(ids)
            # create an n-dimensional search map
            out = self.get_output_npy(idset, dims=dims)
            if len(dims) == 1:
                aset = arg[0]
                key = aset[0]
                vals = np.linspace(*aset[1:])
                for ii, val in enumerate(vals):
                    kw = defaults.copy()
                    kw[key] = val
                    if np.isnan(out[ii]):
                        putargs = []
                        for kk in ggf.get_ggf.__code__.co_varnames:
                            if kk in kw:
                                putargs.append(kw[kk])
                        print(putargs)
                        self.put_arg(tuple([idset, ii] + putargs))
            else:
                raise NotImplementedError("TODO")

    def get_npy_filename(self, idset):
        return "_results_{}.npy".format(idset)

    def get_output_npy(self, idset, dims=None):
        fn = self.get_npy_filename(idset)
        if os.path.exists(fn):
            data = np.load(fn)
        else:
            assert dims is not None, "First call must include 'dims'"
            data = np.nan * np.zeros(dims, dtype=float)
            np.save(fn, data)
        return data
    
    def save_to_output_npy(self, idset, coord, value):
        fn = self.get_npy_filename(idset)
        d = self.get_output_npy(idset)
        d.flat[coord] = value
        np.save(fn, d)

    def process_new_result(self, arg, result):
        _client_ip, _client_user, res = result
        idset = arg[0]
        coord = arg[1]
        self.save_to_output_npy(idset, coord, res)