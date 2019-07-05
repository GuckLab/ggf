"""lossily compress a LUT

Convert to float32 and use gzip compression.

Usage: python shrink_lut path/to/lut.h5 
"""
import sys

import h5py
import numpy as np


path_in = sys.argv[-1]
path_out = path_in[:-3] + "_compressed.h5"

with h5py.File(path_in, mode="r") as h5:
    lut = h5["lut"][:]
    attrs = dict(h5["lut"].attrs)

lut = np.array(lut, dtype=np.float32)

with h5py.File(path_out, "w") as h5:
    h5.create_dataset("lut", data=lut, compression="gzip", fletcher32=True)
    for key in attrs:
        h5["lut"].attrs[key] = attrs[key]