"""Other scientific functions"""
import numpy as np

from .matlab_funcs import legendre


def legendrePlm(m, l, x): 
    a = np.concatenate((np.atleast_1d(np.zeros(m)),
                        np.atleast_1d(np.array(1)),
                        np.atleast_1d(np.zeros(l-m))))
    b = legendre(l, x)
    return np.real_if_close(np.dot(a,b)[0])

