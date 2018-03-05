"""Decomposition of stress in Legendre polynomials

Output of coefficients for expansion in single Legendre Polynomials Pn[np.cos(th)]
"""
import matplotlib.pylab as plt
import numpy as np

from ggf.sci_funcs import legendrePlm
from ggf.core import stress

# compute default stress
theta, sigmarr, coeff = stress(ret_legendre_decomp=True)

# compute stress from coefficients
numpoints = theta.size
sigmarr_c = np.zeros((numpoints,1), dtype=float)
for ii in range(numpoints):
    for jj, cc in enumerate(coeff):
        sigmarr_c[ii] += coeff[jj]*np.real_if_close(legendrePlm(0,jj,np.cos(theta[ii])))

ax = plt.subplot(111, projection="polar")
plt.plot(theta, sigmarr, '-', label="computed stress")
plt.plot(theta, sigmarr_c,':', label="reconstruction from Legendre coefficients")
plt.legend()
plt.tight_layout()
plt.show()
