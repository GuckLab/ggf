"""Decomposition of stress in Legendre polynomials

To compute the GGF, :func:`ggf.globgeomfact.coeff2ggf` uses the
coefficients of the decomposition of the stress into Legendre
polynomials :math:`P_n(\text{cos}(\theta))`. This example visualizes
the small differences between the original stress and the stress
computed from the Legendre coefficients. This plot is automatically
produced by the original Matlab script *StretcherNStress.m*.
"""
import matplotlib.pylab as plt
import numpy as np
import percache

from ggf.sci_funcs import legendrePlm
from ggf.core import stress


@percache.Cache("stress_decomposition.cache", livesync=True)
def compute(**kwargs):
    "Locally cached version of ggf.core.stress"
    return stress(**kwargs)


# compute default stress
theta, sigmarr, coeff = compute(ret_legendre_decomp=True,
                                numpoints=300)

# compute stress from coefficients
numpoints = theta.size
sigmarr_c = np.zeros((numpoints, 1), dtype=float)
for ii in range(numpoints):
    for jj, cc in enumerate(coeff):
        sigmarr_c[ii] += coeff[jj] * \
            np.real_if_close(legendrePlm(0, jj, np.cos(theta[ii])))

# polar plot
plt.figure(figsize=(8, 8))
ax = plt.subplot(111, projection="polar")
plt.plot(theta, sigmarr, '-b', label="computed stress")
plt.plot(theta + np.pi, sigmarr[::-1], '-b')
plt.plot(theta, sigmarr_c, ':r', label="from Legendre coefficients")
plt.plot(theta + np.pi, sigmarr_c[::-1], ':r')
plt.legend()

plt.tight_layout()
plt.show()
