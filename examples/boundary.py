"""Object boundary: stretching and Poisson's ratio

This example illustrates how the parameters Poisson's ratio
:math:`\nu` and stretch ratio :math:`\epsilon` influence
the object boundary used in :func:`ggf.core.stress` and
defined in :func:`ggf.core.boundary`.

Note that the boundary function was not defined correctly prior to
version 0.3.0 (`issue #1 <https://github.com/GuckLab/ggf/issues/1>`__).
Since version 0.3.0, the semi-minor axis is equivalent to the
keyword argument `a` (1 by default).
"""
import numpy as np
import matplotlib.pylab as plt

from ggf.stress.boyde2009.core import boundary

theta = np.linspace(0, 2*np.pi, 300)
costheta = np.cos(theta)

# change epsilon
eps = [.0, .05, .10, .15, .20]
b1s = []
for ep in eps:
    b1s.append(boundary(costheta=costheta,
                        epsilon=ep,
                        nu=.0))

# change Poisson's ratio
nus = [.0, .25, .5]
b2s = []
for nu in nus:
    b2s.append(boundary(costheta=costheta,
                        epsilon=.1,
                        nu=nu))

# plot
plt.figure(figsize=(8, 4))

ax1 = plt.subplot(121, projection="polar")
for ep, bi in zip(eps, b1s):
    ax1.plot(theta, bi, label="ϵ={:.2f}, ν=0".format(ep))
ax1.legend()

ax2 = plt.subplot(122, projection="polar")
for nu, bi in zip(nus, b2s):
    ax2.plot(theta, bi, label="ϵ=.1, ν={:.1f}".format(nu))
ax2.legend()

plt.tight_layout()
plt.show()
