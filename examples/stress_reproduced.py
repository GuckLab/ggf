"""Radial stresses of a prolate spheroid

:cite:`Boyde2009`
"""
import matplotlib.pylab as plt
import numpy as np
import percache

from ggf.core import stress

alpha = 47
apbp = 1.1
radius = alpha*1064e-9 / (2*np.pi)

stretch_ratio = 0.1
poisson_ratio= 0.0


kwargs = {"stretch_ratio": stretch_ratio,
          "object_index": 1.375,
          "medium_index": 1.335,
          "wavelength": 1064e-9,
          "beam_waist": 3 * 1064e-9,
          "dist": 120e-6,
          "radius": radius,
          "power_left":1,  # [W]
          "power_right":1,  # [W]
          "poisson_ratio": poisson_ratio,
          "numpoints": 200,
          }


kwargs1 = kwargs.copy()
kwargs1["power_right"] = 0
kwargs1["stretch_ratio"] = 0
kwargs1["dist"] = 90e-6

kwargs2 = kwargs.copy()
kwargs2["power_right"] = 0
kwargs2["stretch_ratio"] = .05
kwargs2["dist"] = 90e-6

kwargs3 = kwargs.copy()
kwargs3["power_right"] = 0
kwargs3["stretch_ratio"] = .1
kwargs3["dist"] = 90e-6

kwargs4 = kwargs.copy()
kwargs4["dist"] = 60e-6

kwargs5 = kwargs.copy()

kwargs6 = kwargs.copy()
kwargs6["dist"] = 200e-6


@percache.Cache("cache_test", livesync=True)
def compute(**kwargs):
    return stress(**kwargs)

plt.figure(figsize=(15,7))

th1, sigma1 = compute(**kwargs1)
ax1 = plt.subplot(231, projection='polar')
ax1.plot(th1, sigma1)

th2, sigma2 = compute(**kwargs2)
ax2 = plt.subplot(232, projection='polar')
ax2.plot(th2, sigma2)

th3, sigma3 = compute(**kwargs3)
ax3 = plt.subplot(233, projection='polar')
ax3.plot(th3, sigma3)

for ax in [ax1, ax2, ax3]:
    ax.set_rticks([0, 1.5, 3, 4.5])

th4, sigma4 = compute(**kwargs4)
ax4 = plt.subplot(234, projection='polar')
ax4.plot(th4, sigma4)
ax4.set_rticks([0, 4, 8, 12])

th5, sigma5 = compute(**kwargs5)
ax5 = plt.subplot(235, projection='polar')
ax5.plot(th5, sigma5)
ax5.set_rticks([0, 1.5, 3, 4.5])

th6, sigma6 = compute(**kwargs6)
ax6 = plt.subplot(236, projection='polar')
ax6.plot(th6, sigma6)
ax6.set_rticks([0, 0.6, 1.2, 1.8])

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_thetagrids(np.linspace(0, 360, 12, endpoint=False))

plt.show()
