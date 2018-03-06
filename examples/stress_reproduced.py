"""Radial stresses of a prolate spheroid

This examples computes radial stress profiles for spheroidal
objects in the optical stretcher, reproducing figures (9) and (10)
of :cite:`Boyde2009`.
"""
import matplotlib.pylab as plt
import numpy as np
import percache

from ggf.core import stress


@percache.Cache("stress_reproduced.cache", livesync=True)
def compute(**kwargs):
    "Locally cached version of ggf.core.stress"
    return stress(**kwargs)


# variables from the publication
alpha = 47
wavelength = 1064e-9
radius = alpha * wavelength / (2 * np.pi)

kwargs = {"stretch_ratio": .1,
          "object_index": 1.375,
          "medium_index": 1.335,
          "wavelength": wavelength,
          "beam_waist": 3,
          "radius": radius,
          "power_left": 1,
          "power_right": 1,
          "poisson_ratio": 0,
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
kwargs5["dist"] = 120e-6

kwargs6 = kwargs.copy()
kwargs6["dist"] = 200e-6


# polar plots
plt.figure(figsize=(8, 5))

th1, sigma1 = compute(**kwargs1)
ax1 = plt.subplot(231, projection='polar')
ax1.plot(th1, sigma1, "k")
ax1.plot(th1 + np.pi, sigma1[::-1], "k")

th2, sigma2 = compute(**kwargs2)
ax2 = plt.subplot(232, projection='polar')
ax2.plot(th2, sigma2, "k")
ax2.plot(th2 + np.pi, sigma2[::-1], "k")

th3, sigma3 = compute(**kwargs3)
ax3 = plt.subplot(233, projection='polar')
ax3.plot(th3, sigma3, "k")
ax3.plot(th3 + np.pi, sigma3[::-1], "k")

for ax in [ax1, ax2, ax3]:
    ax.set_rticks([0, 1.5, 3, 4.5])
    ax.set_rlim(0, 4.5)

th4, sigma4 = compute(**kwargs4)
ax4 = plt.subplot(234, projection='polar')
ax4.plot(th4, sigma4, "k")
ax4.plot(th4 + np.pi, sigma4[::-1], "k")
ax4.set_rticks([0, 4, 8, 12])
ax4.set_rlim(0, 12)

th5, sigma5 = compute(**kwargs5)
ax5 = plt.subplot(235, projection='polar')
ax5.plot(th5, sigma5, "k")
ax5.plot(th5 + np.pi, sigma5[::-1], "k")
ax5.set_rticks([0, 1.5, 3, 4.5])
ax5.set_rlim(0, 4.5)

th6, sigma6 = compute(**kwargs6)
ax6 = plt.subplot(236, projection='polar')
ax6.plot(th6, sigma6, "k")
ax6.plot(th6 + np.pi, sigma6[::-1], "k")
ax6.set_rticks([0, 0.6, 1.2, 1.8])
ax6.set_rlim(0, 1.8)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_thetagrids(np.linspace(0, 360, 12, endpoint=False))

plt.tight_layout()
plt.show()
