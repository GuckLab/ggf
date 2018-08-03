"""GGF as a function of beam waist radius

Note: This example is not reproducing a previously published figure, which
could be due to false simulation parameters.

This example attempts reproduce figure 7 in :cite:`Bellini2012`.
Two implementations of the optical stretcher are described, the assembled
optical stretcher (AOS) and the monolithic optical stretcher (MOS).
In the publication, these have two different mode field diameters (MFD),
6.2µm (AOS) and 9µm (MOS).
The MFD `is equivalent to twice the gaussian beam waist radius
<https://www.rp-photonics.com/mode_radius.html>`_.
In the original publication, the GGF was plotted over the beam waist
diameter (2w₀) and the AOS and MOS mode field diameters were indicated in the
plots. Here, the data are plotted over the beam waist radius (w₀) and half
of the mode field diameters are indicated in the plot. In comparison to the
original figure, it appears as if a factor of two is missing in
the definition of the plotted beam waist radii.
"""
import matplotlib.pylab as plt
import numpy as np
import percache

from ggf.stress.boyde2009 import stress
from ggf import legendre2ggf
from ggf.stress.boyde2009.stretcher import distance_capillary

medium_index = 1.338

dist1 = distance_capillary(gel_dist=5e-6,
                           glass_dist=40e-6,
                           medium_dist=40e-6,
                           gel_index=1.449,
                           glass_index=1.474,
                           medium_index=medium_index)

dist2 = distance_capillary(gel_dist=0,
                           glass_dist=25e-6,
                           medium_dist=50e-6,
                           gel_index=1.449,
                           glass_index=1.474,
                           medium_index=medium_index)

kwargs = {"object_index": 1.37108,
          "medium_index": medium_index,
          "poisson_ratio": 0,
          "semi_minor": 7.51e-6,
          "stretch_ratio": .041,
          "wavelength": 1064e-9,
          "power_left": 1.5,
          "power_right": 1.5,
          "n_points": 181,
          }

kwargs_aos = kwargs.copy()
kwargs_aos.update({"beam_waist": 3,
                   "dist": dist1,
                   })

kwargs_mos = kwargs.copy()
kwargs_mos.update({"beam_waist": 3,
                   "dist": dist2,
                   })


@percache.Cache("ggf_reproduced.cache", livesync=True)
def compute_ggf(**kw):
    coeff = stress(ret_legendre_decomp=True, **kw)[0]
    return legendre2ggf(coeff, poisson_ratio=kw["poisson_ratio"])


aos_ggf = []
mos_ggf = []

w0s = np.linspace(3e-6, 10e-6, 30)
for w0 in w0s:
    kwargs_aos["beam_waist"] = w0/kwargs["wavelength"]
    aos_ggf.append(compute_ggf(**kwargs_aos))
    kwargs_mos["beam_waist"] = w0/kwargs["wavelength"]
    mos_ggf.append(compute_ggf(**kwargs_mos))

plt.figure(figsize=(8, 4))
ax = plt.subplot(111)
ax.plot(w0s[:len(aos_ggf)]*1e6, aos_ggf,
        marker="^", color="red", label="AOS", lw=2)
ax.plot(w0s[:len(aos_ggf)]*1e6, mos_ggf,
        marker="d", color="blue", label="MOS", lw=2)
ax.axvline(6.2/2, color="red", ls=":", label="AOS beam waist radius")
ax.axvline(9/2, color="blue", ls=":", label="MOS beam waist radius")
ax.legend()
ax.set_ylabel("GGF [Pa]")
ax.set_xlabel("gaussian beam waist radius w₀ [µm]")
ax.grid()

plt.tight_layout()
plt.show()
