"""Creep compliance analysis

This example uses the contour data of an HL60 cell in the OS to
compute its GGF and creep compliance. The `contour data
<_static/creep_compliance_data.h5>`__ were determined
from `this phase-contrast video <_static/creep_compliance.mp4>`__
(prior to video compression). During stretching, the total laser
power was increased from 0.2W to 1.3W (reflexes due to second
harmonic effects appear as white spots).
"""
import ggf
import h5py
import lmfit
import matplotlib.pylab as plt
import numpy as np


def ellipse_fit(radius, theta):
    '''Fit a centered ellipse to data in polar coordinates

    Parameters
    ----------
    radius: 1d ndarray
        radial coordinates
    theta: 1d ndarray
        angular coordinates [rad]

    Returns
    -------
    a, b: floats
        semi-axes of the ellipse; `a` is aligned with theta=0.
    '''
    def residuals(params, radius, theta):
        a = params["a"].value
        b = params["b"].value
        r = a*b / np.sqrt(a**2 * np.sin(theta)**2 + b**2 * np.cos(theta)**2)
        return r - radius

    parms = lmfit.Parameters()
    parms.add(name="a", value=radius.mean())
    parms.add(name="b", value=radius.mean())

    res = lmfit.minimize(residuals, parms, args=(radius, theta))

    return res.params["a"].value, res.params["b"].value


# load the contour data (stored in polar coordinates)
with h5py.File("data/creep_compliance_data.h5", "r") as h5:
    radius = h5["radius"].value * 1e-6  # [µm] to [m]
    theta = h5["theta"].value
    time = h5["time"].value
    meta = dict(h5.attrs)


factors = np.zeros(len(radius), dtype=float)
semimaj = np.zeros(len(radius), dtype=float)
semimin = np.zeros(len(radius), dtype=float)
strains = np.zeros(len(radius), dtype=float)
complnc = np.zeros(len(radius), dtype=float)

for ii in range(len(radius)):
    # determine semi-major and semi-minor axes
    smaj, smin = ellipse_fit(radius[ii], theta[ii])
    semimaj[ii] = smaj
    semimin[ii] = smin
    # compute GGF
    if (time[ii] > meta["time_stretch_begin [s]"]
            and time[ii] < meta["time_stretch_end [s]"]):
        power_per_fiber = meta["power_per_fiber_stretch [W]"]
        f = ggf.get_ggf(
            model="boyde2009",
            semi_major=smaj,
            semi_minor=smin,
            object_index=meta["object_index"],
            medium_index=meta["medium_index"],
            effective_fiber_distance=meta["effective_fiber_distance [m]"],
            mode_field_diameter=meta["mode_field_diameter [m]"],
            power_per_fiber=power_per_fiber,
            wavelength=meta["wavelength [m]"],
            poisson_ratio=.5)
    else:
        power_per_fiber = meta["power_per_fiber_trap [W]"]
        f = np.nan
    factors[ii] = f

# compute compliance
strains = (semimaj-semimaj[0]) / semimaj[0]
complnc = strains / factors
compl_ival = (time > meta["time_stretch_begin [s]"]) * \
    (time < meta["time_stretch_end [s]"])
stretch_index = np.where(compl_ival)[0][0]
complnc_1 = strains/factors[stretch_index]

# plots
plt.figure(figsize=(8, 7))

ax1 = plt.subplot(221, title="ellipse fit semi-axes")
ax1.plot(time, semimaj*1e6, label="semi-major axis")
ax1.plot(time, semimin*1e6, label="semi-minor axis")
ax1.legend()
ax1.set_xlabel("time [s]")
ax1.set_ylabel("axis radius [µm]")

ax2 = plt.subplot(222, title="GGF")
ax2.plot(time, factors)
ax2.set_xlabel("time [s]")
ax2.set_ylabel("global geometric factor [Pa]")

ax3 = plt.subplot(223, title="strain")
ax3.plot(time, (strains)*100)
ax3.set_xlabel("time [s]")
ax3.set_ylabel("deformation $w(t)/r_0$ [%]")

ax4 = plt.subplot(224, title="creep compliance")
ax4.plot(time[compl_ival], complnc[compl_ival])
ax4.set_xlabel("time [s]")
ax4.set_ylabel("compliance $J(t)$ [Pa⁻¹]")

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlim(0, np.round(time.max()))
    ax.axvline(x=meta["time_stretch_begin [s]"], c="r")
    ax.axvline(x=meta["time_stretch_end [s]"], c="r")

plt.tight_layout()
plt.show()
