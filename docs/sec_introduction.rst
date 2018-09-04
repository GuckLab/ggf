============
Introduction
============

.. toctree::
  :maxdepth: 2


What is the package "ggf" used for?
-----------------------------------
It is a Python implementation of two Matlab scripts by
Lars Boyde, *StretcherNStress.m* and *GGF.m*, which are used in
the Guck lab to compute optical stress distributions and resulting
global geometric factors for dielectric, elastic, and spheroidal
objects in the optical stretcher.


What is an optical stretcher?
-----------------------------
The optical stretcher consists of a dual beam laser trap, in its original
configuration built from two opposing optical fibers :cite:`Guck2001`.
When increasing the trapping power, compliant objects such as cells
are stretched along the axis of the trap. Using video analysis, the
measured shape change can be translated into physical properties of the
object.


What is the global geometric factor?
------------------------------------
The global geometric factor (GGF) connects (the unknown variable)
compliance :math:`J` (how easy it is to deform a body consisting of a certain
material) and (the measured variable) strain :math:`\epsilon` (how much this
body is deformed). Thus, the GGF is a measure of stress (force acting on the
surface of the object).

.. math::

    J = \frac{\epsilon}{\text{GGF}}

In an optical stretcher (OS) experiment, the strain :math:`\epsilon`
of an object can be measured by analyzing its deformation (e.g. via a
contour in the intensity image). Using object radius and the measured
change in eccentricity, as well as several parameters of the OS
setup itself, :mod:`ggf` can be used to compute the optical stress
:math:`\sigma` from which the GGF is computed.


How should I migrate my Matlab pipeline to Python?
--------------------------------------------------

To reproduce data
.................
You can access the computations performed in *StretcherNStress.m* via
:func:`ggf.stress.boyde2009.core.stress`.

.. warning:: There was a mistake in the original boundary function (see
    `issue #1 <https://github.com/GuckLab/ggf/issues/1>`__). This affects
    all cases where ``poisson_ratio`` is non-zero. If you would like
    to reproduce exactly the stress profiles of *StretcherNStress.m*,
    please use ggf version 0.2.0.
    
.. code::

    from ggf.stress.boyde2009.core import stress
    theta, sigma, coeff = stress(object_index=1.41,
                                 medium_index=1.3465,
                                 semi_minor=2.8466e-6,   # [m]
                                 poisson_ratio=0.45,
                                 stretch_ratio=0.1,
                                 wavelength=780e-9,      # [m]
                                 beam_waist=3.08,        # [wavelengths]
                                 power_left=.65,         # [W]
                                 power_right=.65,        # [W]
                                 dist=180e-6,            # [m]
                                 field_approx="davis",
                                 ret_legendre_decomp=True)

The GGF can be computed from the coefficients ``coeff`` via
:func:`ggf.legendre2ggf`.

.. code::

    from ggf import legendre2ggf
    GGF = legendre2ggf(coeff, poisson_ratio=.45)

These methods produce the same output as the original Matlab scripts
with an accuracy that is below the standard tolerance of :func:`numpy.allclose`.

For a new project
.................
In general, the method :func:`ggf.get_ggf` is recommended. The difference
to the above method is:

- It makes use of precomputed look-up tables (LUTs) which avoids long
  computation times. The error made by using LUTs maxes at about 1-2%.
- It comes with user-convenient keyword arguments.
- The GGF is computed from 120 Legendre coefficients by default, a number
  that was previously determined automatically and could have potentially
  been too low.

Please note that due to these points, the resulting GGF might
vary from the GGF computed with the original Matlab script.

.. code::

    import ggf
    GGF = ggf.get_ggf(model="boyde2009",
                      semi_major=3.2788e-6,             # [m]
                      semi_minor=2.8466e-6,             # [m]
                      object_index=1.41,
                      medium_index=1.3465,
                      effective_fiber_distance=180e-6,  # [m]
                      mode_field_diameter=4.8e-06,      # [m]
                      power_per_fiber=.65,              # [W]
                      wavelength=780e-9,                # [m]
                      poisson_ratio=0.45)
