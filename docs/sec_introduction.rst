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
global geometric factors for spherical and spheroidal objects
in the optical stretcher.


What is an optical stretcher?
-----------------------------
The optical stretcher consists of a dual beam laser trap, in its original
configuration built from two opposing optical fibers :cite:`Guck2001`.
When increasing the trapping power, compliant objects such as cells
are stretched along the axis of the trap. Using video analysis, the
measured shape change can be translated into physical properties of the
cell.


What is the global geometric factor?
------------------------------------
The global geometric factor (GGF) connects (the unknown variable)
compliance :math:`J` (how easy it is to deform a body consisting of a certain
material) and (the measured variable) strain :math:`\epsilon` (how much this
body is deformed). Thus, the GGF is a measure of stress (force acting on the
surface of the body).

.. math::

    J = \frac{\epsilon}{\text{GGF}}

In an optical stretcher (OS) experiment, the strain :math:`\epsilon`
of a cell can be measured by analyzing its deformation (e.g. via a
contour in the intensity image). Using cell size and the measured
change in eccentricity, as well as several parameters of the OS
setup itself, :mod:`ggf` can be used to compute the optical stress
:math:`\sigma` from which the GGF is computed.


How should I migrate my Matlab pipeline to Python?
--------------------------------------------------
You can access the computations performed in *StretcherNStress.m* via
:func:`ggf.core.stress`.

.. code::

    from ggf.stress.boyde2009 import stress
    theta, sigma, coeff = stress(object_index=1.41,
                                 medium_index=1.3465,
                                 radius=2.8466e-6,    # [m]
                                 poisson_ratio=0.45,
                                 stretch_ratio=0.1,
                                 wavelength=780e-9,   # [m]
                                 beam_waist=3,        # [wavelengths]
                                 power_left=.6,       # [W]
                                 power_right=.6,      # [W]
                                 dist = 100e-6,       # [m]
                                 numpoints=100,
                                 theta_max=np.pi,
                                 field_approx="davis",
                                 ret_legendre_decomp=True)

The GGF can be computed from the coefficients ``coeff`` via
:func:`ggf.globgeomfact.coeff2ggf`.

.. code::

    from ggf import legendre2ggf
    GGF = coeff2ggf(coeff, poisson_ratio=.45)

These methods produce the same output as the original Matlab scripts
with an accuracy that is below the standard tolerance of :func:`numpy.allclose`.
