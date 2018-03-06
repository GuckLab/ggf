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


How should I migrate my Matlab pipeline to Python?
--------------------------------------------------
You can access the computations performed in *StretcherNStress.m* via
:func:`ggf.core.stress`.

.. code::

    from ggf.core import stress
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

    from ggf.globgeomfact import coeff2ggf
    GGF = coeff2ggf(coeff, poisson_ratio=.45)

These methods produce the same output as the original Matlab scripts
with an accuracy that is below the standard tolerance of :func:`numpy.allclose`.
