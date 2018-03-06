============
Theory Notes
============

.. toctree::
  :maxdepth: 2


Optical stretcher
-----------------
The optical stretcher consists of a dual beam laser trap, in its original
configuration built from to opposing optical fibers :cite:`Guck2001`.
When increasing the trapping power, compliant objects such as cells
are stretched along the axis of the trap. Using video analysis, the
measured shape change can be translated into physical properties of the
cell.


Global geometric factor
-----------------------
The global geometric factor (GGF) connects compliance (how easy it is
to deform a body consisting of a certain material) and strain (how
much this body is deformed). Thus, the GGF is a measure of stress
(force acting on the surface of the body).

In an optical stretcher (OS) experiment, the strain :math:`\epsilon`
of a cell can be measured by analyzing its deformation (e.g. via a
contour in the intensity image). Using cell size and the measured
change in eccentricity, as well as several parameters of the OS
setup itself, :mod:`ggf` can be used to compute the optical stress
:math:`\sigma_0` and the geometric factor :math:`F_\text{G}`. Together,
they resemble the GGF

.. math::
    
    \text{GGF} = \sigma_0 F_\text{G}.

The compliance of the cell then computes to

.. math::

    J = \frac{\epsilon}{\text{GGF}}
