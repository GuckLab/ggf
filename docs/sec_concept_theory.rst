==================
Concept and theory
==================

.. toctree::
  :maxdepth: 2


Summary
=======
The computation of the compliance :math:`J` for dielectric, elastic,
spheroidal objects in the OS can be divided into three main tasks:
measuring the deformation :math:`w`, modeling the optical stress
:math:`\sigma_r`, and computing the GGF from the stress.
Several approaches to these problems have been presented in the
related literature and are discussed in the following.


Experimentally quantifying deformation
======================================
The deformation is quantified from video images by fitting an ellipse
to the contour of the stretched object. The deformation can then be
defined as a relation between the semimajor or semiminor axes and the
initial radius (or semiaxes). Note that the prolate spheroidal shape
is only an approximation to the actual shape. The methods in this package
more or less assume that this approximation is valid.


Optical stress profile acting on a prolate spheroid
===================================================
The optical stress :math:`\sigma(\theta)` in dependence of the
angle :math:`\theta` is a result of the optical forces acting on the
surface of the spheroid. The angle :math:`\theta` is defined in the
imaging plane in a typical OS experiment, with :math:`\theta=0`
pointing to the right hand fiber.


:math:`\cos^2\theta` approximation
----------------------------------
Ray optics is used to compute the optical stress acting on a
spheroid and a :math:`\sigma_0 cos^2\theta` model is fitted to the
resulting stress profile with the peak stress
:math:`\sigma_0` :cite:`Guck2001`. The :math:`\sigma_0 cos^2\theta`
approximation simplifies subsequent computations. 

Note that a more general model :math:`\sigma_0 cos^2n\theta` with
larger exponents (e.g. :math:`n` = 2,3,4,...) can also be applied,
e.g. for different fibroblast cell lines :cite:`Ananthakrishnan2006`. 


Semi-analytical perturbation approach (Boyde et al. 2009)
---------------------------------------------------------
- gaussian laser beam
- :math:`a > \lambda`: higher order perturbation theory
- :cite:`Boyde2009`


Generalized Lorentz-Mie theory (Boyde et al. 2012)
--------------------------------------------------
- gaussian laser beam
- spheroidal coordinates
- generalized Lorenz–Mie theory
- not implemented (Matlab sources available upon request)
- :cite:`Boyde2012`


.. _sec_theory_ggf:

Computation of the GGF
======================
The following derivations are based on the theoretical considerations
of Lur'e :cite:`Lure1964` for a rotationally symmetric deformation of
a sphere (which in general does not result in prolate spheroids)
and their application to the OS by Ananthakrishnan
et al. :cite:`Ananthakrishnan2006`. Note that a corrigendum has been published
for this article in 2008 :cite:`Ananthakrishnan2008Corr`.

General approach
----------------
The GGF connects the measured deformation to the shear modulus :math:`G`
which, in OS literature, is usually written in the form

.. math::

    \frac{w}{r_0} = \frac{\text{GGF}}{G}

where :math:`w` is the change in radius of the stretched sphere along
the stretcher axis and :math:`r_0` is the radius of the unstretched sphere.
Note that the quantity :math:`w/r_0` resembles a measure of strain along
the stretcher axis.

The GGF can be computed from the radial stress :math:`\sigma_r(\theta)`
via the radial displacement :math:`u_r(r, \theta)`. These quantities can be
connected via a Legendre decomposition according to
(:cite:`Lure1964`, chapter 6)

.. math::

    u_r(r, \theta) &= \sum_n \left[ A_n r^{n+1} (n+1)(n-2+4\nu)
                               + B_n r^{n-1} n  \right] P_n(\cos \theta)


    \frac{\sigma_r(r, \theta)}{2G} &= \sum_n \left[ A_n r^n (n+1) (n^2 - n - 2 - 2\nu)
                                       + B_n r^{n-2} n (n-1)
                                 \right] P_n(\cos \theta) 

with the Legendre polynomials :math:`P_n` and the Poisson's ratio :math:`\nu`.
The coefficients :math:`A_n` and :math:`B_n` have to be determined from
boundary conditions. For the case of normal loading, which is given by
the electromagnetic boundary conditions in the OS
(:math:`\sigma_\theta=\tau_{r,\theta}=0`), these coefficients compute to:

.. math::

    A_0 = - \frac{s_0}{4G(1+\nu)}
    
    B_0 = A_1 = B_1 = 0

and for :math:`n>=2`:

.. math::

    A_n &= - \frac{s_n}{4Gr_0^n \Delta}

    B_n &= \frac{s_n}{4Gr_0^{n-2} \Delta} \cdot \frac{n^2 + 2n -1 + 2\nu}{n-1}

  \text{with }  \Delta &= n(n-1) + (2n+1) (\nu + 1)

Where :math:`s_n` is the :math:`n\text{th}` component of the Legendre
decomposition of :math:`\sigma_r`

.. math::

    \sigma_r(\theta) = \sum_n s_n P_n(\cos \theta).

The radial displacement then takes the form

.. math::

    u_r(r, \theta) = \frac{r_0}{G} \left[
                     \frac{(1-2\nu) s_0}{2(1+\nu)} +
                     \sum_{n=2}^\infty
                     \frac{2s_n}{2n+1}
                     \left(L_n \left(\frac{r}{r_0}\right)^n + 
                           M_n \left(\frac{r}{r_0}\right)^{n-2} \right)
                     P_n(\cos \theta)
                     \right]

with the coefficients :math:`L_n` and :math:`M_n` given in 
:cite:`Lure1964`, chapter 6.6.
We measure the displacement at the outer perimeter of the stretched sphere
and on the stretcher axis only; Thus, we set :math:`r=r_0` and
:math:`\theta=0` with :math:`w=u_r(r_0, 0)`.

To obtain the GGF, we finally compute

.. math::

    \text{GGF} &= \frac{G}{r_0} u_r(r_0, 0)

        &= \left[
           \frac{(1-2\nu) s_0}{2(1+\nu)} +
           \sum_{n=2}^\infty
           \frac{2s_n}{2n+1}
           \left(L_n + M_n \right)
           P_n(\cos \theta)
           \right].


Notes:

- Due to the fact hat the stress profile in the OS is rotationally
  symmetric w.r.t. the stretcher axis, all odd coefficients
  :math:`s_n` are zero.

- The polar displacement :math:`u_\theta` has been omitted here,
  because it does not represent a quantity measurable in an OS
  experiment.


Special case: :math:`\cos^2\theta` approximation
------------------------------------------------
Following the above approach, the stress profile

.. math::

    \sigma_r(\theta) = \sigma_0 \cos^2\theta

with the peak stress :math:`\sigma_0` can be decomposed into
two Legendre polynomials

.. math:: 

    \sigma_r(\theta) &= s_0 P_0(\cos\theta) + s_2 P_2(\cos\theta)

    s_0 &= \frac{1}{3} \sigma_0

    s_2 &= \frac{2}{3} \sigma_0

Inserting these Legendre coefficients in the above equation for the
GGF yields

.. math::

   \text{GGF} = \frac{\sigma_0}{2(1+\nu)} \left[
         \frac{1}{3} \left( (1-2\nu) + \frac{(-7 + 4\nu)(1+\nu)}{7+5\nu} \right)
         + \frac{(7-4\nu)(1+\nu)}{7+5\nu} \cos^2\theta \right].
   
Historically, the relation between strain, stress, and shear modulus
was written in the form

.. math::
    
    \frac{w}{r_0} = \frac{\sigma_0 F_\text{G}}{G}

with the geometrical factor :math:`F_\text{G}` that does not include
the peak stress :math:`\sigma_0`. Hence the term "global geometrical factor"
:math:`\text{GGF} = \sigma_0 F_\text{G}`.


Computation of compliance
=========================
A typical OS experiment records the deformation :math:`w(t)` over time
:math:`t`. The quantity of interest is the (creep) compliance :math:`J(t)`.
With :math:`J = 1/G`, it computes to

.. math::

    J(t) = \frac{w(t)}{r_0} \cdot \frac{1}{\text{GGF}(t)}. 

Note that the GGF is now time-dependent, because the optical stress
profile :math:`\sigma_r`, from which the GGF is computed, also depends
on the deformation. 
