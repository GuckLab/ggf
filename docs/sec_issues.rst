Known issues
============

Accuracy of the mode field diameter
-----------------------------------
The mode field diameter (MFD) is an important parameter for the computation
of the GGF (see :func:`ggf.get_ggf`). Manufacturers often list the mode field
diameters with large error margins and only for a single wavelength
(e.g. 5.0 ± 0.5 μm @ 850nm for a 
`THORLABS 780HP <https://www.thorlabs.de/thorProduct.cfm?partNumber=780HP>`_).
Thus, an accurate value of the MFD is not given, especially for other wavelengths.

According to thorlabs (personal communication), the
*MFD is not a directly measured value, but a function of wavelength, core radius and the refractive indices of the core and the cladding. The measurement of MFD is accomplished by the Variable Aperture Method in the Far Field (VAMFF). An aperture is placed in the far field of the fiber output, and the intensity is measured. As successively smaller apertures are placed in the beam, the intensity levels are measured for each aperture; the data can then be plotted as power vs. the sine of the aperture half-angle (or the numerical aperture). The MFD is then determined using Petermann's second definition, which is a mathematical model that does not assume a specific shape of power distribution. The MFD in the near field can be determined from this far-field measurement using the Hankel Transform.*

Thus, there are several sources of error that are propagated to the MFD.
It would probably be best to directly measure the MFD and investigate how its
measurement error propagates to the GGF. To our knowledge, this has not yet been done.


Method-specific differences
---------------------------
For some parameter combinations, the methods in :cite:`Boyde2009` and :cite:`Boyde2012`
yield very different stress profile shapes (data not shown). It has not yet been investigated how
these differences affect the GGF and whether they can be explained by the fact
that the generalized Lorentz-Mie theory approach is simply more accurate.

