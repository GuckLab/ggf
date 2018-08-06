import numpy as np
import pathlib

from .core import legendre2ggf, stress2ggf, stress2legendre  # noqa: F401
from . import lut
from . import stress
from .stress.geometry import fiber_distance_capillary  # noqa: F401


def get_ggf(model, semi_major, semi_minor, object_index, medium_index,
            effective_fiber_distance=100e-6, mode_field_diameter=3e-6,
            power_per_fiber=.6, wavelength=1064e-9, poisson_ratio=0.5,
            n_poly=120, use_lut=None, verbose=False):
    """Model the global geometric factor

    Parameters
    ----------
    model: str
        Model to use, one of: `boyde2009`
    semi_major: float
        Semi-major axis of an ellipse fit to the object perimeter [m]
    semi_minor: float
        Semi-minor axis of an ellipse fit to the object perimeter [m]
    object_index: float
        Refractive index of the object
    medium_index: float
        Refractive index of the surrounding medium
    effective_fiber_distance: float
        Effective distance between the two trapping fibers relative
        to the medium refractive index [m]. For an open setup, this is
        the physical distance between the fibers. For a closed setup
        (capillary), this distance takes into account the refractive
        indices and thicknesses of the glass capillary and index
        matching gel. For the closed setup, the convenience function
        :func:`ggf.fiber_distance_capillary` can be used.
    mode_field_diameter: float
        The mode field diameter MFD of the fiber used [m]. Note that
        the MFD is dependent on the wavelength used. If the
        manufacturer did not provide a value for the MFD, the MFD
        can be approximated as ``3*wavelenth`` for a single-mode
        fiber.
    power_per_fiber: float
        The laser power coupled into each of the fibers [W]
    wavelength: float
        The laser wavelength used for the trap [m]
    poisson_ratio: float
        The Poisson's ratio of the stretched material. Set this
        to 0.5 for volume conservation.
    n_poly: int
        Number of Legendre polynomials to use for computing the GGF.
        Note that only even Legendre polynomials are used and thus,
        this number is effectively halved. To reproduce the GGF as
        computed with the Boyde2009 Matlab script, set this value
        to `None`.
    use_lut: None, str, pathlib.Path or bool
        Use look-up tables to compute the GGF. If set to `None`,
        the internal LUTs will be used or the GGF is computed if
        it cannot be found in a LUT. If `True`, the internal LUTs
        will be used and a `NotInLUTError` will be raised if the
        GGF cannot be found there. Alternatively, a path to a
        LUT hdf5 file can be passed.
    verbose: int
        Increases verbosity

    Returns
    -------
    ggf: float
        global geometric factor
    """
    if use_lut not in [False]:
        try:
            if isinstance(use_lut, (str, pathlib.Path)):
                lut_path = use_lut
            else:
                lut_path = None
            ggf = lut.get_ggf_lut(
                model=model,
                semi_major=semi_major,
                semi_minor=semi_minor,
                object_index=object_index,
                medium_index=medium_index,
                effective_fiber_distance=effective_fiber_distance,
                mode_field_diameter=mode_field_diameter,
                power_per_fiber=power_per_fiber,
                wavelength=wavelength,
                poisson_ratio=poisson_ratio,
                n_poly=n_poly,
                lut_path=lut_path,
                verbose=verbose)
        except lut.NotInLUTError:
            if use_lut:  # user specifically defined to ONLY use LUT
                raise
            else:  # force computation of the GGF
                ggf = get_ggf(
                    model=model,
                    semi_major=semi_major,
                    semi_minor=semi_minor,
                    object_index=object_index,
                    medium_index=medium_index,
                    effective_fiber_distance=effective_fiber_distance,
                    mode_field_diameter=mode_field_diameter,
                    power_per_fiber=power_per_fiber,
                    wavelength=wavelength,
                    poisson_ratio=poisson_ratio,
                    n_poly=n_poly,
                    use_lut=False,
                    verbose=verbose)
    else:
        theta, sigma = stress.get_stress(
            model=model,
            semi_major=semi_major,
            semi_minor=semi_minor,
            object_index=object_index,
            medium_index=medium_index,
            effective_fiber_distance=effective_fiber_distance,
            mode_field_diameter=mode_field_diameter,
            power_per_fiber=power_per_fiber,
            wavelength=wavelength,
            verbose=verbose)

        if n_poly is None:
            # number of orders (estimate from Boyde 2009)
            alpha = semi_minor * 2 * np.pi / wavelength  # size parameter
            n_poly = np.int(np.round(2 + alpha + 4 * (alpha)**(1 / 3) + 10))
        else:
            n_poly = int(np.round(n_poly))

        ggf = stress2ggf(stress=sigma,
                         theta=theta,
                         poisson_ratio=poisson_ratio,
                         n_poly=n_poly)

    return ggf
