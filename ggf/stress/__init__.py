from . import boyde2009
from .geometry import fiber_distance_capillary  # noqa: F401

VALID_MODELS = ["boyde2009"]


def get_stress(model, semi_major, semi_minor, object_index, medium_index,
               effective_fiber_distance=100e-6, mode_field_diameter=3e-6,
               power_per_fiber=.6, wavelength=1064e-9, n_points=100,
               verbose=False):
    """Compute the optical stress profile in the optical stretcher

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
    n_points: int
        Number of points to compute.
    verbose: int
        Increases verbosity

    Returns
    -------
    theta: 1d ndarray of length `n_points`
        Polar angles [rad]
    sigma: 1d ndarray of length `n_points`
        Radial stress profile along `theta` [Pa]
    """
    if model not in VALID_MODELS:
        msg = "`model` must be one of {}, got '{}'".format(VALID_MODELS, model)
        raise ValueError(msg)

    if model == "boyde2009":
        func = boyde2009.get_stress

    return func(
        semi_major=semi_major,
        semi_minor=semi_minor,
        object_index=object_index,
        medium_index=medium_index,
        effective_fiber_distance=effective_fiber_distance,
        mode_field_diameter=mode_field_diameter,
        power_per_fiber=power_per_fiber,
        wavelength=wavelength,
        n_points=n_points,
        verbose=verbose)
