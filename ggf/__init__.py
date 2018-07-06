import numpy as np

from .core import legendre2ggf, stress2ggf, stress2legendre
from . import stress
from .stress.geometry import fiber_distance_capillary


def get_ggf(model, semi_major, semi_minor, object_index, medium_index,
            effective_fiber_distance=100e-6, mode_field_diameter=3e-6,
            power_per_fiber=.6, wavelength=1064e-9, poisson_ratio=0.5,
            verbose=False):
    """Model the global geometric factor"""
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

    # number of orders (estimate from Boyde 2009)
    alpha = semi_minor * 2 * np.pi / wavelength  # size parameter
    n_poly = np.int(np.round(2+alpha+4*(alpha)**(1/3) + 10))

    ggf = stress2ggf(stress=sigma, theta=theta, poisson_ratio=poisson_ratio,
                     n_poly=n_poly)

    return ggf
