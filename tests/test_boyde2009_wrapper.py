import numpy as np

import ggf
from ggf.stress import boyde2009


def test_basic():
    """Test that all parameters are passed correctly"""
    radius = 2.8466e-6
    stretch_ratio = 0.1
    poisson_ratio = 0

    semi_minor = radius
    semi_major = radius * (1 + stretch_ratio)

    object_index = 1.3776
    medium_index = 1.335
    effective_fiber_distance = 100e-6
    mode_field_diameter = 4.8e-6
    wavelength = 800e-9
    power_per_fiber = .6

    n_points = 10

    th1, sg1 = boyde2009.core.stress(
        semi_minor=radius,
        stretch_ratio=stretch_ratio,
        poisson_ratio=poisson_ratio,
        object_index=object_index,
        medium_index=medium_index,
        dist=effective_fiber_distance/2,
        wavelength=wavelength,
        beam_waist=mode_field_diameter/2/wavelength,
        power_left=power_per_fiber,
        power_right=power_per_fiber,
        n_points=n_points)

    th2, sg2 = ggf.stress.get_stress(
        model="boyde2009",
        semi_major=semi_major,
        semi_minor=semi_minor,
        object_index=object_index,
        medium_index=medium_index,
        effective_fiber_distance=effective_fiber_distance,
        mode_field_diameter=mode_field_diameter,
        power_per_fiber=power_per_fiber,
        wavelength=wavelength,
        n_points=n_points)

    assert np.allclose(th1, th2)
    assert np.allclose(sg1, sg2)


def test_boundary():
    """Geometry should not depend on Poisson's ratio"""
    radius = 2.8466e-6
    stretch_ratio = 0.1
    poisson_ratio = 0.45

    # general formulas
    semi_minor = radius * (1 - poisson_ratio * stretch_ratio)
    semi_major = radius * (1 + stretch_ratio)

    costheta = np.cos(np.linspace(0, np.pi, 100))
    b1 = boyde2009.core.boundary(
        costheta=costheta,
        a=1,
        epsilon=(semi_major - semi_minor) / semi_minor,
        nu=0)
    b2 = boyde2009.core.boundary(
        costheta=costheta,
        a=1,
        epsilon=stretch_ratio,
        nu=poisson_ratio)

    assert np.allclose(b1, b2)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
