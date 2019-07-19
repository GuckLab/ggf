import numpy as np
import pytest

import ggf


def test_exact_cell():
    """Test the exact value on a LUT grid point"""
    stretch_ratio = 0.065
    semi_minor = 6.724137931034484e-06
    semi_major = semi_minor * (stretch_ratio + 1)
    f = ggf.get_ggf(model="boyde2009",
                    semi_major=semi_major,
                    semi_minor=semi_minor,
                    object_index=1.333*1.025,
                    medium_index=1.333,
                    effective_fiber_distance=175e-6,
                    mode_field_diameter=4.8e-6,
                    power_per_fiber=.65,
                    wavelength=780e-9,
                    poisson_ratio=.5,
                    n_poly=120,
                    use_lut=True)
    exact = 0.7711334992513761
    assert np.allclose(exact, f, rtol=0, atol=1e-7)


@pytest.mark.xfail
def test_exact_fus():
    """Test the exact value on a LUT grid point"""
    f = ggf.get_ggf(model="boyde2009",
                    semi_major=2.5e-6,
                    semi_minor=2.5e-6,
                    object_index=1.340*1.009,
                    medium_index=1.340,
                    effective_fiber_distance=175e-6,
                    mode_field_diameter=4.8e-6,
                    power_per_fiber=.5,
                    wavelength=780e-9,
                    poisson_ratio=.4,
                    n_poly=120,
                    use_lut=True)
    exact = 0.7878645380697753
    # upon compression this becomes: 0.78786457
    assert np.allclose(exact, f, rtol=0, atol=1e-7)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
