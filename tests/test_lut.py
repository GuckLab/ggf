import numpy as np
import pytest

import ggf


@pytest.mark.xfail
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
    assert np.allclose(exact, f, rtol=0, atol=1e-10)


@pytest.mark.xfail
def test_exact_fus():
    """Test the exact value on a LUT grid point"""
    f = ggf.get_ggf(model="boyde2009",
                    semi_major=2.8e-6,
                    semi_minor=2.8e-6,
                    object_index=1.340*1.045,
                    medium_index=1.340,
                    effective_fiber_distance=170e-6,
                    mode_field_diameter=4.8e-6,
                    power_per_fiber=.65,
                    wavelength=780e-9,
                    poisson_ratio=.4,
                    n_poly=120,
                    use_lut=True)
    exact = 1.23706977246589
    assert np.allclose(exact, f, rtol=0, atol=1e-10)


@pytest.mark.xfail
def test_basic_fus_semi_major():
    f = ggf.get_ggf(model="boyde2009",
                    semi_major=2.8e-6*1.0018181818181818182,
                    semi_minor=2.8e-6,
                    object_index=1.340*1.045,
                    medium_index=1.340,
                    effective_fiber_distance=170e-6,
                    mode_field_diameter=4.8e-6,
                    power_per_fiber=.65,
                    wavelength=780e-9,
                    poisson_ratio=.4,
                    n_poly=None,
                    use_lut=True)
    exact = 1.237550828156066
    assert np.allclose(exact, f, rtol=0, atol=1e-5)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
