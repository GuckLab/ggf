import pathlib

import numpy as np

from ggf.stress.boyde2009 import core


def test_basic():
    # Originally, the reference data ere created with Boyde's Matlab
    # script using default parameters. Due to a mistake in the boundary
    # function (issue #1), the reference data were recreated with the
    # Python ggf package.
    rpath = pathlib.Path(__file__).resolve().parent / "data"
    th_ref, sigmarr_ref = np.loadtxt(str(rpath / "stress_default.dat"))
    coeff_ref = np.loadtxt(str(rpath / "coeff_default.dat"))
    th, sigmarr, coeff = core.stress(ret_legendre_decomp=True)
    assert np.allclose(th, th_ref)
    assert np.allclose(sigmarr, sigmarr_ref)
    assert np.allclose(coeff, coeff_ref)


def test_basic_barton():
    # Originally, the reference data ere created with Boyde's Matlab
    # script using default parameters. Due to a mistake in the boundary
    # function (issue #1), the reference data were recreated with the
    # Python ggf package.
    rpath = pathlib.Path(__file__).resolve().parent / "data"
    th_ref, sigmarr_ref = np.loadtxt(str(rpath / "stress_default_barton.dat"))
    coeff_ref = np.loadtxt(str(rpath / "coeff_default_barton.dat"))
    th, sigmarr, coeff = core.stress(field_approx="barton",
                                     ret_legendre_decomp=True)
    assert np.allclose(th, th_ref)
    assert np.allclose(sigmarr, sigmarr_ref)
    assert np.allclose(coeff, coeff_ref)


def test_droplet():
    rpath = pathlib.Path(__file__).resolve().parent / "data"
    coeff_ref = np.loadtxt(str(rpath / "coeff_droplet1.dat"))
    semi_major = 4.62469e-6
    semi_minor = 4.589914e-6
    # We take the average of the two to compute the stress.
    # This is wrong.
    radius = (semi_major + semi_minor) / 2
    th, sigmarr, coeff = core.stress(semi_minor=radius,
                                     object_index=1.41,
                                     medium_index=1.3465,
                                     ret_legendre_decomp=True)
    assert np.allclose(coeff, coeff_ref)


def test_barton_davis_difference():
    th1, sigmarr1 = core.stress(semi_minor=1e-6, field_approx="davis")
    th2, sigmarr2 = core.stress(semi_minor=1e-6, field_approx="barton")
    assert np.allclose(sigmarr1, sigmarr2, rtol=0, atol=3e-4)


if __name__ == "__main__":
    test_droplet()
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
