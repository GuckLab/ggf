import pathlib

import numpy as np

from ggf.stress.boyde2009.core import stress_legendre_decomp
from ggf.stress.boyde2009.globgeomfact import coeff2ggf


def test_basic():
    rpath = pathlib.Path(__file__).resolve().parent / "data"
    coeff1 = np.loadtxt(str(rpath / "coeff_droplet1.dat"))
    coeff2 = np.loadtxt(str(rpath / "coeff_droplet2.dat"))

    ggf1 = coeff2ggf(coeff1, poisson_ratio=.45)
    ggf2 = coeff2ggf(coeff2, poisson_ratio=.45)

    assert np.allclose(ggf1, 0.75019897326485)
    assert np.allclose(ggf2, 0.752157535649718)


def test_stress2ggf():
    sigma_0 = 1.341
    theta = np.linspace(0, np.pi, 100, endpoint=True)
    stress = sigma_0 * (np.cos(theta))**2
    poisson_ratio = 0.45
    n_poly = 10

    coeff = stress_legendre_decomp(th=theta, sigmarr=stress, n_poly=n_poly)
    ggf = coeff2ggf(coeff, poisson_ratio=poisson_ratio)
    
    # analytical solution (see Ananthakrishnan 2006, Appendix)
    nu = poisson_ratio
    fg = 1 / (2*(1+nu)) \
         * (1/3 * ((1-2*nu) + (4*nu-7)*(1+nu) / (5*nu+7)) \
            + (-4*nu + 7)*(1+nu) / (5*nu + 7) * 1) 
    
    assert np.allclose(ggf, fg*sigma_0)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()