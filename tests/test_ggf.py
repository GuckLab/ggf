import pathlib

import numpy as np

from ggf.globgeomfact import coeff2ggf
from scipy.stats._discrete_distns import poisson


def test_basic():
    rpath = pathlib.Path(__file__).resolve().parent / "data"
    coeff1 = np.loadtxt(str(rpath / "coeff_droplet1.dat"))
    coeff2 = np.loadtxt(str(rpath / "coeff_droplet2.dat"))

    ggf1 = coeff2ggf(coeff1, poisson_ratio=.45)
    ggf2 = coeff2ggf(coeff2, poisson_ratio=.45)

    assert np.allclose(ggf1, 0.75019897326485)
    assert np.allclose(ggf2, 0.752157535649718)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
