import numpy as np

from ggf import matlab_funcs


def test_besselj():
    # https://de.mathworks.com/help/matlab/ref/besselj.html
    z = np.arange(6) * .2
    req = matlab_funcs.besselj(1, z)
    ans = [0,
           0.0995,
           0.1960,
           0.2867,
           0.3688,
           0.4401,
           ]
    assert np.allclose(req, ans, rtol=0, atol=5e-5)


def test_besselh():
    z = np.arange(1, 6) * .2
    req = matlab_funcs.besselh(2, z)
    ansr = [0.0050,
            0.0197,
            0.0437,
            0.0758,
            0.1149,
            ]
    ansi = [-32.1571,
            -8.2983,
            -3.8928,
            -2.3586,
            -1.6507,
            ]
    assert np.allclose(req.real, ansr, rtol=0, atol=5e-5)
    assert np.allclose(req.imag, ansi, rtol=0, atol=5e-5)


def test_gammaln():
    x = np.array([1/5, 1/2, 2/3, 8/7, 3])
    req = matlab_funcs.gammaln(x)
    ans = [1.5241, 0.5724, 0.3032, -0.0667, 0.6931]
    assert np.allclose(req, ans, rtol=0, atol=5e-5)


def test_legendre():
    z = np.arange(6) * .2
    req = matlab_funcs.legendre(3, z)
    ans = [[0.000, -.2800, -.440, -.3600, .0800, 1.0000],
           [1.5000, 1.1758, 0.2750, -0.9600, -1.980, 0.0000],
           [0.0000, 2.8800, 5.0400, 5.7600, 4.3200, 0.0000],
           [-15.0000, -14.1091, -11.5481, -7.6800, -3.2400, 0.0000]
           ]
    assert np.allclose(req, ans, rtol=0, atol=5e-5)


def test_lscov():
    # https://de.mathworks.com/help/matlab/ref/lscov.html
    x1 = np.array([.2, .5, .6, .8, 1.0, 1.1]) 
    x2 = np.array([.1, .3, .4, .9, 1.1, 1.4]) 
    X = np.array([np.ones(x1.size), x1, x2]).T
    y = np.array([.17, .26, .28, .23, .27, .34])

   
    w = np.ones_like(x1)
    req = matlab_funcs.lscov(X, y, w)
    ans = [0.1203, 0.3284, -0.1312]
    assert np.allclose(req, ans, rtol=0, atol=5e-5)
    
    w2 = np.array([1, 1, 1, 1, 1, .1])
    req2 = matlab_funcs.lscov(X, y, w2)
    ans2 = [0.1046, 0.4614, -0.2621]
    assert np.allclose(req2, ans2, rtol=0, atol=5e-5)


def test_quadl():
    myfunc = lambda x: 1./(x**3-2*x-5)
    Q = matlab_funcs.quadl(myfunc, 0, 2)
    assert np.allclose(Q, -0.4605, rtol=0, atol=5e-5)


if __name__ == "__main__":
    # Run all tests
    loc = locals()
    for key in list(loc.keys()):
        if key.startswith("test_") and hasattr(loc[key], "__call__"):
            loc[key]()
