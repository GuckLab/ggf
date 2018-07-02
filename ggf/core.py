import numpy as np

from .matlab_funcs import lscov, legendre
from .sci_funcs import legendrePlm



def legendre2ggf(coeff, poisson_ratio):
    m = 1/poisson_ratio
    Delta = lambda n: n*(n-1) + (2*n+1) * (m+1)/m 
    L_n = lambda n: -1/Delta(n) * (2*n+1) * (n+1) * (n-2+4/m)
    M_n = lambda n: 1/Delta(n) * (2*n+1) * (n**2 + 2*n - 1 + 2/m) * n / (n-1)
    # Q_n = lambda n: -1/Delta(n) * (2*n+1) * (n + 5 - 4/m)
    # S_n = lambda n: M_n(n) / n
    x = 1  # evaluate displacements at the boundary of the sphere
    theta = 0  # evaluate displacements only on the trapping axis

    # We use the notation: u_r(theta=0, R=radius) / radius = GGF / G
    # Thus, in Lur'e eq. (6.6.8), we move radius and G to the left.
    # ggf = u_r * G / radius

    ggf = 0

    for n, sn in enumerate(coeff):
        if n == 0:
            # n=0 contribution:
            ggf += (m-2) * sn / (2*(m+1))
        elif n % 2:
            if not np.allclose(sn, 0):
                msg = "Odd coeffecient n={} is non-zero: {}".format(n, sn)
                raise ValueError(msg)
        else:
            ggf += 1/8 * 2*sn / (2*n+1) \
                   * (L_n(n) * x**n + M_n(n) * x**(n-2)) \
                   * np.real_if_close(legendre(n, np.cos(theta))[0][0])

    # Note that u_theta is not considered here!    
    return ggf


def stress2legendre(stress, theta, n_poly):
    """Decompose stress into even Legendre Polynomials"""
    # Sigma = Sum_n [Coeff(n) P_n(np.cos(theta))]
    nmax = n_poly                    # number of Legendre polynomials used in fit

    # transfer data from stress plot into pair of corresponding variables
    # [Theta,Sigma]
    numpoints = theta.shape[0]
    theta = theta.reshape(-1, 1)
    sigma = stress.reshape(-1, 1)

    # Write set of linear equations for stresses in terms of Legendre functions
    legmat = np.zeros((numpoints,nmax), dtype=float)
    for ii in range(numpoints):
        for jj in np.arange(nmax)[::2]: # skip odd Legendre Polynomials since stress is an even function (symmetrical)
            legmat[ii, jj] = np.real_if_close(legendrePlm(0, jj, np.cos(theta[ii])))

    coeff = lscov(legmat, sigma)

    return coeff


def stress2ggf(stress, theta, poisson_ratio, n_poly):
    """
    """
    coeff = stress2legendre(stress=stress, theta=theta, n_poly=n_poly)

    ggf = legendre2ggf(coeff=coeff, poisson_ratio=poisson_ratio)
    
    return ggf
