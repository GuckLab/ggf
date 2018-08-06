import numpy as np

from .matlab_funcs import lscov, legendre
from .sci_funcs import legendrePlm


def legendre2ggf(coeff, poisson_ratio):
    """Compute the global geometric factor from Legendre coefficients

    The definition of the Legendre coefficients is given in
    :ref:`the theory section <sec_theory_ggf>`.

    Parameters
    ----------
    coeff: 1d ndarray
        Legendre coefficients as defined in :cite:`Lure1964`
    poisson_ratio: float
        Poisson's ratio of the stretched material. Set this
        to 0.5 for volume conservation.

    Returns
    -------
    ggf: float
        Global geometric factor

    Notes
    -----
    All odd Legendre coefficients are assumed to be zero, because the
    stress profile is symmetric with respect to the stretcher axis.
    """
    m = 1 / poisson_ratio

    def Delta(n): return n * (n - 1) + (2 * n + 1) * (m + 1) / m

    def L_n(n): return -1 / Delta(n) * (2 * n + 1) * (n + 1) * (n - 2 + 4 / m)

    def M_n(n): return 1 / Delta(n) * (2 * n + 1) * \
        (n**2 + 2 * n - 1 + 2 / m) * n / (n - 1)
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
            ggf += (m - 2) * sn / (2 * (m + 1))
        elif n % 2:
            if not np.allclose(sn, 0):
                msg = "Odd coeffecient n={} is non-zero: {}".format(n, sn)
                raise ValueError(msg)
        else:
            ggf += 1 / 8 * 2 * sn / (2 * n + 1) \
                * (L_n(n) * x**n + M_n(n) * x**(n - 2)) \
                * np.real_if_close(legendre(n, np.cos(theta))[0][0])

    # Note that u_theta is not considered here!
    return ggf


def stress2legendre(stress, theta, n_poly):
    """Decompose stress into even Legendre Polynomials

    The definition of the Legendre decomposition is given in
    :ref:`the theory section <sec_theory_ggf>`.

    Parameters
    ----------
    stress: 1d ndarray
        Radial stress profile (in imaging plane)
    theta: 1d ndarray
        Polar angles corresponding to `stress`
    n_poly: int
        Number of Legendre polynomials to use

    Returns
    -------
    coeff: 1d ndarray
        Legendre coefficients as defined in :cite:`Lure1964`

    Notes
    -----
    All odd Legendre coefficients are assumed to be zero, because the
    stress profile is symmetric with respect to the stretcher axis.
    Therefore, only `n_poly/2` polynomials are considered.
    """
    # Sigma = Sum_n [Coeff(n) P_n(np.cos(theta))]
    # number of Legendre polynomials used in fit
    nmax = n_poly
    # transfer data from stress plot into pair of corresponding variables
    # [Theta,Sigma]
    numpoints = theta.shape[0]
    theta = theta.reshape(-1, 1)
    sigma = stress.reshape(-1, 1)

    # Write set of linear equations for stresses in terms of Legendre functions
    legmat = np.zeros((numpoints, nmax), dtype=float)
    for ii in range(numpoints):
        # skip odd Legendre Polynomials since stress is an even function
        # (symmetrical)
        for jj in np.arange(nmax)[::2]:
            legmat[ii, jj] = np.real_if_close(
                legendrePlm(0, jj, np.cos(theta[ii])))

    coeff = lscov(legmat, sigma)

    return coeff


def stress2ggf(stress, theta, poisson_ratio, n_poly=120):
    """Compute the GGf from radial stress using Legendre decomposition

    Parameters
    ----------
    stress: 1d ndarray
        Radial stress profile (in imaging plane)
    theta: 1d ndarray
        Polar angles corresponding to `stress`
    poisson_ratio: float
        Poisson's ratio of the stretched material. Set this
        to 0.5 for volume conservation.
    n_poly: int
        Number of Legendre polynomials to use

    Returns
    -------
    ggf: float
        Global geometric factor

    Notes
    -----
    All odd Legendre coefficients are assumed to be zero, because the
    stress profile is symmetric with respect to the stretcher axis.
    Therefore, only `n_poly/2` polynomials are considered.
    """
    coeff = stress2legendre(stress=stress, theta=theta, n_poly=n_poly)

    ggf = legendre2ggf(coeff=coeff, poisson_ratio=poisson_ratio)

    return ggf
