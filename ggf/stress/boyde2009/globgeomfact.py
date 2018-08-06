# flake8: noqa
"""Computation of the global geometric factor"""
import numpy as np

from ...matlab_funcs import lscov, legendre


def coeff2ggf(coeff, poisson_ratio=.45):
    """Compute the global geometric factor from stress coefficients

    The radial displacements of an elastic sphere can be expressed in
    terms of Legendre polynomials (see :cite:`Lure1964` equation 6.2.9)
    whose coefficients are computed from the Legendre decomposition of
    the radial stress.

    Notes
    -----

    - For a :math:`\\sigma_0 \cos^n(\\theta)` stress profile,
      the GGF already includes the peak stress according to:

      .. math::

          \\text{GGF} = \\sigma_0 F_\\text{G}.

    - This is a conversion of the Matlab script GGF.m to Python. The
      code solves a linear system of equations to determine all
      Legendre coefficients. The new implementation in
      :func:`ggf.legendre2ggf` uses the direct solution and
      thus should be preferred.
    """
    # Parameters used in the program

    Theta1 = 0
    # Nu =input('enter the value of nu (Poisson s ratio) \n')
    Nu = poisson_ratio
    ##

    SigmaMat = coeff
    n = SigmaMat.size - 1
    # amendement due to the correction of the paper by Ananthakrishnan
    SigmaMatamended = SigmaMat / 2
    Sigmatot = np.zeros((2 * n + 2, 1))

    for k in range(n):
        Sigmatot[k] = SigmaMatamended[k]
        #Sigmatot[k+n+1] = 0

    Mat = np.zeros((2 * n + 2, 2 * n + 2))
    r = 1  # r is set to 1 for the calculation.
    N = n

    for d in range(n + 1):
        Mat[d, d] = (r**d) * (d + 1) * (d**2 - d - 2 - 2 * Nu)
        Mat[d, d + N + 1] = (r**(d - 2)) * d * (d - 1)

    for d in range(2, n):
        # i pour le coefficient et j pour la ligne de sigma pour le coefficient
        # P(d)
        Mat[d + N, d + 1] = (r**(d + 1)) * (-1 + 2 * (d + 1) + (d + 1)**2 +
                                            2 * Nu) * (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))
        Mat[d + N, d - 1] = (r**(d - 1)) * (-1 + 2 * (d - 1) + ((d - 1)**2) +
                                            2 * Nu) * (d) * (d - 1) / (1 + 2 * (d - 1))
        Mat[d + N, d + N + 1 + 1] = (r**(d - 1)) * \
            (d) * (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))
        Mat[d + N, d + N] = (r**(d - 3)) * (d - 2) * \
            (d) * (d - 1) / (1 + 2 * (d - 1))

    d = 0
    Mat[d + N, d + 1] = (r**(d + 1)) * (-1 + 2 * (d + 1) + (d + 1)**2 + 2 * Nu) * \
        (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))
    Mat[d + N, d + N + 1 + 1] = (r**(d - 1)) * \
        (d) * (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))

    d = 1
    Mat[d + N, d + 1] = (r**(d + 1)) * (-1 + 2 * (d + 1) + ((d + 1)**2) + 2 * Nu) * \
        (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))
    Mat[d + N, d + N + 1 + 1] = (r**(d - 1)) * \
        (d) * (-(d + 1)) * (d + 2) / (1 + 2 * (d + 1))

    d = n
    Mat[d + N, d - 1] = (r**(d - 1)) * (-1 + 2 * (d - 1) + ((d - 1)**2) + 2 * Nu) * \
        (d) * (d - 1) / (1 + 2 * (d - 1))
    Mat[d + N, d + N] = (r**(d - 3)) * (d - 2) * (d) * \
        (d - 1) / (1 + 2 * (d - 1))

    d = n + 1
    Mat[d + N, d - 1] = (r**(d - 1)) * (-1 + 2 * (d - 1) + ((d - 1)**2) + 2 * Nu) * \
        (d) * (d - 1) / (1 + 2 * (d - 1))
    Mat[d + N, d + N] = (r**(d - 3)) * (d - 2) * (d) * \
        (d - 1) / (1 + 2 * (d - 1))

    Mat1 = Mat
    sol = lscov(Mat1, Sigmatot)

    PL = np.zeros(n + 1, dtype=float)
    w0 = np.zeros(n + 1, dtype=float)

    # Calculation of the Global Geometrical Factor including sigma0
    ww = 0
    for k in range(0, n):
        PLA = legendre(k, np.cos(Theta1))
        PL[k] = np.real_if_close(PLA[0][0])  # legendre polynomial with m=0!
        w0[k] = ((sol[k]) * (r**(k + 1)) * (k + 1) * (k - 2 + 4 * Nu) +
                 (sol[n + k + 1]) * r**(k - 1) * k) * PL[k]
        ww = ww + w0[k]

    GF = ww

    return GF
