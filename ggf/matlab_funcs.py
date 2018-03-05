import numpy as np
from scipy import special
from scipy import integrate


def besselh(n, z):
    """ Hankel function with k = 1
        n - real order
        z - complex argument
    """
    return special.hankel1(n,z)


def besselj(n, z):
    """ Bessel function of first kind
        n - real order
        z - complex argument
    """
    return special.jv(n,z)


def gammaln(x):
    return special.gammaln(x)


def legendre(n, x):
    """ This function behaves as the the matlab legendre function.
    
        x is treated always as a row vector.
    
        The statement legendre(2,0:0.1:0.2) returns the matrix
                x = 0   x = 0.1  x = 0.2
        m = 0  -0.5000  -0.4850     -0.4400
        m = 1   0       -0.2985  -0.5879
        m = 2   3.0000   2.9700  2.8800
    """
    x = np.real_if_close(x).flatten() # flattened
    x = np.array(x, dtype=float) # warning if not real
    result = np.zeros((x.shape[0], n+1), dtype=complex)
    for ii in range(x.shape[0]):
        # Gives us row vector
        a = special.lpmn(n, n, x[ii])[0].transpose()[-1]
        result[ii] = a
    return result.transpose()


def lscov(A, B, w=None):
    """
    A*x = B, that is, x minimizes (B - A*x)'*diag(w)*(B - A*x).
    
    w typically contains either counts or inverse variances.
    """
    # https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python
    # https://de.mathworks.com/help/matlab/ref/lscov.html
    if w is None:
        Aw = A.copy()
        Bw = B.T.copy()
    else:
        W = np.sqrt(np.diag(np.array(w).flatten()))
        Aw = np.dot(W, A)
        Bw = np.dot(B.T, W)
    
    # set rcond=1e-10 to prevent diverging odd indices in x
    # (problem specific to ggf/stress computation)
    x, residuals, rank, s = np.linalg.lstsq(Aw, Bw.T, rcond=1e-10)
    return np.array(x).flatten()


def quadl(func, a, b):
    wrapreal = lambda x: func(x).real
    wrapimag = lambda x: func(x).imag
    
    rp = integrate.quad(wrapreal,
                        np.real_if_close(a),
                        np.real_if_close(b))[0]

    ri = integrate.quad(wrapimag,
                        np.real_if_close(a),
                        np.real_if_close(b))[0]
                      
    return rp + 1j*ri

