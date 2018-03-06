from pkg_resources import resource_filename
from functools import lru_cache
import warnings

import numpy as np

from .matlab_funcs import besselh, besselj, gammaln, lscov, quadl
from .sci_funcs import legendrePlm



def boundary(costheta, a=1, epsilon=.1, nu=0):
    """Projected boundary of a prolate spheroid
    
    Compute the boundary according to equation (4) in
    :cite:`Boyde2009` with the addition to of the
    Poisson's ratio of the object.

    .. math::

       B(\\theta) = a (1+\\epsilon) (1-\\nu \\epsilon)
       \\left[ (1+\\epsilon)^2 - \\epsilon (1+\\nu)
       (2+\\epsilon (1-\\nu)) \\cos^2 \\theta \\right]^{-1/2}
    
    Parameters
    ----------
    costheta: float or np.ndarray
        Cosine of polar coordinates :math:`\\theta`
        at which to compute the boundary.
    a: float
        Equatorial radii of prolate spheroid (semi-minor axis)
        :math:`b' = c' \\equiv a`.
    epsilon: float
        Stretch ratio; defines size of semi-major axis:
        :math:`a' = (1+\\epsilon) a`. Note that this is not
        the eccentricity of the prolate spheroid.
    nu: float
        Poisson's ratio :math:`\\nu` of the material.
    
    Returns
    -------
    B: 1d ndarray
        Radial object boundary in dependence of theta
        :math:`B(\\theta)`.

    Notes
    -----
    For :math:`\\nu=0`, the above equation becomes
    equation (4) in :cite:`Boyde2009`.
    """
    x = costheta
    B = a*(1+epsilon) \
        * (1-nu*epsilon) \
        / ((1+epsilon)**2 - epsilon*(1+nu)*(2+epsilon*(1-nu))*x**2)**.5
    return B


@lru_cache(maxsize=32)
def get_hgc():
    """Load hypergeometric coefficients from *hypergeomdata2.dat*.

    These coefficients were computed by Lars Boyde
    using Wolfram Mathematica.
    """
    hpath = resource_filename("ggf", "hypergeomdata2.dat")
    hgc = np.loadtxt(hpath)
    return hgc


def stress(object_index=1.41, medium_index=1.3465, poisson_ratio=0.45,
           radius=2.8466e-6, stretch_ratio=0.1, wavelength=780e-9,
           beam_waist=3, power_left=.6, power_right=.6, dist=100e-6,
           numpoints=100, theta_max=np.pi, field_approx="davis",
           ret_legendre_decomp=False, verbose=False):
    """Compute the stress acting on a prolate spheroid

    Parameters
    ----------
    object_index: float
        Refractive index of the spheroid
    medium_index: float
        Refractive index of the surrounding medium
    poisson_ratio: float
        Poisson's ratio of the spheroid material
    radius: float
        TODO
    stretch_ratio: float
        TODO
    wavelength: float
        Wavelenth of the gaussian beam [m]
    beam_waist: float
        Beam waist radius of the gaussian beam [wavelengths]
    power_left: float
        Laser power of the left beam [W]
    power_right: float
        Laser power of the right beam [W]
    dist: float
        Distance between beam waist and object center [m]
    numpoints: int
        Number of points to compute stresses for
    theta_max: float
        Maximum angle to compute stressed for
    field_approx: str
        TODO
    ret_legendre_decomp: bool
        If True, return coefficients of decomposition of stress
        into Legendre polynomials
    verbose: int
        Increase verbosity

    Returns
    -------
    theta: 1d ndarray
        Angles for which stresses are computed
    sigma_rr: 1d ndarray
        Radial stress corresponding to angles
    coeff: 1d ndarray
        If `ret_legendre_decomp` is True, return compositions
        of decomposition of stress into Legendre polynomials.

    Notes
    -----
    - The angles `theta` are computed on a grid that does not
      include zero and `theta_max`.
    - This implementation was first presented in :cite:`Boyde2009`.
    """
    if field_approx not in ["davis", "barton"]:
        raise ValueError("`field_approx` must be 'davis' or 'barton'")
    object_index = complex(object_index)
    medium_index = complex(medium_index)
    W0 = beam_waist * wavelength
    epsilon = stretch_ratio
    nu = poisson_ratio
    
    #ZRL = 0.5*medium_index*2*np.pi/wavelength*W0**2 # Rayleigh range [m]
    #WZ  = W0*(1+(beam_pos+d)**2/ZRL**2)**0.5   # beam waist at specified position [m]
    
    K0 = 2*np.pi/wavelength         # wave vector [m]
    Alpha = radius*K0                # size parameter  
    C = 3e8                 # speed of light [m/s]


    lmax= np.int(np.round(2+Alpha+4*(Alpha)**(1/3) + 10)) # maximum number of orders
    if lmax > 120:
        msg = 'Required number of orders for accurate expansion exceeds allowed maximum!' \
              + 'Reduce size of trapped particle!'
        raise ValueError(msg)

    if epsilon == 0:
        mmax = 3                # spherical object, no point-matching needed (mmax = 0)
    else:
        if (epsilon > 0.15):
            warnings.warn('Warning, cell stretching ratio is pretty high!') 
        mmax = 6*lmax           # spheroidal object, point-matching required (mmax has to be divisible by 3)

    EpsilonI  = medium_index**2           # permittivity in surrounding medium [1]  
    EpsilonII = object_index**2          # permittivity in within cell [1]

    MuI  = 1.000               # permeability in surrounding medium [1]
    MuII = 1.000               # permeability within cell [1]
    K1I  = 1j*K0*EpsilonI      # wave constant in Maxwell's equations (surrounding medium) [1/m]
    K1II = 1j*K0*EpsilonII     # wave constant in Maxwell's equations (within cell) [1/m]
    K2I =  1j*K0               # wave constant in Maxwell's equations (surrounding medium) [1/m]
    K2II = 1j*K0               # wave constant in Maxwell's equations (within cell) [1/m]
    KI = (-K1I*K2I)**0.5       # wave vector (surrounding medium) [1/m]
    KII = (-K1II*K2II)**0.5    # wave vector (within cell) [1/m]

    # dimensionless parameters
    k0  = 1                    # wave vector
    a   = radius * K0               # internal radius of stretched cell
    d   = dist * K0               # distance from cell centre to optical stretcher
    #ap = a*(1+stretch_ratio)             # semi-major axis (after stretching)
    #bp = a*(1-poisson_ratio*stretch_ratio)          # semi-minor axis (after stretching) 
    w0  = W0 * K0              # Gaussian width 
    k1I  = K1I/K0              # wave constant in Maxwell's equations (surrounding medium)
    k1II = K1II/K0             # wave constant in Maxwell's equations (within cell)
    k2I  = K2I/K0              # wave constant in Maxwell's equations (surrounding medium)
    k2II = K2II/K0             # wave constant in Maxwell's equations (within cell)
    kI   = KI/K0               # wave vector (surrounding medium)
    kII  = KII/K0              # wave vector (within cell)
    beta = kI                  # wave vector of Gaussian beam
    
    # other definitions
    EL = np.sqrt(power_left/(medium_index*C*W0**2))  # amplitude of electric field of left laser [kg m/(s**2 C)]
    ER = np.sqrt(power_right/(medium_index*C*W0**2))  # amplitude of electric field of right laser [kg m/(s**2 C)]
    HL = beta/k0*EL            # left laser amplitude of magnetic field
    HR = beta/k0*ER            # right laser amplitude of magnetic field
    
    zR = beta*w0**2/2          # definition of Rayleigh range
    S  = (1+1j*d/zR)**(-1)     # complex amplitude for Taylor expansion
    s  = 1/(beta*w0)           # expansion parameter for Gaussian (Barton)

    ## Functions
    # object boundary function: r(th) = a*B1(x) x= cos(th)
    B1 = lambda x: boundary(costheta=x, a=1, epsilon=epsilon, nu=nu)

    # Riccati Bessel functions and their derivatives
    psi  = lambda l, z: (np.pi/2*z)**(1/2)*besselj(l+1/2,z)                             # Riccati Bessel function (psi)
    psi1 = lambda l, z: (np.pi/(2.*z))**(1/2)*(z*besselj(l-1/2,z)-l*besselj(l+1/2,z)) # first derivative (psi')
    psi2 = lambda l, z: (np.pi/2)**(1/2)*(l+l**2-z**2)*besselj (l+1/2,z)*z**(-3/2)        # second derivative (psi'')

    # First order Taylor expansion of psi is too inaccurate for larger values of k*a*Eps.
    # Hence, to match 1-st and higher order terms in Eps, subtract the 0-th order terms (no angular dependence)
    # from the exact function psi (including angular dependence) 
    psiex  = lambda l,z,x: psi(l,z*B1(x))  - psi(l,z)    # Riccati Bessel function excluding angular dependence in 0-th order (psiex)
    psi1ex = lambda l,z,x: psi1(l,z*B1(x)) - psi1(l,z)  # first derivative of psiex
    psi2ex = lambda l,z,x: psi2(l,z*B1(x)) - psi2(l,z)  # second derivative of psi

    # defined for abbreviation 
    psixx  = lambda l,z,x: psi(l,z*B1(x)) 
    psi1xx = lambda l,z,x: psi1(l,z*B1(x)) 
    psi2xx = lambda l,z,x: psi2(l,z*B1(x)) 

    # Hankel function and its derivative
    xi   = lambda l, z: (np.pi/2*z)**(1/2)*besselh(l+1/2,z) 
    xi1  = lambda l, z: (np.pi/(2*z))**(1/2)*((l+1)*besselh(l+1/2,z)-z*besselh(l+3/2,z))
    xi2  = lambda l, z: (np.pi/2)**(1/2)/z**(3/2)*(l+l**2-z**2)*besselh(l+1/2,z) 

    # Comments: see above for psiex 
    xiex  = lambda l,z,x: xi(l,z*B1(x))  - xi(l,z) 
    xi1ex = lambda l,z,x: xi1(l,z*B1(x)) - xi1(l,z) 
    xi2ex = lambda l,z,x: xi2(l,z*B1(x)) - xi2(l,z) 
    
    xixx  = lambda l,z,x: xi(l,z*B1(x)) 
    xi1xx = lambda l,z,x: xi1(l,z*B1(x)) 
    xi2xx = lambda l,z,x: xi2(l,z*B1(x)) 
        
    
    #% Associated Legendre functions P(m)_l(x) and their derivatives 
    #% select mth component of vector 'legendre'  [P**(m)_l(x)]
    # [zeros(m,1);1;zeros(l-m,1)].'*legendre(l,x)

    
    #% legendre polynomial [P**(m)_l(x)] 
    legendrePl = lambda l, x: legendrePlm(1, l, x)
    #% legendre polynomial [P**(1)_l(x)]
    legendrePlm1 = lambda m, l, x: ((l-m+1.)*legendrePlm(m, l+1, x)-(l+1.)*x*legendrePlm(m, l, x))/(x**2 - 1)
    #% derivative d/dx[P**(m)_l(x)]
    legendrePl1 = lambda l, x: legendrePlm1(1,l,x)

    # defined to avoid division by zero (which can occur for x=1 in legendrePl1...
    legendrePlmex1  = lambda m,l,x: -((l-m+1)*legendrePlm(m,l+1,x)-(l+1)*x*legendrePlm(m,l,x))   
    legendrePlex1   = lambda l, x: legendrePlmex1(1,l,x)                                            

    # Hypergeometric and Gamma functions 
    hypergeomcoeff = get_hgc()


    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    # Gaussian beam (incident fields) - in Cartesian basis (either Davis first order or Barton fifth order fields)
    if field_approx == "davis": # electric and magnetic fields according to Davis (first order)
        # left
        eExiL = lambda r,th,phi: EL *(1+1j*(r*np.cos(th)+d)/zR)**(-1)*np.exp (-r**2*np.sin(th)**2/(w0**2*(1+1j*(r*np.cos(th)+d)/zR)))*np.exp(1j*beta*(r*np.cos(th)+d))
        eEyiL = lambda r,th,phi: 0                                                                                                                          
        eEziL = lambda r,th,phi: -1j*(1+1j*(r*np.cos(th)+d)/zR)**(-1) *r*np.sin(th)*np.cos(phi)/zR * eExiL(r,th,phi)                                           
        eHxiL = lambda r,th,phi: 0                                                                                                                          
        eHyiL = lambda r,th,phi: HL *(1+1j*(r*np.cos(th)+d)/zR)**(-1)*np.exp(-r**2*np.sin(th)**2/(w0**2*(1+1j*(r*np.cos(th)+d)/zR))) *np.exp(1j*beta*(r*np.cos(th)+d))  
        eHziL = lambda r,th,phi: -1j*(1+1j*(r*np.cos(th)+d)/zR)**(-1) *r*np.sin(th)*np.sin(phi)/zR * eHyiL(r,th,phi)                                           
    
        # right
        eExiR = lambda r,th,phi: ER *(1-1j*(r*np.cos(th)-d)/zR)**(-1)*np.exp(-r**2*np.sin(th)**2/(w0**2*(1-1j*(r*np.cos(th)-d)/zR))) *np.exp(-1j*beta*(r*np.cos(th)-d)) 
        eEyiR = lambda r,th,phi: 0                                                                                                                          
        eEziR = lambda r,th,phi: +1j*(1-1j*(r*np.cos(th)-d)/zR)**(-1) *r*np.sin(th)*np.cos(phi)/zR * eExiR(r,th,phi)                                           
        eHxiR = lambda r,th,phi: 0                                                                                                                          
        eHyiR = lambda r,th,phi: -HR *(1-1j*(r*np.cos(th)-d)/zR)**(-1)*np.exp(-r**2*np.sin(th)**2/(w0**2*(1-1j*(r*np.cos(th)-d)/zR))) *np.exp(-1j*beta*(r*np.cos(th)-d))
        eHziR = lambda r,th,phi: +1j*(1-1j*(r*np.cos(th)-d)/zR)**(-1) *r*np.sin(th)*np.sin(phi)/zR * eHyiR(r,th,phi)    

    else: # electric and magnetic fields according to Barton (fifth order)
        # Note: in Barton propagation is: np.exp(i (omega*t-k*z)) and not np.exp(i (k*z-omega*t)).
        # Hence, take complex conjugate of equations or make the changes z -> -z, Ez -> -Ez, Hz -> -Hz while x and y are not affected.
        Xi   = lambda r,th,phi: r*np.sin(th)*np.cos (phi)/w0
        Eta  = lambda r,th,phi: r*np.sin(th)*np.sin(phi)/w0
        ZetaL = lambda r,th:     (r*np.cos(th)+d)/(beta*w0**2)
        Rho  = lambda r,th,phi: np.sqrt(Xi (r,th,phi)**2 + Eta (r,th,phi)**2)
        QL    = lambda r,th:     1./(1j+2*ZetaL(r,th))
        Psi0L = lambda r,th,phi: 1j*QL (r,th)*np.exp(-1j*(Rho(r,th,phi))**2*QL(r,th))

        eExiL = lambda r,th,phi: np.conj( EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               ( 1+s**2*(-Rho(r,th,phi)**2*QL (r,th)**2 + 1j*Rho (r,th,phi)**4*QL (r,th)**3 - 2*QL (r,th)**2*Xi (r,th,phi)**2) +  \
                                 s**4*(  2*Rho (r,th,phi)**4*QL (r,th)**4 -3*1j*Rho (r,th,phi)**6*QL (r,th)**5 -0.5*Rho (r,th,phi)**8*QL (r,th)**6 +  \
                                (8*Rho (r,th,phi)**2*QL (r,th)**4 -2*1j*Rho (r,th,phi)**4*QL (r,th)**5)*Xi (r,th,phi)**2  ) ) )  
        eEyiL = lambda r,th,phi: np.conj( EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               ( s**2*(-2*QL (r,th)**2*Xi (r,th,phi)*Eta(r,th,phi)) +  \
                                 s**4*((8*Rho (r,th,phi)**2*QL (r,th)**4 -2*1j*Rho (r,th,phi)**4*QL (r,th)**5)*Xi (r,th,phi)*Eta(r,th,phi)) ) )
        eEziL = lambda r,th,phi: np.conj( EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               ( s*(-2*QL (r,th)*Xi(r,th,phi)) + s**3*(6*Rho (r,th,phi)**2*QL (r,th)**3 - 2*1j*Rho (r,th,phi)**4*QL (r,th)**4)*Xi(r,th,phi) +  \
                                 s**5*(-20*Rho (r,th,phi)**4*QL (r,th)**5 + 10*1j*Rho (r,th,phi)**6*QL (r,th)**6 + Rho (r,th,phi)**8*QL (r,th)**7)*Xi(r,th,phi) ) )
        eHxiL = lambda r,th,phi: np.conj( +np.sqrt(EpsilonI)*EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               ( s**2*(-2*QL (r,th)**2*Xi (r,th,phi)*Eta(r,th,phi)) + \
                                 s**4*((8*Rho (r,th,phi)**2*QL (r,th)**4 -2*1j*Rho (r,th,phi)**4*QL (r,th)**5)*Xi (r,th,phi)*Eta(r,th,phi)) ) )
        eHyiL = lambda r,th,phi: np.conj( +np.sqrt(EpsilonI)*EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               ( 1+s**2*(-Rho(r,th,phi)**2*QL (r,th)**2 + 1j*Rho (r,th,phi)**4*QL (r,th)**3 - 2*QL (r,th)**2*Eta (r,th,phi)**2) +  \
                                 s**4*(  2*Rho (r,th,phi)**4*QL (r,th)**4 -3*1j*Rho (r,th,phi)**6*QL (r,th)**5 -0.5*Rho (r,th,phi)**8*QL (r,th)**6 +  \
                                (8*Rho (r,th,phi)**2*QL (r,th)**4 -2*1j*Rho (r,th,phi)**4*QL (r,th)**5)*Eta (r,th,phi)**2  ) ) )  
        eHziL = lambda r,th,phi: np.conj( np.sqrt(EpsilonI)*EL *Psi0L (r,th,phi)*np.exp(-1j*ZetaL (r,th)/s**2)* \
                               (s*(-2*QL (r,th)*Eta(r,th,phi)) + s**3*(6*Rho (r,th,phi)**2*QL (r,th)**3 - 2*1j*Rho (r,th,phi)**4*QL (r,th)**4)*Eta(r,th,phi) +  \
                                s**5*(-20*Rho (r,th,phi)**4*QL (r,th)**5 + 10*1j*Rho (r,th,phi)**6*QL (r,th)**6 + Rho (r,th,phi)**8*QL (r,th)**7)*Eta(r,th,phi) ) )

        # right
        # Take left fiber fields and make coordinate changes (x,y,z,d) -> (x,-y,-z,d) and amplitude changes (Ex,Ey,Ez) -> (Ex,-Ey,-Ez).
        ZetaR = lambda r,th:     -(r*np.cos(th)-d)/(beta*w0**2)
        QR    = lambda r,th:     1./(1j+2*ZetaR(r,th))
        Psi0R = lambda r,th,phi: 1j*QR (r,th)*np.exp(-1j*(Rho(r,th,phi))**2*QR(r,th))

        eExiR = lambda r,th,phi: np.conj( ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               ( 1+s**2*(-Rho(r,th,phi)**2*QR (r,th)**2 + 1j*Rho (r,th,phi)**4*QR (r,th)**3 - 2*QR (r,th)**2*Xi (r,th,phi)**2) +  \
                                 s**4*(  2*Rho (r,th,phi)**4*QR (r,th)**4 -3*1j*Rho (r,th,phi)**6*QR (r,th)**5 -0.5*Rho (r,th,phi)**8*QR (r,th)**6 +  \
                                (8*Rho (r,th,phi)**2*QR (r,th)**4 -2*1j*Rho (r,th,phi)**4*QR (r,th)**5)*Xi (r,th,phi)**2  ) ) )  
        eEyiR = lambda r,th,phi: np.conj( ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               ( s**2*(-2*QR (r,th)**2*Xi (r,th,phi)*Eta(r,th,phi)) +  \
                                 s**4*((8*Rho (r,th,phi)**2*QR (r,th)**4 -2*1j*Rho (r,th,phi)**4*QR (r,th)**5)*Xi (r,th,phi)*Eta(r,th,phi)) ) )
        eEziR = lambda r,th,phi: - np.conj( ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               ( s*(-2*QR (r,th)*Xi(r,th,phi)) + s**3*(6*Rho (r,th,phi)**2*QR (r,th)**3 - 2*1j*Rho (r,th,phi)**4*QR (r,th)**4)*Xi(r,th,phi) +  \
                                 s**5*(-20*Rho (r,th,phi)**4*QR (r,th)**5 + 10*1j*Rho (r,th,phi)**6*QR (r,th)**6 + Rho (r,th,phi)**8*QR (r,th)**7)*Xi(r,th,phi) ) )
        eHxiR = lambda r,th,phi: - np.conj( +np.sqrt(EpsilonI)*ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               ( s**2*(-2*QR (r,th)**2*Xi (r,th,phi)*Eta(r,th,phi)) + \
                                 s**4*((8*Rho (r,th,phi)**2*QR (r,th)**4 -2*1j*Rho (r,th,phi)**4*QR (r,th)**5)*Xi (r,th,phi)*Eta(r,th,phi)) ) )
        eHyiR = lambda r,th,phi: - np.conj( +np.sqrt(EpsilonI)*ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               ( 1+s**2*(-Rho(r,th,phi)**2*QR (r,th)**2 + 1j*Rho (r,th,phi)**4*QR (r,th)**3 - 2*QR (r,th)**2*Eta (r,th,phi)**2) +  \
                                 s**4*(  2*Rho (r,th,phi)**4*QR (r,th)**4 -3*1j*Rho (r,th,phi)**6*QR (r,th)**5 -0.5*Rho (r,th,phi)**8*QR (r,th)**6 +  \
                                (8*Rho (r,th,phi)**2*QR (r,th)**4 -2*1j*Rho (r,th,phi)**4*QR (r,th)**5)*Eta (r,th,phi)**2  ) ) )  
        eHziR = lambda r,th,phi: np.conj( np.sqrt(EpsilonI)*ER *Psi0R (r,th,phi)*np.exp(-1j*ZetaR (r,th)/s**2)* \
                               (s*(-2*QR (r,th)*Eta(r,th,phi)) + s**3*(6*Rho (r,th,phi)**2*QR (r,th)**3 - 2*1j*Rho (r,th,phi)**4*QR (r,th)**4)*Eta(r,th,phi) +  \
                                s**5*(-20*Rho (r,th,phi)**4*QR (r,th)**5 + 10*1j*Rho (r,th,phi)**6*QR (r,th)**6 + Rho (r,th,phi)**8*QR (r,th)**7)*Eta(r,th,phi) ) )

    # Gaussian beam (incident fields) - in spherical polar coordinates basis
    # left
    eEriL   = lambda r, th, phi: np.sin(th)*np.cos(phi)*eExiL(r,th,phi) + np.sin(th)*np.sin(phi)*eEyiL(r,th,phi) + np.cos(th)*eEziL(r,th,phi)     
    eEthiL  = lambda r, th, phi: np.cos(th)*np.cos(phi)*eExiL(r,th,phi) + np.cos(th)*np.sin(phi)*eEyiL(r,th,phi) - np.sin(th)*eEziL(r,th,phi)   
    eEphiiL = lambda r, th, phi:        - np.sin(phi)*eExiL(r,th,phi)          + np.cos(phi)*eEyiL(r,th,phi)                                                          
    eHriL   = lambda r, th, phi: np.sin(th)*np.cos(phi)*eHxiL(r,th,phi) + np.sin(th)*np.sin(phi)*eHyiL(r,th,phi) + np.cos(th)*eHziL(r,th,phi)     
    eHthiL  = lambda r, th, phi: np.cos(th)*np.cos(phi)*eHxiL(r,th,phi) + np.cos(th)*np.sin(phi)*eHyiL(r,th,phi) - np.sin(th)*eHziL(r,th,phi)   
    eHphiiL = lambda r, th, phi:        - np.sin(phi)*eHxiL(r,th,phi)          + np.cos(phi)*eHyiL(r,th,phi)                                                          
    # right
    eEriR   = lambda r, th, phi: np.sin(th)*np.cos(phi)*eExiR(r,th,phi) + np.sin(th)*np.sin(phi)*eEyiR(r,th,phi) + np.cos(th)*eEziR(r,th,phi)     
    eEthiR  = lambda r, th, phi: np.cos(th)*np.cos(phi)*eExiR(r,th,phi) + np.cos(th)*np.sin(phi)*eEyiR(r,th,phi) - np.sin(th)*eEziR(r,th,phi)   
    eEphiiR = lambda r, th, phi:        - np.sin(phi)*eExiR(r,th,phi)          + np.cos(phi)*eEyiR(r,th,phi)                                                          
    eHriR   = lambda r, th, phi: np.sin(th)*np.cos(phi)*eHxiR(r,th,phi) + np.sin(th)*np.sin(phi)*eHyiR(r,th,phi) + np.cos(th)*eHziR(r,th,phi)     
    eHthiR  = lambda r, th, phi: np.cos(th)*np.cos(phi)*eHxiR(r,th,phi) + np.cos(th)*np.sin(phi)*eHyiR(r,th,phi) - np.sin(th)*eHziR(r,th,phi)   
    eHphiiR = lambda r, th, phi:        - np.sin(phi)*eHxiR(r,th,phi)          + np.cos(phi)*eHyiR(r,th,phi)                                                          


    eBL = np.zeros(lmax, dtype=complex)
    mBL = np.zeros(lmax, dtype=complex)
    eBR = np.zeros(lmax, dtype=complex)
#    eBR_test = np.zeros(lmax, dtype=complex)
    mBR = np.zeros(lmax, dtype=complex)

    if field_approx == "davis": # Davis 
        tau = np.zeros(lmax, dtype=complex)
        for ii in range(lmax):
            l = ii + 1
            tau[ii] = 0                             # sum over m and n
            niimax = int(np.floor((l-1)/3))
            for n in range(niimax+1):
                miimax = int(np.floor((l-1)/2 - 3/2*n))
                for m in range(miimax+1):
                    if l < (2*m + 3*n + 2):
                        Delta = 0
                    else:
                        Delta = 1
                    tau[ii]=tau[ii]+ np.exp(gammaln(l/2-m-n)+gammaln(m+n+2)-gammaln(l-2*m-3*n)-gammaln(l/2+2)-gammaln(m+1)-gammaln(n+1))*hypergeomcoeff[l-1,m+n] \
                                    *(-1)**(m+1)/(S**m *(beta*zR)**n) *(S/(1j*beta*w0))**(2*m+2*n)*(1-Delta*2*S/(beta*zR)*(l-1-2*m-3*n))
#                    print(tau)
                                   
        # calculate expansion coefficients of order l for electric and magnetic Debye potentials 
        emB = lambda l: S*np.exp(1j*beta*d)*(1j*beta/(2*kI))**(l-1)*(l+1/2)**2/(l+1)*np.exp(gammaln(2*l)-gammaln(l+1))*tau[l-1]
        for ii in range(lmax):
            l = ii + 1
            # left
            eBL[ii] = EL*emB(l)            
            mBL[ii] = HL*emB(l)            
            # right
            eBR[ii] = ER*emB(l) # should include factor of (-1)**(l-1) for symmetry reasons, this is taken into account only after least square fitting (see below)
            mBR[ii] = HR*emB(l) # should include factor of (-1)**l for symmetry reasons, this is taken into account only after least square fitting (see 
    else: # Barton 
        
        for ii in range(lmax):
            l = ii + 1
            eBL[ii] = np.sqrt(2)*quadl(lambda th: eEriL(a,th,np.pi/4)*legendrePl(l,np.cos(th))*np.sin(th),0,np.pi)*(2*l+1)/(2*l*(l+1))*a**2/psi(l,kI*a)*kI**2/(l*(l+1))
            mBL[ii] = np.sqrt(2)*quadl(lambda th: eHriL(a,th,np.pi/4)*legendrePl(l,np.cos(th))*np.sin(th),0,np.pi)*(2*l+1)/(2*l*(l+1))*a**2/psi(l,kI*a)*kI**2/(l*(l+1))
            eBR[ii] = np.sqrt(2)*quadl(lambda th: eEriR(a,th,np.pi/4)*legendrePl(l,np.cos(th))*np.sin(th),0,np.pi)*(2*l+1)/(2*l*(l+1))*a**2/psi(l,kI*a)*kI**2/(l*(l+1))
            mBR[ii] = np.sqrt(2)*quadl(lambda th: eHriR(a,th,np.pi/4)*legendrePl(l,np.cos(th))*np.sin(th),0,np.pi)*(2*l+1)/(2*l*(l+1))*a**2/psi(l,kI*a)*kI**2/(l*(l+1))

            eBR[ii] = eBR[ii] * (-1)**(l-1)  # make symetrical with left expansion coefficients
            mBR[ii] = mBR[ii] *(-1)**l      # make symetrical with left expansion coefficients

    # coefficients for internal fields (eCl, mCl) and scattered fields (eDl, mDl) 
    eCL = np.zeros(lmax, dtype=complex)
    mCL = np.zeros(lmax, dtype=complex)
    eCR = np.zeros(lmax, dtype=complex)
    mCR = np.zeros(lmax, dtype=complex)

    eDL = np.zeros(lmax, dtype=complex)
    mDL = np.zeros(lmax, dtype=complex)
    eDR = np.zeros(lmax, dtype=complex)
    mDR = np.zeros(lmax, dtype=complex)
    for ii in range(lmax):
        l = ii + 1
        # internal (left and right)   
        eCL[ii] = k1I/kI*(kII)**2*(xi(l,kI*a)*psi1(l,kI*a)-xi1(l,kI*a)*psi(l,kI*a))/(k1I*kII*xi(l,kI*a)*psi1(l,kII*a)-k1II*kI*xi1(l,kI*a)*psi(l,kII*a))*eBL[ii]
        mCL[ii] = k2I/kI*(kII)**2*(xi(l,kI*a)*psi1(l,kI*a)-xi1(l,kI*a)*psi(l,kI*a))/(k2I*kII*xi(l,kI*a)*psi1(l,kII*a)-k2II*kI*xi1(l,kI*a)*psi(l,kII*a))*mBL[ii]
        eCR[ii] = k1I/kI*(kII)**2*(xi(l,kI*a)*psi1(l,kI*a)-xi1(l,kI*a)*psi(l,kI*a))/(k1I*kII*xi(l,kI*a)*psi1(l,kII*a)-k1II*kI*xi1(l,kI*a)*psi(l,kII*a))*eBR[ii]
        mCR[ii] = k2I/kI*(kII)**2*(xi(l,kI*a)*psi1(l,kI*a)-xi1(l,kI*a)*psi(l,kI*a))/(k2I*kII*xi(l,kI*a)*psi1(l,kII*a)-k2II*kI*xi1(l,kI*a)*psi(l,kII*a))*mBR[ii]
        # scattered (left and right)
        eDL[ii] = (k1I*kII*psi(l,kI*a)*psi1(l,kII*a)-k1II*kI*psi1(l,kI*a)*psi(l,kII*a))/(k1II*kI*xi1(l,kI*a)*psi(l,kII*a)-k1I*kII*xi(l,kI*a)*psi1(l,kII*a))*eBL[ii]
        mDL[ii] = (k2I*kII*psi(l,kI*a)*psi1(l,kII*a)-k2II*kI*psi1(l,kI*a)*psi(l,kII*a))/(k2II*kI*xi1(l,kI*a)*psi(l,kII*a)-k2I*kII*xi(l,kI*a)*psi1(l,kII*a))*mBL[ii]
        eDR[ii] = (k1I*kII*psi(l,kI*a)*psi1(l,kII*a)-k1II*kI*psi1(l,kI*a)*psi(l,kII*a))/(k1II*kI*xi1(l,kI*a)*psi(l,kII*a)-k1I*kII*xi(l,kI*a)*psi1(l,kII*a))*eBR[ii]
        mDR[ii] = (k2I*kII*psi(l,kI*a)*psi1(l,kII*a)-k2II*kI*psi1(l,kI*a)*psi(l,kII*a))/(k2II*kI*xi1(l,kI*a)*psi(l,kII*a)-k2I*kII*xi(l,kI*a)*psi1(l,kII*a))*mBR[ii]

    ## First Order Expansion Coefficients 
    # coefficients for internal fields (eCcl, mCcl) and scattered fields (eDdl, mDdl) 
    eLambda1L = {}
    mLambda1L = {}
    eLambda2L = {}
    mLambda2L = {}    
    eLambda3L = {}
    mLambda3L = {}    
    eLambda1R = {}
    mLambda1R = {}
    eLambda2R = {}
    mLambda2R = {}    
    eLambda3R = {}
    mLambda3R = {}
    
    for jj in range(lmax):
        l = jj + 1
        # left       
        eLambda1L[l] = lambda x, l=l, jj=jj: (eBL[jj]/kI*psi1ex(l,kI*a,x)-eCL[jj]/kII*psi1ex(l,kII*a,x)+eDL[jj]/kI*xi1ex(l,kI*a,x))                   # electric parameter1 left
        mLambda1L[l] = lambda x, l=l, jj=jj: (mBL[jj]/kI*psi1ex(l,kI*a,x)-mCL[jj]/kII*psi1ex(l,kII*a,x)+mDL[jj]/kI*xi1ex(l,kI*a,x))                   # magnetic parameter1 left
        eLambda2L[l] = lambda x, l=l, jj=jj: (k1I/kI**2*eBL[jj]*psiex(l,kI*a,x)-k1II/kII**2*eCL[jj]*psiex(l,kII*a,x)+k1I/kI**2*eDL[jj]*xiex(l,kI*a,x))   # electric parameter2 left
        mLambda2L[l] = lambda x, l=l, jj=jj: (k2I/kI**2*mBL[jj]*psiex(l,kI*a,x)-k2II/kII**2*mCL[jj]*psiex(l,kII*a,x)+k2I/kI**2*mDL[jj]*xiex(l,kI*a,x))   # magnetic parameter2 left
        eLambda3L[l] = lambda x, l=l, jj=jj: (eBL[jj]*(psiex(l,kI*a,x)+psi2ex(l,kI*a,x))*k1I-eCL[jj]*(psiex(l,kII*a,x)+psi2ex(l,kII*a,x))*k1II \
                                             + eDL[jj]*(xiex(l,kI*a,x)+xi2ex(l,kI*a,x))*k1I)                                                       # electric parameter3 left    
        mLambda3L[l] = lambda x, l=l, jj=jj: (mBL[jj]*(psiex(l,kI*a,x)+psi2ex(l,kI*a,x))*MuI-mCL[jj]*(psiex(l,kII*a,x)+psi2ex(l,kII*a,x))*MuII \
                                             + mDL[jj]*(xiex(l,kI*a,x)+xi2ex(l,kI*a,x))*MuI)                                                       # magnetic parameter3 left                         
        # right
        eLambda1R[l] = lambda x, l=l, jj=jj: (eBR[jj]/kI*psi1ex(l,kI*a,x)-eCR[jj]/kII*psi1ex(l,kII*a,x)+eDR[jj]/kI*xi1ex(l,kI*a,x))                   # electric parameter1 right
        mLambda1R[l] = lambda x, l=l, jj=jj: (mBR[jj]/kI*psi1ex(l,kI*a,x)-mCR[jj]/kII*psi1ex(l,kII*a,x)+mDR[jj]/kI*xi1ex(l,kI*a,x))                   # magnetic parameter1 right
        eLambda2R[l] = lambda x, l=l, jj=jj: (k1I/kI**2*eBR[jj]*psiex(l,kI*a,x)-k1II/kII**2*eCR[jj]*psiex(l,kII*a,x)+k1I/kI**2*eDR[jj]*xiex(l,kI*a,x))   # electric parameter2 right
        mLambda2R[l] = lambda x, l=l, jj=jj: (k2I/kI**2*mBR[jj]*psiex(l,kI*a,x)-k2II/kII**2*mCR[jj]*psiex(l,kII*a,x)+k2I/kI**2*mDR[jj]*xiex(l,kI*a,x))   # magnetic parameter2 right
        eLambda3R[l] = lambda x, l=l, jj=jj: (eBR[jj]*(psiex(l,kI*a,x)+psi2ex(l,kI*a,x))*k1I-eCR[jj]*(psiex(l,kII*a,x)+psi2ex(l,kII*a,x))*k1II \
                                             + eDR[jj]*(xiex(l,kI*a,x)+xi2ex(l,kI*a,x))*k1I)                                                       # electric parameter3 right                              
        mLambda3R[l] = lambda x, l=l, jj=jj: (mBR[jj]*(psiex(l,kI*a,x)+psi2ex(l,kI*a,x))*MuI-mCR[jj]*(psiex(l,kII*a,x)+psi2ex(l,kII*a,x))*MuII \
                                             + mDR[jj]*(xiex(l,kI*a,x)+xi2ex(l,kI*a,x))*MuI)                                                          # magnetic parameter3 right                                 

    # define points for least square fitting [similar to operating on equations with int_ {-1}**(+1} dx legendrePl(m,x)...]
    x = {}
    for m in range(1, int(mmax/3)+1):
        x[int(m)] = (-1 + 0.1*(m-1)/(mmax/3))
        x[int(m+mmax/3)] = (-0.9 + 1.8*m/(mmax/3))
        x[int(m+2*mmax/3)] = (+0.9 + 0.1*m/(mmax/3))

    
    efun1aL = np.zeros((mmax, lmax), dtype=complex)
    efun1bL = np.zeros((mmax, lmax), dtype=complex)
    efun2aL = np.zeros((mmax, lmax), dtype=complex)
    efun2bL = np.zeros((mmax, lmax), dtype=complex)
    efun3L = np.zeros((mmax, lmax), dtype=complex)

    mfun1aL = np.zeros((mmax, lmax), dtype=complex)
    mfun1bL = np.zeros((mmax, lmax), dtype=complex)
    mfun2aL = np.zeros((mmax, lmax), dtype=complex)
    mfun2bL = np.zeros((mmax, lmax), dtype=complex)
    mfun3L = np.zeros((mmax, lmax), dtype=complex)

    efun1aR = np.zeros((mmax, lmax), dtype=complex)
    efun1bR = np.zeros((mmax, lmax), dtype=complex)
    efun2aR = np.zeros((mmax, lmax), dtype=complex)
    efun2bR = np.zeros((mmax, lmax), dtype=complex)
    efun3R = np.zeros((mmax, lmax), dtype=complex)

    mfun1aR = np.zeros((mmax, lmax), dtype=complex)
    mfun1bR = np.zeros((mmax, lmax), dtype=complex)
    mfun2aR = np.zeros((mmax, lmax), dtype=complex)
    mfun2bR = np.zeros((mmax, lmax), dtype=complex)
    mfun3R = np.zeros((mmax, lmax), dtype=complex)

    for ii in range(mmax):
        m = ii + 1
        # define points for least square fitting [similar to operating on equations with int_ {-1}**(+1} dx legendrePl(m,x)...]
        # the points x[m] lie in the range [-1,1] or similarly th(m) is in the range [0,pi]
        # x[m] = (-1 + 2*(m-1)/(mmax-1))
        for jj in range(lmax):
            l = jj + 1
            # left
            efun1aL[ii,jj] = eLambda1L[l](x[m])*legendrePl(l,x[m])
            efun1bL[ii,jj] = eLambda1L[l](x[m])*legendrePlex1(l,x[m])
            mfun1aL[ii,jj] = mLambda1L[l](x[m])*legendrePl(l,x[m])
            mfun1bL[ii,jj] = mLambda1L[l](x[m])*legendrePlex1(l,x[m])
            mfun2aL[ii,jj] = mLambda2L[l](x[m])*legendrePl(l,x[m])
            mfun2bL[ii,jj] = mLambda2L[l](x[m])*legendrePlex1(l,x[m])
            efun2aL[ii,jj] = eLambda2L[l](x[m])*legendrePl(l,x[m])
            efun2bL[ii,jj] = eLambda2L[l](x[m])*legendrePlex1(l,x[m])
            efun3L[ii,jj]  = eLambda3L[l](x[m])*legendrePl(l,x[m])
            mfun3L[ii,jj]  = mLambda3L[l](x[m])*legendrePl(l,x[m])
    
            # right
            efun1aR[ii,jj] = eLambda1R[l](x[m])*legendrePl(l,x[m])
            efun1bR[ii,jj] = eLambda1R[l](x[m])*legendrePlex1(l,x[m])
            mfun1aR[ii,jj] = mLambda1R[l](x[m])*legendrePl(l,x[m])
            mfun1bR[ii,jj] = mLambda1R[l](x[m])*legendrePlex1(l,x[m])
            mfun2aR[ii,jj] = mLambda2R[l](x[m])*legendrePl(l,x[m])
            mfun2bR[ii,jj] = mLambda2R[l](x[m])*legendrePlex1(l,x[m])
            efun2aR[ii,jj] = eLambda2R[l](x[m])*legendrePl(l,x[m])
            efun2bR[ii,jj] = eLambda2R[l](x[m])*legendrePlex1(l,x[m])
            efun3R[ii,jj]  = eLambda3R[l](x[m])*legendrePl(l,x[m])
            mfun3R[ii,jj]  = mLambda3R[l](x[m])*legendrePl(l,x[m])

    # first order BC can be written in form: M11*eCc + M12*eDd + M13*mCc + M14*mDd = N1
    #                                        ...
    #                                        M61*eCc + M62*eDd + M63*mCc + M64*mDd = N6
    # M11...M66 are matrices including the pre-factors for the individual terms in the sums of eCc[jj], etc
    # eCc...mDd are vectors including all the expansion coefficients (eCc = [eCc(1),...,eCc(lmax)])
    # N1...N6 include the 0-th order terms and angular-dependent Legendre and Bessel functions
 
    N1L = np.zeros(mmax, dtype=complex)
    N2L = np.zeros(mmax, dtype=complex)
    N3L = np.zeros(mmax, dtype=complex)
    N4L = np.zeros(mmax, dtype=complex)
    N5L = np.zeros(mmax, dtype=complex)
    N6L = np.zeros(mmax, dtype=complex)

    N1R = np.zeros(mmax, dtype=complex)
    N2R = np.zeros(mmax, dtype=complex)
    N3R = np.zeros(mmax, dtype=complex)
    N4R = np.zeros(mmax, dtype=complex)
    N5R = np.zeros(mmax, dtype=complex)
    N6R = np.zeros(mmax, dtype=complex)
    
    for ii in range(mmax):
        # left
        N1L[ii] = np.sum(efun1aL[ii,:]) - np.sum(mfun2bL[ii,:]) 
        N2L[ii] = np.sum(efun1bL[ii,:]) - np.sum(mfun2aL[ii,:]) 
        N3L[ii] = np.sum(mfun1aL[ii,:]) - np.sum(efun2bL[ii,:]) 
        N4L[ii] = np.sum(mfun1bL[ii,:]) - np.sum(efun2aL[ii,:]) 
        N5L[ii] = np.sum(efun3L[ii,:])
        N6L[ii] = np.sum(mfun3L[ii,:])
        # right
        N1R[ii] = np.sum(efun1aR[ii,:]) - np.sum(mfun2bR[ii,:]) 
        N2R[ii] = np.sum(efun1bR[ii,:]) - np.sum(mfun2aR[ii,:]) 
        N3R[ii] = np.sum(mfun1aR[ii,:]) - np.sum(efun2bR[ii,:]) 
        N4R[ii] = np.sum(mfun1bR[ii,:]) - np.sum(efun2aR[ii,:]) 
        N5R[ii] = np.sum(efun3R[ii,:])
        N6R[ii] = np.sum(mfun3R[ii,:])

    ##
    M11 = np.zeros((mmax, lmax), dtype=complex)
    M12 = np.zeros((mmax, lmax), dtype=complex)
    M13 = np.zeros((mmax, lmax), dtype=complex)
    M14 = np.zeros((mmax, lmax), dtype=complex)

    M21 = np.zeros((mmax, lmax), dtype=complex)
    M22 = np.zeros((mmax, lmax), dtype=complex)
    M23 = np.zeros((mmax, lmax), dtype=complex)
    M24 = np.zeros((mmax, lmax), dtype=complex)
    
    M31 = np.zeros((mmax, lmax), dtype=complex)
    M32 = np.zeros((mmax, lmax), dtype=complex)
    M33 = np.zeros((mmax, lmax), dtype=complex)
    M34 = np.zeros((mmax, lmax), dtype=complex)

    M41 = np.zeros((mmax, lmax), dtype=complex)
    M42 = np.zeros((mmax, lmax), dtype=complex)
    M43 = np.zeros((mmax, lmax), dtype=complex)
    M44 = np.zeros((mmax, lmax), dtype=complex)

    M51 = np.zeros((mmax, lmax), dtype=complex)
    M52 = np.zeros((mmax, lmax), dtype=complex)
    M53 = np.zeros((mmax, lmax), dtype=complex)
    M54 = np.zeros((mmax, lmax), dtype=complex)

    M61 = np.zeros((mmax, lmax), dtype=complex)
    M62 = np.zeros((mmax, lmax), dtype=complex)
    M63 = np.zeros((mmax, lmax), dtype=complex)
    M64 = np.zeros((mmax, lmax), dtype=complex)
    
    for ii in range(mmax):
        m = ii + 1
        for jj in range(lmax):
            l = jj + 1
    
            M11[ii,jj] = +1/kII  *psi1xx (l,kII*a,x[m])*legendrePl(l,x[m]) 
            M12[ii,jj] = -1/kI   *xi1xx (l,kI*a,x[m])  *legendrePl(l,x[m]) 
            M13[ii,jj] = +1/k1II *psixx (l,kII*a,x[m]) *legendrePlex1(l,x[m]) 
            M14[ii,jj] = -1/k1I  *xixx (l,kI*a,x[m])   *legendrePlex1(l,x[m]) 
    
            M21[ii,jj] = +1/kII  *psi1xx (l,kII*a,x[m])*legendrePlex1(l,x[m]) 
            M22[ii,jj] = -1/kI   *xi1xx (l,kI*a,x[m])  *legendrePlex1(l,x[m]) 
            M23[ii,jj] = +1/k1II *psixx (l,kII*a,x[m]) *legendrePl(l,x[m]) 
            M24[ii,jj] = -1/k1I  *xixx (l,kI*a,x[m])   *legendrePl(l,x[m]) 
    
            M31[ii,jj] = +1/k2II *psixx (l,kII*a,x[m]) *legendrePlex1(l,x[m]) 
            M32[ii,jj] = -1/k2I  *xixx (l,kI*a,x[m])   *legendrePlex1(l,x[m]) 
            M33[ii,jj] = M11[ii,jj] 
            M34[ii,jj] = M12[ii,jj] 
    
            M41[ii,jj] = +1/k2II *psixx (l,kII*a,x[m]) *legendrePl(l,x[m]) 
            M42[ii,jj] = -1/k2I  *xixx (l,kI*a,x[m])   *legendrePl(l,x[m]) 
            M43[ii,jj] = M21[ii,jj] 
            M44[ii,jj] = M22[ii,jj] 
    
            M51[ii,jj] = + k1II*(psixx(l,kII*a,x[m])+psi2xx(l,kII*a,x[m])) *legendrePl(l,x[m]) 
            M52[ii,jj] = - k1I *(xixx(l,kI*a,x[m])+xi2xx(l,kI*a,x[m]))     *legendrePl(l,x[m]) 
            M53[ii,jj] = 0 
            M54[ii,jj] = 0 
    
            M61[ii,jj] = 0
            M62[ii,jj] = 0 
            M63[ii,jj] = + MuII *(psixx(l,kII*a,x[m])+psi2xx(l,kII*a,x[m])) *legendrePl(l,x[m])  
            M64[ii,jj] = - MuI  *(xixx(l,kI*a,x[m])+xi2xx(l,kI*a,x[m]))     *legendrePl(l,x[m]) 

    Matrix = np.zeros((6*mmax, 6*lmax), dtype=complex)
    
    for ii in range(mmax):
        m = ii + 1
        for jj in range(lmax):
            l = jj + 1
            Matrix[ii,jj]                 = M11[ii,jj]
            Matrix[ii,jj+lmax]            = M12[ii,jj]
            Matrix[ii,jj+2*lmax]          = M13[ii,jj]
            Matrix[ii,jj+3*lmax]          = M14[ii,jj]

            Matrix[ii+mmax,jj]            = M21[ii,jj]
            Matrix[ii+mmax,jj+lmax]       = M22[ii,jj]
            Matrix[ii+mmax,jj+2*lmax]     = M23[ii,jj]
            Matrix[ii+mmax,jj+3*lmax]     = M24[ii,jj]
            
            Matrix[ii+2*mmax,jj]          = M31[ii,jj]
            Matrix[ii+2*mmax,jj+lmax]     = M32[ii,jj]
            Matrix[ii+2*mmax,jj+2*lmax]   = M33[ii,jj]
            Matrix[ii+2*mmax,jj+3*lmax]   = M34[ii,jj]
    
            Matrix[ii+3*mmax,jj]          = M41[ii,jj]
            Matrix[ii+3*mmax,jj+lmax]     = M42[ii,jj]
            Matrix[ii+3*mmax,jj+2*lmax]   = M43[ii,jj]
            Matrix[ii+3*mmax,jj+3*lmax]   = M44[ii,jj]
    
            Matrix[ii+4*mmax,jj]          = M51[ii,jj]
            Matrix[ii+4*mmax,jj+lmax]     = M52[ii,jj]
            Matrix[ii+4*mmax,jj+2*lmax]   = M53[ii,jj]
            Matrix[ii+4*mmax,jj+3*lmax]   = M54[ii,jj]
    
            Matrix[ii+5*mmax,jj]          = M61[ii,jj]
            Matrix[ii+5*mmax,jj+lmax]     = M62[ii,jj]
            Matrix[ii+5*mmax,jj+2*lmax]   = M63[ii,jj]
            Matrix[ii+5*mmax,jj+3*lmax]   = M64[ii,jj]

    VectorL = np.zeros(6*mmax, dtype=complex)
    VectorR = np.zeros(6*mmax, dtype=complex)
    for ii in range(mmax):
        # left and right
        VectorL[ii]        = N1L[ii]
        VectorL[ii+mmax]   = N2L[ii]
        VectorL[ii+2*mmax] = N3L[ii]
        VectorL[ii+3*mmax] = N4L[ii]
        VectorL[ii+4*mmax] = N5L[ii]
        VectorL[ii+5*mmax] = N6L[ii]
    
        VectorR[ii]        = N1R[ii]
        VectorR[ii+mmax]   = N2R[ii]
        VectorR[ii+2*mmax] = N3R[ii]
        VectorR[ii+3*mmax] = N4R[ii]
        VectorR[ii+4*mmax] = N5R[ii]
        VectorR[ii+5*mmax] = N6R[ii]

    Weight = np.zeros(6*mmax, dtype=complex)
    # weights 
    for ii in range(mmax):
        m = ii + 1
        Weight[ii]        = 1   
        Weight[ii+mmax]   = 1   
        Weight[ii+2*mmax] = 1
        Weight[ii+3*mmax] = 1            
        Weight[ii+4*mmax] = 1/(100*m**0.8)  
        Weight[ii+5*mmax] = 1/(22*m**0.20)  

    # Weights were commented out?:
    # xL = lscov (Matrix,VectorL.',Weight.').'
    xL = lscov(np.matrix(Matrix), np.matrix(VectorL).T, np.matrix(Weight).T)
    xR = lscov(np.matrix(Matrix), np.matrix(VectorR).T, np.matrix(Weight).T)

    if verbose:
        print('If Eps=0, ignore previous error messages regarding rank deficiency!')
    eCcL = np.zeros(lmax, dtype=complex)
    eDdL = np.zeros(lmax, dtype=complex)
    mCcL = np.zeros(lmax, dtype=complex)
    mDdL = np.zeros(lmax, dtype=complex)
    
    eCcR = np.zeros(lmax, dtype=complex)
    eDdR = np.zeros(lmax, dtype=complex)
    mCcR = np.zeros(lmax, dtype=complex)
    mDdR = np.zeros(lmax, dtype=complex)
    
    for jj in range(lmax):
        l = jj + 1
        # left and right first order expansion coefficients
        eCcL[jj] = xL[jj]
        eDdL[jj] = xL[jj+lmax]
        mCcL[jj] = xL[jj+2*lmax]
        mDdL[jj] = xL[jj+3*lmax]
    
        eCcR[jj] = xR[jj]        *(-1)**(l-1)     
        eDdR[jj] = xR[jj+lmax]   *(-1)**(l-1)
        mCcR[jj] = xR[jj+2*lmax] *(-1)**l
        mDdR[jj] = xR[jj+3*lmax] *(-1)**l
    
        # corrected expansion coefficients for right optical fiber*  
        eBR[jj] = eBR[jj] *(-1)**(l-1)
        eCR[jj] = eCR[jj] *(-1)**(l-1)
        eDR[jj] = eDR[jj] *(-1)**(l-1)
        mBR[jj] = mBR[jj] *(-1)**l
        mCR[jj] = mCR[jj] *(-1)**l
        mDR[jj] = mDR[jj] *(-1)**l
        #*additional factors [(-1)**(l-1),(-1)**l] should have been included when eB, mB are calculated however due to least squares approach slight antisymmetry arises...
        # hence include factors only after least square approach  

    ## Field Expansions For Incident Fields
    # sums of radial, zenithal and azimuthal field terms from l=1 to l=lmax 
    # initialisation of sums 
    # Paul: These are not used?!
    #EriL   = lambda r, th, phi: 0
    #EriR   = lambda r, th, phi: 0 
    #EthiL  = lambda r, th, phi: 0
    #EthiR  = lambda r, th, phi: 0 
    #EphiiL = lambda r, th, phi: 0
    #EphiiR = lambda r, th, phi: 0 
    #HriL   = lambda r, th, phi: 0
    #HriR   = lambda r, th, phi: 0 
    #HthiL  = lambda r, th, phi: 0
    #HthiR  = lambda r, th, phi: 0 
    #HphiiL = lambda r, th, phi: 0
    #HphiiR = lambda r, th, phi: 0 
    # left
    #EriL   = lambda r, th, phi: EriL(r,th,phi)   + eBL[l-1]*(psi(l,kI*r)+psi2(l,kI*r))*legendrePl(l,np.cos(th))*np.cos(phi)         
    #EthiL  = lambda r, th, phi: EthiL(r,th,phi)  -(eBL[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                                        + mBL[l-1]*psi (l,kI*r)/(k1I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)    
    #EphiiL = lambda r, th, phi: EphiiL(r,th,phi) -(eBL[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                                        + mBL[l-1]*psi (l,kI*r)/(k1I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.sin(phi)   
    #HriL   = lambda r, th, phi: HriL(r,th,phi)   + mBL[l-1]*(psi(l,kI*r)+psi2(l,kI*r))*legendrePl(l,np.cos(th))*np.sin(phi)         
    #HthiL  = lambda r, th, phi: HthiL(r,th,phi)  -(mBL[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                                        + eBL[l-1]*psi (l,kI*r)/(k2I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi)    
    #HphiiL = lambda r, th, phi: HphiiL(r,th,phi) +(mBL[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                                        + eBL[l-1]*psi (l,kI*r)/(k2I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.cos(phi)   
    # right
    #EriR   = lambda r, th, phi: EriR(r,th,phi)   + eBR[l-1]*(psi(l,kI*r)+psi2(l,kI*r))*legendrePl(l,np.cos(th))*np.cos(phi)         
    #EthiR  = lambda r, th, phi: EthiR(r,th,phi)  -(eBR[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                                        + mBR[l-1]*psi (l,kI*r)/(k1I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)    
    #EphiiR = lambda r, th, phi: EphiiR(r,th,phi) -(eBR[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                                        + mBR[l-1]*psi (l,kI*r)/(k1I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.sin(phi)   
    #HriR   = lambda r, th, phi: HriR(r,th,phi)   + mBR[l-1]*(psi(l,kI*r)+psi2(l,kI*r))*legendrePl(l,np.cos(th))*np.sin(phi)         
    #HthiR  = lambda r, th, phi: HthiR(r,th,phi)  -(mBR[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                                        + eBR[l-1]*psi (l,kI*r)/(k2I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi)    
    #HphiiR = lambda r, th, phi: HphiiR(r,th,phi) +(mBR[l-1]*psi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                                        + eBR[l-1]*psi (l,kI*r)/(k2I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.cos(phi)   

    ## Field Expansions For Scattered Fields
    # sums of radial, zenithal and azimuthal field terms from l=1 to l=lmax 
    # Paul: workaround to recursing into lambda functions
    def wrapper_expansion(lambda_func):
        def wrapped(r, th, ph):
            result = 0
            for jj in range(lmax):
                l = jj + 1
                result += lambda_func(r, th, ph, l)
            return result
        return wrapped
    # left
    ErsL_it   = lambda r,th,phi,l: (eDL[l-1] + eDdL[l-1])*(xi(l,kI*r)+xi2(l,kI*r))*legendrePl (l,np.cos(th))*np.cos(phi)
    EthsL_it  = lambda r,th,phi,l: -((eDL[l-1] + eDdL[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl1 (l,np.cos(th))*np.sin(th) \
                                  +(mDL[l-1] + mDdL[l-1])*xi (l,kI*r)/(k1I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)   
    EphisL_it = lambda r,th,phi,l: -((eDL[l-1] + eDdL[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
                                  +(mDL[l-1] + mDdL[l-1])*xi (l,kI*r)/(k1I*r)*legendrePl1 (l,np.cos(th))*np.sin(th))*np.sin(phi) 
    HrsL_it   = lambda r,th,phi,l: +(mDL[l-1] + mDdL[l-1])*(xi(l,kI*r)+xi2(l,kI*r))*legendrePl (l,np.cos(th))*np.sin(phi)         
    HthsL_it  = lambda r,th,phi,l: -((mDL[l-1] + mDdL[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl1 (l,np.cos(th))*np.sin(th) \
                                  +(eDL[l-1] + eDdL[l-1])*xi (l,kI*r)/(k2I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi)   
    HphisL_it = lambda r,th,phi,l: +((mDL[l-1] + mDdL[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
                                  +(eDL[l-1] + eDdL[l-1])*xi (l,kI*r)/(k2I*r)*legendrePl1 (l,np.cos(th))*np.sin(th))*np.cos(phi) 
    # right
    ErsR_it   = lambda r,th,phi,l: +(eDR[l-1] + eDdR[l-1])*(xi(l,kI*r)+xi2(l,kI*r))*legendrePl(l,np.cos(th))*np.cos(phi)         
    EthsR_it  = lambda r,th,phi,l: -((eDR[l-1] + eDdR[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
                                  +(mDR[l-1] + mDdR[l-1])*xi (l,kI*r)/(k1I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)   
    EphisR_it = lambda r,th,phi,l: -((eDR[l-1] + eDdR[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
                                  +(mDR[l-1] + mDdR[l-1])*xi (l,kI*r)/(k1I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.sin(phi)  
    HrsR_it   = lambda r,th,phi,l: +(mDR[l-1] + mDdR[l-1])*(xi(l,kI*r)+xi2(l,kI*r))*legendrePl(l,np.cos(th))*np.sin(phi)         
    HthsR_it  = lambda r,th,phi,l: -((mDR[l-1] + mDdR[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
                                  +(eDR[l-1] + eDdR[l-1])*xi (l,kI*r)/(k2I*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi) 
    HphisR_it = lambda r,th,phi,l: +((mDR[l-1] + mDdR[l-1])*xi1 (l,kI*r)/(kI*r)*legendrePl (l,np.cos(th))/np.sin(th) \
                                  +(eDR[l-1] + eDdR[l-1])*xi (l,kI*r)/(k2I*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.cos(phi)

    ErsL = wrapper_expansion(ErsL_it)
    EthsL = wrapper_expansion(EthsL_it)
    EphisL = wrapper_expansion(EphisL_it)
    HrsL = wrapper_expansion(HrsL_it)
    HthsL = wrapper_expansion(HthsL_it)
    HphisL = wrapper_expansion(HphisL_it)
    
    ErsR = wrapper_expansion(ErsR_it)
    EthsR = wrapper_expansion(EthsR_it)
    EphisR = wrapper_expansion(EphisR_it)
    HrsR = wrapper_expansion(HrsR_it)
    HthsR = wrapper_expansion(HthsR_it)
    HphisR = wrapper_expansion(HphisR_it)

    ## Field Expansions For Internal Fields (Within)
    # sums of radial, zenithal and azimuthal field terms from l=1 to l=lmax 
    # initialisation of sums
    ErwL_it   = lambda r,th,phi,l: +(eCL[l-1] + eCcL[l-1])*(psi(l,kII*r)+psi2(l,kII*r))*legendrePl(l,np.cos(th))*np.cos(phi)      
    EthwL_it  = lambda r,th,phi,l: -((eCL[l-1] + eCcL[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
                                  + (mCL[l-1] + mCcL[l-1])*psi (l,kII*r)/(k1II*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)  
    EphiwL_it = lambda r,th,phi,l: -((eCL[l-1] + eCcL[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl (l,np.cos(th))/np.sin(th) \
                                  + (mCL[l-1] + mCcL[l-1])*psi (l,kII*r)/(k1II*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.sin(phi) 
    # Paul: These are not used?!
    #HrwL_it   = lambda r,th,phi,l: +(mCL[l-1] + mCcL[l-1])*(psi(l,kII*r)+psi2(l,kII*r))*legendrePl(l,np.cos(th))*np.sin(phi)       
    #HthwL_it  = lambda r,th,phi,l: -((mCL[l-1] + mCcL[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                              + (eCL[l-1] + eCcL[l-1])*psi (l,kII*r)/(k2II*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi)
    #HphiwL_it = lambda r,th,phi,l: +((mCL[l-1] + mCcL[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                              + (eCL[l-1] + eCcL[l-1])*psi (l,kII*r)/(k2II*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.cos(phi)
    ## right
    #ErwR_it   = lambda r,th,phi,l: +(eCR[l-1] + eCcR[l-1])*(psi(l,kII*r)+psi2(l,kII*r))*legendrePl(l,np.cos(th))*np.cos(phi)       
    #EthwR_it  = lambda r,th,phi,l: -((eCR[l-1] + eCcR[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                              + (mCR[l-1] + mCcR[l-1])*psi (l,kII*r)/(k1II*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.cos(phi)
    #EphiwR_it = lambda r,th,phi,l: -((eCR[l-1] + eCcR[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                              + (mCR[l-1] + mCcR[l-1])*psi (l,kII*r)/(k1II*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.sin(phi)
    #HrwR_it   = lambda r,th,phi,l: +(mCR[l-1] + mCcR[l-1])*(psi(l,kII*r)+psi2(l,kII*r))*legendrePl(l,np.cos(th))*np.sin(phi)       
    #HthwR_it  = lambda r,th,phi,l: -((mCR[l-1] + mCcR[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl1(l,np.cos(th))*np.sin(th) \
    #                              + (eCR[l-1] + eCcR[l-1])*psi (l,kII*r)/(k2II*r)*legendrePl (l,np.cos(th))/np.sin(th))*np.sin(phi)
    #HphiwR_it = lambda r,th,phi,l: +((mCR[l-1] + mCcR[l-1])*psi1 (l,kII*r)/(kII*r)*legendrePl (l,np.cos(th))/np.sin(th) \
    #                              + (eCR[l-1] + eCcR[l-1])*psi (l,kII*r)/(k2II*r)*legendrePl1(l,np.cos(th))*np.sin(th))*np.cos(phi)

    ErwL = wrapper_expansion(ErwL_it)
    EthwL = wrapper_expansion(EthwL_it)
    EphiwL = wrapper_expansion(EphiwL_it)
    # Paul: These are not used?!
    #HrwL = wrapper_expansion(HrwL_it)
    #HthwL = wrapper_expansion(HthwL_it)
    #HphiwL = wrapper_expansion(HphiwL_it)
    #
    #ErwR = wrapper_expansion(ErwR_it)
    #EthwR = wrapper_expansion(EthwR_it)
    #EphiwR = wrapper_expansion(EphiwR_it)
    #HrwR = wrapper_expansion(HrwR_it)
    #HthwR = wrapper_expansion(HthwR_it)
    #HphiwR = wrapper_expansion(HphiwR_it)

    ## Overall Fields
    # sum of the electric fields in region I (independent of time) [separated for left and right]
    ErL   = lambda r, th, phi: ErsL(r,th,phi)   + eEriL(r,th,phi)
    ErR   = lambda r, th, phi: ErsR(r,th,phi)   + eEriR(r,th,phi)
    EthL  = lambda r, th, phi: EthsL(r,th,phi)  + eEthiL(r,th,phi)
    EthR  = lambda r, th, phi: EthsR(r,th,phi)  + eEthiR(r,th,phi)
    EphiL = lambda r, th, phi: EphisL(r,th,phi) + eEphiiL(r,th,phi)
    EphiR = lambda r, th, phi: EphisR(r,th,phi) + eEphiiR(r,th,phi)
    HrL   = lambda r, th, phi: HrsL(r,th,phi)   + eHriL(r,th,phi)
    HrR   = lambda r, th, phi: HrsR(r,th,phi)   + eHriR(r,th,phi)
    HthL  = lambda r, th, phi: HthsL(r,th,phi)  + eHthiL(r,th,phi)
    HthR  = lambda r, th, phi: HthsR(r,th,phi)  + eHthiR(r,th,phi)
    HphiL = lambda r, th, phi: HphisL(r,th,phi) + eHphiiL(r,th,phi)
    HphiR = lambda r, th, phi: HphisR(r,th,phi) + eHphiiR(r,th,phi)

    ## Checking Boundary Conditions 
    # generate points on boundary
    r = np.zeros(22, dtype=complex)
    th = np.zeros(22, dtype=complex)
    phi = np.pi/4
    for ii in range(22):
        th[ii]   = np.pi/20.003*(ii+1-0.999)
        r[ii]    = a*B1(np.cos(th[ii]))

    FEr1 = np.zeros(22, dtype=float)
    FEr2 = np.zeros(22, dtype=float)
    FEth1 = np.zeros(22, dtype=float)
    FEth2 = np.zeros(22, dtype=float) 
    FEphi1 = np.zeros(22, dtype=float)
    FEphi2 = np.zeros(22, dtype=float)
    
    FEr = np.zeros(22, dtype=float)
    FEth = np.zeros(22, dtype=float)
    FEphi = np.zeros(22, dtype=float)

    for ii in range(22):
        # internal and external fields 
        FEr1[ii]    = np.real(EpsilonI *  ErL(r[ii],th[ii],phi))
        FEr2[ii]    = np.real(EpsilonII* ErwL(r[ii],th[ii],phi))
        FEth1[ii]   = np.real( EthL(r[ii],th[ii],phi))
        FEth2[ii]   = np.real(EthwL(r[ii],th[ii],phi))
        FEphi1[ii]  = np.real( EphiL(r[ii],th[ii],phi))
        FEphi2[ii]  = np.real(EphiwL(r[ii],th[ii],phi))
        # ratios
        FEr[ii]     = FEr1[ii]/FEr2[ii]
        FEth[ii]    = FEth1[ii]/FEth2[ii]
        FEphi[ii]   = FEphi1[ii]/FEphi2[ii]

    if verbose:
        print('Following output indicates the accuracy of the matching of the boundary conditions.')
        print('theta, Er, Er_ratio, Eth, Eth_ratio, Ephi, Ephi_ratio.')
        print(np.array([th, FEr1, FEr, FEth1, FEth, FEphi1, FEphi]))
        print('The ratio should be close to 1 for good matching of the BCs. However if the fields tend to zero the ratio is not that relevant.')

    ## Data for Plotting
    # time-averaged, radial stress on infinitesimal area dA = r**2 np.sin(th) dr dth dphi
    sigmarr = lambda r, th, phi: -1/(16*np.pi)*np.real(
                EpsilonI*(ErL(r,th,phi)*np.conj(ErL(r,th,phi))-EthL(r,th,phi)*np.conj(EthL(r,th,phi))-EphiL(r,th,phi)*np.conj(EphiL(r,th,phi))) \
              + EpsilonI*(ErR(r,th,phi)*np.conj(ErR(r,th,phi))-EthR(r,th,phi)*np.conj(EthR(r,th,phi))-EphiR(r,th,phi)*np.conj(EphiR(r,th,phi))) \
              + MuI*     (HrL(r,th,phi)*np.conj(HrL(r,th,phi))-HthL(r,th,phi)*np.conj(HthL(r,th,phi))-HphiL(r,th,phi)*np.conj(HphiL(r,th,phi))) \
              + MuI*     (HrR(r,th,phi)*np.conj(HrR(r,th,phi))-HthR(r,th,phi)*np.conj(HthR(r,th,phi))-HphiR(r,th,phi)*np.conj(HphiR(r,th,phi))))

    th = np.zeros(numpoints, dtype=float)
    r = np.zeros(numpoints, dtype=float)
    sigma = np.zeros(numpoints, dtype=float)
    
    for ii in range(numpoints):
        th[ii]     = theta_max/(numpoints - 0.998) *(ii+1-0.999)
        r[ii]      = a*B1(np.cos(th[ii]))
        sigma[ii]  = sigmarr(r[ii],th[ii],np.pi/4)
        

    res = [th, sigma]
    
    if ret_legendre_decomp:
        coeff = stress_legendre_decomp(th=th,
                                       sigmarr=sigma,
                                       n_poly=lmax)
        res.append(coeff)

    return res


def stress_legendre_decomp(th, sigmarr, n_poly):
    """Decompose stress into Legendre Polynomials"""
    # Sigma = Sum_n [Coeff(n) P_n(np.cos(theta))]
    nmax = n_poly                    # number of Legendre polynomials used in fit

    # transfer data from stress plot into pair of corresponding variables
    # [Theta,Sigma]
    numpoints = th.shape[0]
    theta = th.reshape(-1, 1)
    sigma = sigmarr.reshape(-1, 1)

    # Write set of linear equations for stresses in terms Legendre functions
    legmat = np.zeros((numpoints,nmax), dtype=float)
    for ii in range(numpoints):
        for jj in np.arange(nmax)[::2]: # skip odd Legendre Polynomials since stress is an even function (symmetrical)
            legmat[ii, jj] = np.real_if_close(legendrePlm(0, jj, np.cos(theta[ii])))

    coeff = lscov(legmat, sigma)

    return coeff
