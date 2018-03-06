"""Computation of the global geometric factor"""
import numpy as np

from .matlab_funcs import lscov, legendre
from .stretcher import stress_capillary, stress_open


def coeff2ggf(coeff, poisson_ratio=.45):
    """Compute the global geometric factor from stress coefficients

    The original Matlab script states that the computed GGF
    already includes the peak stress:

    .. math::
        
        \text{GGF} = \sigma_0 F_\text{G}.
    """
    ## Parameters used in the program
    
    Theta1 =0
    # Nu =input('enter the value of nu (Poisson s ratio) \n')
    Nu = poisson_ratio
    ##

    SigmaMat= coeff
    n=SigmaMat.size-1
    SigmaMatamended = SigmaMat/2 # amendement due to the correction of the paper by Ananthakrishnan
    Sigmatot=np.zeros((2*n+2,1))
    
    for k in range(n):
        Sigmatot[k]=SigmaMatamended[k]
        #Sigmatot[k+n+1] = 0
    
    Mat=np.zeros((2*n+2,2*n+2))
    r=1# r is set to 1 for the calculation. 
    N=n
    
    for d in range(n+1):
        Mat[d,d]=(r**d)*(d + 1)*(d**2-d-2-2*Nu)# 
        Mat[d,d+N+1]=(r**(d - 2))*d*(d - 1)# 

    for d in range(2,n):
        Mat[d+N,d+1]=(r**(d+1))*(-1+2*(d+1)+(d+1)**2+2*Nu)*(-(d+1))*(d+2)/(1+2*(d+1)) #  i pour le coefficient et j pour la ligne de sigma pour le coefficient P(d)
        Mat[d+N,d-1]=(r**(d-1))*(-1+2*(d-1)+((d-1)**2)+2*Nu)*(d)*(d-1)/(1+2*(d-1))
        Mat[d+N,d+N+1+1]=(r**(d-1))*(d)*(-(d+1))*(d+2)/(1+2*(d+1))
        Mat[d+N,d+N]=(r**(d-3))*(d-2)*(d)*(d-1)/(1+2*(d-1))
    
    d=0
    Mat[d+N,d+1]=(r**(d+1))*(-1+2*(d+1)+(d+1)**2+2*Nu)*(-(d+1))*(d+2)/(1+2*(d+1))
    Mat[d+N,d+N+1+1]=(r**(d-1))*(d)*(-(d+1))*(d+2)/(1+2*(d+1))
   
    d=1
    Mat[d+N,d+1]=(r**(d+1))*(-1+2*(d+1)+((d+1)**2)+2*Nu)*(-(d+1))*(d+2)/(1+2*(d+1))
    Mat[d+N,d+N+1+1]=(r**(d-1))*(d)*(-(d+1))*(d+2)/(1+2*(d+1))

    d=n
    Mat[d+N,d-1]=(r**(d-1))*(-1+2*(d-1)+((d-1)**2)+2*Nu)*(d)*(d-1)/(1+2*(d-1))
    Mat[d+N,d+N]=(r**(d-3))*(d-2)*(d)*(d-1)/(1+2*(d-1))
    
    d=n+1
    Mat[d+N,d-1]=(r**(d-1))*(-1+2*(d-1)+((d-1)**2)+2*Nu)*(d)*(d-1)/(1+2*(d-1))
    Mat[d+N,d+N]=(r**(d-3))*(d-2)*(d)*(d-1)/(1+2*(d-1))

    Mat1=Mat
    sol=lscov(Mat1,Sigmatot)
    
    PL = np.zeros(n+1, dtype=complex)
    w0 = np.zeros(n+1, dtype=complex)

    ##  Calculation of the Global Geometrical Factor including sigma0
    ww=0
    for k in range(0,n):
        PLA=legendre(k,np.cos(Theta1))
        PL[k]=PLA[0][0] # legendre polynomial with m=0!
        w0[k]=((sol[k])*(r**(k+1))*(k+1 )*(k-2 +4*Nu) +(sol[n+k+1])*r**(k-1)*k)*PL[k]
        ww=ww+w0[k]

    GF=ww

    return GF


def ggf_capillary(semi_major, semi_minor, object_index,
                  fiber_core_radius=2.2e-6, wavelength=780e-9, 
                  poisson_ratio=.45, 
                  gel_dist=2e-6, glass_dist=40e-6, medium_dist=40e-6,
                  gel_index=1.449, glass_index=1.474, medium_index=1.335,
                  power_left=.6, power_right=.6,
                  numpoints=100, theta_max=np.pi,
                  verbose=False):
    th, sigma, coeff = stress_capillary(semi_major=semi_major,
                                        semi_minor=semi_minor,
                                        object_index=object_index,
                                        fiber_core_radius=fiber_core_radius,
                                        wavelength=wavelength,
                                        poisson_ratio=poisson_ratio,
                                        gel_dist=gel_dist,
                                        glass_dist=glass_dist,
                                        medium_dist=medium_dist,
                                        gel_index=gel_index,
                                        glass_index=glass_index,
                                        medium_index=medium_index,
                                        power_left=power_left,
                                        power_right=power_right,
                                        numpoints=numpoints,
                                        theta_max=theta_max,
                                        ret_legendre_decomp=True,
                                        verbose=verbose)

    return coeff2ggf(coeff, poisson_ratio=poisson_ratio)


def ggf_open(semi_major, semi_minor, object_index, poisson_ratio=.45,
             medium_index=1.335, wavelength=780e-9,
             dist_object_fiber=100e-6, fiber_core_radius=2.2e-6,
             power_left=.6, power_right=.6,
             numpoints=100, theta_max=np.pi,
             verbose=False):
    
    th, sigma, coeff = stress_open(semi_major=semi_major,
                                   semi_minor=semi_minor,
                                   object_index=object_index,
                                   poisson_ratio=poisson_ratio,
                                   medium_index=medium_index,
                                   wavelength=wavelength,
                                   dist_object_fiber=dist_object_fiber,
                                   fiber_core_radius=fiber_core_radius,
                                   power_left=power_left,
                                   power_right=power_right,
                                   numpoints=numpoints,
                                   theta_max=theta_max,
                                   ret_legendre_decomp=True,
                                   verbose=verbose)
    
    return coeff2ggf(coeff, poisson_ratio=poisson_ratio)
