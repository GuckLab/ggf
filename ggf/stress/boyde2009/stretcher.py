"""Convenience functions for optical stretcher computations"""
import numpy as np

from .core import stress



def distance_capillary(gel_dist=2e-6, glass_dist=40e-6, medium_dist=40e-6,
                       gel_index=1.449, glass_index=1.474, medium_index=1.335):
    """Effective distance between optical fiber and channel center
    
    When the optical stretcher is combined with a microfluidic
    channel ("closed setup"), then the effective distance between
    the optical fiber and the channel center (location of the
    stretched object) is defined by the refractive indices of
    the optical components: index matching gel between fiber and
    channel wall, microfluidic glass channel wall, and medium
    inside the channel.
    
    Parameters
    ----------
    gel_dist: float
        Thickness of index matching gel (distance between fiber and glass wall) [m]
    glass_dist: float
        Thickness of glass wall [m]
    medium_dist: float
        Distance between glass wall (side that is in contact with cell medium) and cell center [m]
    gel_index: float
        Refractive index of index matching gel
    glass_index: float
        Refractive index of channel glass wall
    medium_index: float
        Refractive index of index medium inside channel
    
    Returns
    -------
    eff_dist: float
        Effective distance between fiber and channel center
    
    Notes
    -----
    The effective distance is computed relative to the medium,
    i.e. if `gel_index` == `glass_index` == `medium_index`, then
    `eff_dist` = `gel_dist` + `glass_dist` + `medium_dist`.
    """
    eff_dist = medium_index / gel_index * gel_dist \
               + medium_index / glass_index * glass_dist \
               + medium_dist
    return eff_dist


def semiax_to_radrat(semi_major, semi_minor, poisson_ratio=.45):
    """Convert semi-major and semi-minor axes to radius and stretch ratio
    
    Notes
    -----
    Example: a stretched spheroidal cell with
    stretch_ratio=0.1 and poisson_ratio=0.5 has 
    semi-minor axes: b = c = radius (1-poisson_ratio stretch_ratio) = 0.95 radius
    semi-major axis:     a = radius (1+stretch_ratio) = 1.10 radius

    """
    stretch_ratio = (semi_major - semi_minor) / (semi_minor + poisson_ratio * semi_major)
    radius = semi_major / (1+ stretch_ratio)

    assert np.allclose(semi_major, radius*(1+stretch_ratio))
    assert np.allclose(semi_minor, radius*(1-poisson_ratio*stretch_ratio))
    return radius, stretch_ratio


def stress_capillary(semi_major, semi_minor, object_index,
                     fiber_core_radius=2.2e-6, wavelength=780e-9, 
                     poisson_ratio=.45, 
                     gel_dist=2e-6, glass_dist=40e-6, medium_dist=40e-6,
                     gel_index=1.449, glass_index=1.474, medium_index=1.335,
                     power_left=.6, power_right=.6,
                     numpoints=100, theta_max=np.pi,
                     ret_legendre_decomp=False,
                     verbose=False):                        

    dist = distance_capillary(gel_dist=gel_dist,
                              glass_dist=glass_dist,
                              medium_dist=medium_dist,
                              gel_index=gel_index,
                              glass_index=glass_index,
                              medium_index=medium_index)
    
    return stress_open(semi_major=semi_major,
                       semi_minor=semi_minor,
                       object_index=object_index,
                       poisson_ratio=poisson_ratio,
                       medium_index=medium_index,
                       wavelength=wavelength,
                       dist_object_fiber=dist,
                       fiber_core_radius=fiber_core_radius,
                       power_left=power_left,
                       power_right=power_right,
                       numpoints=numpoints,
                       theta_max=theta_max,
                       ret_legendre_decomp=ret_legendre_decomp,
                       verbose=verbose)


def stress_open(semi_major, semi_minor, object_index, poisson_ratio=.45,
                medium_index=1.335, wavelength=780e-9,
                dist_object_fiber=100e-6, fiber_core_radius=2.2e-6,
                power_left=.6, power_right=.6,
                numpoints=100, theta_max=np.pi,
                ret_legendre_decomp=False,
                verbose=False):

    radius, stretch_ratio = semiax_to_radrat(semi_major=semi_major,
                                             semi_minor=semi_minor,
                                             poisson_ratio=poisson_ratio)

    beam_waist = fiber_core_radius / wavelength

    return stress(object_index=object_index,
                  medium_index=medium_index,
                  radius=radius,
                  poisson_ratio=poisson_ratio,
                  stretch_ratio=stretch_ratio,
                  wavelength=wavelength,
                  beam_waist=beam_waist,
                  power_left=power_left,
                  power_right=power_right,
                  dist=dist_object_fiber,
                  numpoints=numpoints,
                  theta_max=theta_max,
                  ret_legendre_decomp=ret_legendre_decomp,
                  verbose=verbose)
