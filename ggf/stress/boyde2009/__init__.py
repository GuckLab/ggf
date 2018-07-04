from . import core


def get_stress(semi_major, semi_minor, object_index, medium_index,
               dist_object_fiber=100e6, beam_waist_radius=3e-6, 
               power_per_fiber=.6, wavelength=1064e-9,
               verbose=False):
    """Wraps around core.stress without poisson_ratio"""
    return core.stress(object_index=object_index,
                       medium_index=medium_index,
                       poisson_ratio=0,
                       radius=semi_minor,
                       stretch_ratio=semi_major/semi_minor,
                       wavelength=wavelength,
                       beam_waist=beam_waist_radius/wavelength,
                       power_left=power_per_fiber,
                       power_right=power_per_fiber,
                       dist=dist_object_fiber,
                       verbose=verbose)
