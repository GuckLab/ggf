from . import boyde2009
from .geometry import distance_capillary

VALID_MODELS = ["boyde2009"]


def get_stress(model, semi_major, semi_minor, object_index, medium_index,
               dist_object_fiber=100e-6, beam_waist_radius=3e-6,
               power_per_fiber=.6, wavelength=1064e-9, verbose=False):
    """
    
    """
    if model not in VALID_MODELS:
        msg = "`model` must be one of {}, got '{}'".format(VALID_MODELS, model)
        raise ValueError(msg)

    if model == "boyde2009":
        return boyde2009.get_stress(semi_major=semi_major,
                                    semi_minor=semi_minor,
                                    object_index=object_index,
                                    medium_index=medium_index,
                                    dist_object_fiber=dist_object_fiber,
                                    beam_waist_radius=beam_waist_radius,
                                    power_per_fiber=power_per_fiber,
                                    wavelength=wavelength,
                                    verbose=verbose)