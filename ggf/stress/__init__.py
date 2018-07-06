from . import boyde2009
from .geometry import fiber_distance_capillary

VALID_MODELS = ["boyde2009"]


def get_stress(model, semi_major, semi_minor, object_index, medium_index,
               effective_fiber_distance=100e-6, mode_field_diameter=3e-6,
               power_per_fiber=.6, wavelength=1064e-9, n_points=100,
               verbose=False):
    """
    
    """
    if model not in VALID_MODELS:
        msg = "`model` must be one of {}, got '{}'".format(VALID_MODELS, model)
        raise ValueError(msg)

    if model == "boyde2009":
        func = boyde2009.get_stress
    
    return func(
        semi_major=semi_major,
        semi_minor=semi_minor,
        object_index=object_index,
        medium_index=medium_index,
        effective_fiber_distance=effective_fiber_distance,
        mode_field_diameter=mode_field_diameter,
        power_per_fiber=power_per_fiber,
        wavelength=wavelength,
        n_points=n_points,
        verbose=verbose)
