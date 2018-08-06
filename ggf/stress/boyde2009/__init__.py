from . import core


def get_stress(semi_major, semi_minor, object_index, medium_index,
               effective_fiber_distance=100e6, mode_field_diameter=4.8e-6,
               power_per_fiber=.6, wavelength=1064e-9, n_points=100,
               verbose=False):
    """Wraps around core.stress without poisson_ratio

    See :func:`ggf.stress.get_stress` for parameter descriptions.
    """
    return core.stress(
        object_index=object_index,
        medium_index=medium_index,
        poisson_ratio=0,
        semi_minor=semi_minor,
        stretch_ratio=(semi_major - semi_minor) / (semi_minor),
        wavelength=wavelength,
        beam_waist=mode_field_diameter / 2 / wavelength,
        power_left=power_per_fiber,
        power_right=power_per_fiber,
        dist=effective_fiber_distance / 2,
        n_points=n_points,
        verbose=verbose)
