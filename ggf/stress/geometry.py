

def fiber_distance_capillary(gel_thickness=2e-6, glass_thickness=40e-6,
                             channel_width=40e-6, gel_index=1.449,
                             glass_index=1.474, medium_index=1.335):
    """Effective distance between the two optical fibers

    When the optical stretcher is combined with a microfluidic
    channel (closed setup), then the effective distance between
    the two optical fibers (with the the stretched object at the
    channel center) is defined by the refractive indices of
    the optical components: index matching gel between fiber and
    channel wall, microfluidic glass channel wall, and medium
    inside the channel.

    Parameters
    ----------
    gel_thickness: float
        Thickness of index matching gel (distance between fiber
        and glass wall) [m]
    glass_thickness: float
        Thickness of glass wall [m]
    channel_width: float
        Width of the microfluidic channed [m]
    gel_index: float
        Refractive index of index matching gel
    glass_index: float
        Refractive index of channel glass wall
    medium_index: float
        Refractive index of index medium inside channel

    Returns
    -------
    eff_dist: float
        Effective distance between the fibers

    Notes
    -----
    The effective distance is computed relative to the medium,
    i.e. if `gel_index` == `glass_index` == `medium_index`, then
    `eff_dist` = 2*`gel_dist` + 2*`glass_dist` + `channel_width`.
    """
    eff_dist = 2 * medium_index / gel_index * gel_thickness \
        + 2 * medium_index / glass_index * glass_thickness \
        + channel_width
    return eff_dist
