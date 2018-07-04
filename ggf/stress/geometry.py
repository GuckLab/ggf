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
