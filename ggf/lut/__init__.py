"""Look-up table for GGF computation

Keyword arguments are treated as follows:
- semi_major axis is converted to stretch ratio
- object_index is normalized by medium_index

- No interpolation: wavelength, mode_field_diameter
- Linear dependency (only 2 points are computed): 



"""
import pathlib
from pkg_resources import resource_filename

import h5py


def get_lut_paths():
    lutpath = pathlib.Path(resource_filename("ggf", "lut"))
    paths = lutpath.glob("*.h5")
    return sorted(paths)



def get_ggf_lut(model, semi_major, semi_minor, object_index, medium_index,
                effective_fiber_distance, mode_field_diameter,
                power_per_fiber, wavelength, poisson_ratio,
                n_poly=None, verbose=False):
    """Linear interpolation of the GGF from a look-up table
    
    Notes
    -----
    - To avoid invalid values in the look-up table (LUT), such as
      `semi_major < semi_minor` or `object_index < medium_index`,
      the LUT is not built using the exact same keyword arguments
      as this method:
    
      - `object_index` is stored as
        ``relative_object_index = object_index / medium_index``
      - `semi_major` is stored as
        ``stretch_ratio = (semi_major - semi_minor) / semi_minor``
    - The following keywords are not interpolated in the LUT:
      - `model`
      - `wavelength`: the OS uses a fixed wavelength
      - `mode_field_diameter`: the fiber geometry is fixed
      - `power_per_fiber`: usually fixed for reproducibility
      - `n_poly`: set to a high number (e.g. 120)
    - The following are approximate guiding values for when a keyword
      can be considered linear:
      - stretch_ratio: linear only within interval of 0.004
      - semi_minor: linear only within interval of 0.05µm
      - relative_object_index: linear only within interval of 0.003
      - medium_index: linear only within interval of 0.005
      - poisson_ratio: good linearity
      - power_per_fiber: good linearity
      - effective_fiber_distance: linear only within interval of 15µm
    """
    # convert major_axis to stretch ratio
    stretch_ratio = (semi_major - semi_minor) / semi_minor
    # normalize object index with medium_index
    relative_object_index = object_index / medium_index
    # determine the correct look-up table
