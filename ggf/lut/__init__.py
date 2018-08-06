"""Look-up table for GGF computation"""
import pathlib
from pkg_resources import resource_filename

import h5py
import numpy as np
from scipy.interpolate import interpn


class NotInLUTError(BaseException):
    pass


def _match_lut(kwargs, lut_path=None):
    """Find matching LUT for a set of keyword arguments

    Parameters
    ----------
    kwargs: dict
        keyword arguments for ggf.get_ggf
    lut_path: str, pathlib.Path or None
        path to an hdf5 file containing a LUT; if None,
        the default LUTs will be searched

    Returns
    -------
    lut_path: pathlib.Path
        path to matching LUT hdf5 file
    """
    if lut_path:
        paths = [pathlib.Path(lut_path)]
    else:
        paths = get_lut_paths()
    for path in paths:
        with h5py.File(path, mode="r") as h5:
            for key in kwargs:
                if _param_in_lut(key, kwargs[key], h5["lut"].attrs):
                    pass  # everything OK so far
                else:
                    break  # try next LUT
            else:
                break  # this is the right LUT
    else:
        msg = "No matching LUT found for: {} in {}".format(kwargs, paths)
        raise NotInLUTError(msg)
    return path


def _param_in_lut(key, value, h5_attrs):
    """Test whether the given parameter is accessible via a LUT

    Parameters
    ----------
    key: str
        keyword name (see ggf.get_ggf)
    value: str or float
        keyword value to check
    h5_attrs: dict
        attributes from a LUT hdf5 file

    Returns
    -------
    inlut: bool
        True, if the key-value pair is covered by the LUT
    """
    inlut = False
    if key == "n_poly" and value is None:  # default value is given in LUT
        inlut = True
    elif key in h5_attrs:  # fixed parameter
        if h5_attrs[key] == value or value is None:
            inlut = True
    elif isinstance(value, str):  # model
        inlut = False
    else:  # range
        vmin = h5_attrs["{} min".format(key)]
        vmax = h5_attrs["{} max".format(key)]
        if value <= vmax and value >= vmin:
            inlut = True
        else:
            inlut = False
    return inlut


def get_lut_paths():
    """Return a list of look-up table hdf5 files in ggf"""
    lutpath = pathlib.Path(resource_filename("ggf", "lut"))
    paths = lutpath.glob("*.h5")
    return sorted(paths)


def get_ggf_lut(model, semi_major, semi_minor, object_index, medium_index,
                effective_fiber_distance, mode_field_diameter,
                power_per_fiber, wavelength, poisson_ratio,
                n_poly=None, lut_path=None, verbose=False):
    """Linear interpolation of the GGF from a look-up table

    Parameters
    ----------
    model: str
        Model to use, one of: `boyde2009`
    semi_major: float
        Semi-major axis of an ellipse fit to the object perimeter [m]
    semi_minor: float
        Semi-minor axis of an ellipse fit to the object perimeter [m]
    object_index: float
        Refractive index of the object
    medium_index: float
        Refractive index of the surrounding medium
    effective_fiber_distance: float
        Effective distance between the two trapping fibers relative
        to the medium refractive index [m]. For an open setup, this is
        the physical distance between the fibers. For a closed setup
        (capillary), this distance takes into account the refractive
        indices and thicknesses of the glass capillary and index
        matching gel. For the closed setup, the convenience function
        :func:`ggf.fiber_distance_capillary` can be used.
    mode_field_diameter: float
        The mode field diameter MFD of the fiber used [m]. Note that
        the MFD is dependent on the wavelength used. If the
        manufacturer did not provide a value for the MFD, the MFD
        can be approximated as ``3*wavelenth`` for a single-mode
        fiber.
    power_per_fiber: float
        The laser power coupled into each of the fibers [W]
    wavelength: float
        The laser wavelength used for the trap [m]
    poisson_ratio: float
        The Poisson's ratio of the stretched material. Set this
        to 0.5 for volume conservation.
    n_poly: int
        Number of Legendre polynomials to use for computing the GGF.
        Note that only even Legendre polynomials are used and thus,
        this number is effectively halved.
    lut_path: str or pathlib.Path
        Path to a LUT hdf5 file. If `None`, the internal LUTs are
        used.
    verbose: int
        Increases verbosity

    Returns
    -------
    ggf: float
        Linearly interpolated global geometric factor

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
      - semi_minor: linear only within interval of 0.08µm
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
    # determine the correct LUT
    kwargs = {"model": model,
              "stretch_ratio": stretch_ratio,
              "semi_minor": semi_minor,
              "relative_object_index": relative_object_index,
              "medium_index": medium_index,
              "effective_fiber_distance": effective_fiber_distance,
              "mode_field_diameter": mode_field_diameter,
              "power_per_fiber": power_per_fiber,
              "wavelength": wavelength,
              "poisson_ratio": poisson_ratio,
              "n_poly": n_poly}
    # get the right LUT
    lut_path = _match_lut(kwargs, lut_path)
    if verbose:
        print("Using LUT path: {}".format(lut_path))
    # get LUT data
    with h5py.File(lut_path, mode="r") as h5:
        values = h5["lut"].value
        meta = dict(h5["lut"].attrs)
    # order of interpolation dimensions
    order = meta["dimension_order"].split(",")
    # grid points
    points = []
    for label in order:
        points.append(np.linspace(meta["{} min".format(label)],
                                  meta["{} max".format(label)],
                                  meta["{} num".format(label)]))
    # interpolation coordinates
    xi = [kwargs[kk] for kk in order]
    # perform interpolation
    ggf = interpn(points=points, values=values, xi=xi)
    return np.asscalar(ggf)
