"""Used for quantifying LUT error by randomly sampling kwargs"""
import pathlib

import h5py
import numpy as np
import percache

import ggf

NUM_ERRS = 1000  # number of random error values to compute

mycache = percache.Cache("lut_test.cache", livesync=True)


@mycache
def compute_ggf(**kwargs):
    return ggf.get_ggf(use_lut=False, **kwargs)


def get_kwarg_ranges(lut_path):
    with h5py.File(lut_path, mode="r") as h5:
        attrs = dict(h5["lut"].attrs)

    fixed = {}
    ranges = {}

    for kw in ["model", "stretch_ratio", "semi_minor", "relative_object_index",
               "medium_index", "effective_fiber_distance",
               "mode_field_diameter", "power_per_fiber",
               "wavelength", "poisson_ratio", "n_poly"]:
        if kw in attrs:
            fixed[kw] = attrs[kw]
        else:
            ranges[kw] = (attrs["{} min".format(kw)],
                          attrs["{} max".format(kw)])

    return fixed, ranges


def map_lut2geom(kwargs):
    kw = kwargs.copy()
    relative_object_index = kw.pop("relative_object_index")
    kw["object_index"] = relative_object_index * kw["medium_index"]

    stretch_ratio = kw.pop("stretch_ratio")
    kw["semi_major"] = stretch_ratio * kw["semi_minor"] + kw["semi_minor"]
    return kw


# get lut paths
paths = pathlib.Path(__file__).parent.glob("*.h5")
paths = sorted(paths)[::-1]


for path in paths:
    # make everything reproducible
    np.random.set_state(np.random.RandomState(42).get_state())
    fixed, ranges = get_kwarg_ranges(path)
    errs = []
    # dice out a few parameters
    for ii in range(NUM_ERRS):
        kwargs = fixed.copy()
        for rr in ranges:
            kwargs[rr] = np.random.uniform(low=ranges[rr][0],
                                           high=ranges[rr][1])
        kw_ggf = map_lut2geom(kwargs)
        # get LUT value
        ggf1 = ggf.get_ggf(use_lut=path, **kw_ggf)
        # compute value (cached)
        ggf2 = compute_ggf(**kw_ggf)
        erri = (ggf2-ggf1)/ggf2
        errs.append(erri)
        print(path, "{:03d} error: {:.1f}%".format(ii, erri*100))
        np.savetxt(path.with_suffix(".errors.txt"), errs)
