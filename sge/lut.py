"""LUT computation on a distributed system"""
import sys

from lut_manager import PM_Client, PM_Server


# Each top-level key defines a different LUT. The keys within
# the top-level dictionaries define the range of parameters used
# to compute LUT. A float means 1-dimentsional. A List defines
# the spacing (`np.linspace`). E.g. '"stretch_ratio": [0, .20, 47],'
# means that stretch ratios are computed for all values in
# `np.linspce(1, .2, 47, endpoint=True)`. To get an idea how large
# the spacing should be, you can plot 1D line-scans through the LUT
# with scan.py and plot_scan.py.
luts = {
    "guck_open_fus":
        {"model": "boyde2009",
         "stretch_ratio": [0, .20, 63],
         "semi_minor": [2.5e-6, 8.5e-6, 17],
         "relative_object_index": [1.009, 1.042, 23],
         "medium_index": [1.340, 1.347, 3],
         "effective_fiber_distance": [175e-6, 215e-6, 2],
         "mode_field_diameter": 4.8e-6,
         "power_per_fiber": [.5, .7, 2],
         "wavelength": 780e-9,
         "poisson_ratio": [0.4, .5, 2],
         "n_poly": 120,
         },
    "guck_all_cell":
        {"model": "boyde2009",
         "stretch_ratio": [0, .13, 33],
         "semi_minor": [5e-6, 10e-6, 30],
         "relative_object_index": [1.015, 1.035, 7],
         "medium_index": [1.333, 1.338, 2],
         "effective_fiber_distance": [158e-6, 175e-6, 2],
         "mode_field_diameter": 4.8e-6,
         "power_per_fiber": [0.5, 1.0, 2],
         "wavelength": 780e-9,
         "poisson_ratio": [0.4, .5, 2],
         "n_poly": 120,
         },
    #    "guck_open_cell":
    #        {"model": "boyde2009",
    #         "stretch_ratio": [0, .15, 40],
    #         "semi_minor": [6.5e-6, 8e-6, 15],
    #         "relative_object_index": [1.02, 1.03, 5],
    #         "medium_index": [1.333, 1.338, 2],
    #         "effective_fiber_distance": [170e-6, 190e-6, 2],
    #         "mode_field_diameter": 4.8e-6,
    #         "power_per_fiber": .65,
    #         "wavelength": 780e-9,
    #         "poisson_ratio": [0.4, .5, 2],
    #         "n_poly": 120,
    #         },
}


if __name__ == "__main__":
    mode = sys.argv[-1]

    if mode == "server":
        with PM_Server(server_args=luts) as server:
            server.start()
    elif mode == "client":
        client = PM_Client()
        client.start()
    else:
        print("Unknown mode:", mode)
        print("Specify 'server' or 'client'")
