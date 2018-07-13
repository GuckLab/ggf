import sys

from lut_manager import PM_Client, PM_Server


luts = {
    "guck_open_cell":
        {"model": "boyde2009",
         "stretch_ratio": [0, .15, 40],
         "semi_minor": [6.5e-6, 8e-6, 15],
         "relative_object_index": [1.02, 1.03, 5],
         "medium_index": [1.333, 1.338, 2],
         "effective_fiber_distance": [170e-6, 190e-6, 2],
         "mode_field_diameter": 4.8e-6,
         "power_per_fiber": .65,
         "wavelength": 780e-9,
         "poisson_ratio": [0.4, .5, 2],
         "n_poly": 120,
         },

#    "guck_open_fus":
#        {"model": "boyde2009",
#         "stretch_ratio": [0, .15, 40],
#         "semi_minor": [6.5e-6, 8e-6, 15],
#         "relative_object_index": [1.00001, 1.10, 50],
#         "medium_index": [1.340, 1.347, 3],
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
