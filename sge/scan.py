import sys

from scan_manager import PM_Client, PM_Server


case = 1

if case == 1:
    # default parameters
    kwargs = dict(model="boyde2009",
                  stretch_ratio=.06,
                  semi_minor=8e-6,
                  relative_object_index=1.0319101123595504,
                  medium_index=1.335,
                  effective_fiber_distance=180e-6,
                  mode_field_diameter=4.8e-6,
                  power_per_fiber=.65,
                  wavelength=780e-9,
                  poisson_ratio=0.5,
                  n_poly=None,
                  verbose=False)

    # changed parameters (each item will produce a new scan)
    server_args = [[["stretch_ratio", 0, .15, 50]],
                   [["semi_minor", 6.5e-6, 8e-6, 50]],
                   [["relative_object_index", 1.01, 1.10, 50]],
                   [["medium_index", 1.33, 1.43, 20]],
                   [["effective_fiber_distance", 170e-6, 190e-6, 10]],
                   [["mode_field_diameter", 3e-6, 6e-6, 10]],
                   [["power_per_fiber", .5, .8, 10]],
                   [["wavelength", 700e-9, 820e-9, 10]],
                   [["poisson_ratio", 0.4, .5, 10]],
                   [["n_poly", 20, 140, 10]],
                   ]


if __name__ == "__main__":
    mode = sys.argv[-1]
    
    if mode == "server":
        with PM_Server(server_args, defaults=kwargs) as server:
            server.start()
    elif mode == "client":
        client = PM_Client()
        client.start()
    else:
        print("Unknown mode:", mode)
        print("Specify 'server' or 'client'")

