import sys

from manager import PM_Client, PM_Server


if __name__ == "__main__":
    # default parameters
    kwargs = dict(model="boyde2009",
                  semi_major=10e-6,
                  semi_minor=10e-6,
                  object_index=1.3776,
                  medium_index=1.335,
                  effective_fiber_distance=180e-6,
                  mode_field_diameter=2.4e-6,
                  power_per_fiber=.65,
                  wavelength=780e-9,
                  poisson_ratio=0.5,
                  n_poly=None,
                  verbose=False)

    # changed parameters (each item will produce a new scan)
    server_args = [[["semi_major", 5e-6, 20e-6, 10]],
                   [["object_index", 1.34, 1.39, 10]],
                   [["medium_index", 1.33, 1.34, 10]],
                   [["effective_fiber_distance", 150e-6, 200e-6, 10]],
                   [["mode_field_diameter", 2e-6, 3e-6, 10]],
                   [["power_per_fiber", .5, .8, 10]],
                   [["wavelength", 700e-9, 820e-9, 10]],
                   [["poisson_ratio", 0.4, .5, 10]],
                   [["n_poly", 20, 200, 10]],
                   ]
    
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

