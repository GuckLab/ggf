from lmfit.models import LinearModel
import matplotlib.pylab as plt
import numpy as np

from scan import server_args

for arg in server_args:
    ids = []
    dims = []
    for aset in arg:
        ids.append("{}-{:.3e}-{:.3e}-{:05d}".format(aset[0][:3],
                                                    aset[1],
                                                    aset[2],
                                                    aset[3]))
    idset = "_".join(ids)
    data = np.load("_results_{}.npy".format(idset))
    
    aset = arg[0]
    vals = np.linspace(*aset[1:])
    ax = plt.subplot(111)
    ax.plot(vals, data, "o")
    ax.set_xlabel(aset[0])
    ax.set_ylabel("GGF")

    # lmfit
    try:
        y = data
        x = vals
        mod = LinearModel(nan_policy="omit")
        pars = mod.guess(y, x=x)
        out = mod.fit(y, pars, x=x)
        ax.plot(x, out.eval(params=out.params, x=x), "-")
    except:
        pass


    if np.nanmin(data) > .6 and np.nanmax(data) < .9:
        ax.set_ylim(.6, .9)


    plt.grid()
    plt.tight_layout()
    plt.savefig("_plot_{}.png".format(idset), dpi=100)
    plt.close()
    

    