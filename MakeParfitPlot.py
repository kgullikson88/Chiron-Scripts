import pandas
import sys
import os
from scipy.interpolate import InterpolatedUnivariateSpline as spline

import FittingUtilities
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import constants, units as u
import numpy as np

import HelperFunctions
import StellarModel


def make_models(chainfile, orders, model_getter=None, modeldir=None):
    names = ['rv', 'temperature', 'metal', 'vsini', 'logg', 'alpha']
    chain = pandas.read_csv(chainfile, names=names, sep="\t", header=False, skipinitialspace=True)
    n = len(chain)
    c = constants.c.cgs.to(u.km / u.s).value

    if model_getter is None:
        # Need to make a model_getter instance
        if modeldir is None:
            sys.exit("Must give either a model_getter instance or a model directory!")
        T_min = min(chain['temperature']) - 1000
        T_max = max(chain['temperature']) + 1000
        logg_min = min(chain['logg']) - 0.5
        logg_max = max(chain['logg']) + 0.5
        metal_min = min(chain['metal']) - 0.5
        metal_max = max(chain['metal']) + 0.5
        alpha_min = 0.0
        alpha_max = 0.4
        model_getter = StellarModel.KuruczGetter(modeldir,
                                                 T_min=T_min,
                                                 T_max=T_max,
                                                 logg_min=logg_min,
                                                 logg_max=logg_max,
                                                 metal_min=metal_min,
                                                 metal_max=metal_max,
                                                 alpha_min=alpha_min,
                                                 alpha_max=alpha_max,
                                                 wavemin=350.0)

    # Make a model for chain
    models = [[] for order in orders]
    for i in range(n):
        temperature = chain.iloc[i].temperature
        logg = chain.iloc[i].logg
        metal = chain.iloc[i].metal
        alpha = chain.iloc[i].alpha
        vsini = chain.iloc[i].vsini
        rv = chain.iloc[i].rv
        model = model_getter(temperature, logg, metal, alpha, vsini=vsini)
        modelfcn = spline(model.x, model.y)
        model_orders = []
        for j, order in enumerate(orders):
            m = modelfcn(order.x * (1.0 - rv / c))
            #model_orders.append(m)
            models[j].append(m)
            #models.append(np.array(model_orders))
    return models


if __name__ == "__main__":
    fileList = []
    par_dir = "{:s}/Dropbox/School/Research/AstarStuff/Parameters/".format(os.environ['HOME'])
    modeldir = None
    for arg in sys.argv[1:]:
        if "-model" in arg:
            modeldir = arg.split("=")[1]
        else:
            fileList.append(arg)

    for fname in fileList:
        orders = HelperFunctions.ReadExtensionFits(fname)
        header = fits.getheader(fname)
        star = header['object']
        date = header['date']
        stardir = "{:s}{:s}/".format(par_dir, star.replace(" ", "_"))
        datedir = "{:s}{:s}/".format(stardir, date)
        chain_filename = "{:s}chain.dat".format(datedir)

        if os.path.isfile(chain_filename):
            models = make_models(chain_filename, orders, modeldir=modeldir)

            n = len(models)
            for i, order in enumerate(orders):
                median_model = np.median(models[i], axis=0)
                ratio = order.y / median_model
                order.cont = FittingUtilities.Continuum(order.x, ratio, fitorder=5, lowreject=2, highreject=2)
                plt.plot(order.x, order.y / order.cont, 'r-', lw=2)
                for model in models:
                    plt.plot(order.x, model, 'k-', alpha=0.1)
            plt.show()
