import sys
import FittingUtilities

from scipy.interpolate import InterpolatedUnivariateSpline as spline
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u, constants

import HelperFunctions
import Broaden
import StellarModel
import triangle
import os
from astropy.io import fits
from math import floor
import GenericSearch

# ########################################################################
# ########################################################################


def ErrorFunction2(params, data, model_getter):
    model_orders = MakeModel(params, data, model_getter)
    loglikelihood = []
    N = 0.0
    for i, (o, m) in enumerate(zip(data, model_orders)):
        ratio = o.y / m
        cont = FittingUtilities.Continuum(o.x, ratio, fitorder=5, lowreject=2, highreject=2)
        # cont = np.poly1d(np.polyfit(o.x, ratio, 5))(o.x)
        #cont = o.cont
        loglikelihood.append((o.y - cont * m) / o.err)
        N += o.size()
        data[i].cont = cont
    loglikelihood = np.hstack(loglikelihood)
    print "RV = {:g}, vsini = {:g} --> X^2 = {:g}".format(params['rv'].value, params['vsini'].value,
                                                          np.sum(loglikelihood ** 2) / N)
    return loglikelihood


def ErrorFunction(params, data, model_getter):
    model_orders = MakeModel(params, data, model_getter)
    loglikelihood = []
    N = 0.0
    for i, (o, m) in enumerate(zip(data, model_orders)):
        ratio = o.y / m
        cont = FittingUtilities.Continuum(o.x, ratio, fitorder=5, lowreject=2, highreject=2)
        #cont = np.poly1d(np.polyfit(o.x, ratio, 5))(o.x)

        loglikelihood.append((o.y - o.cont * m) / o.err)
        N += o.size()
    loglikelihood = np.hstack(loglikelihood)
    for key in params.keys():
        print "{:s} = {:.10f}".format(key, params[key].value)
    print "X^2 = {:g}\n\n".format(np.sum(loglikelihood ** 2) / N)
    #print "RV = {:g}, vsini = {:g} --> X^2 = {:g}".format(params['rv'].value, params['vsini'].value, np.sum(loglikelihood**2)/N)
    return loglikelihood


def MakeModel(pars, data, model_getter):
    vsini = pars['vsini'].value * u.km.to(u.cm)
    rv = pars['rv'].value
    T = pars['temperature'].value
    logg = pars['logg'].value
    metal = pars['metal'].value
    alpha = pars['alpha'].value

    c = constants.c.cgs.to(u.km / u.s).value

    #Get the model from the ModelGetter instance
    model = model_getter(T, logg, metal, alpha)

    # Next, broaden the model
    broadened = Broaden.RotBroad(model, vsini, linear=True)
    modelfcn = spline(broadened.x, broadened.y)
    model_orders = []
    print rv /c
    for order in data:
        model_orders.append(modelfcn(order.x * (1 - rv / c)))
    return model_orders


def LM_Model(x, vsini, rv, temperature, logg, metal, alpha, model_getter=None):
    if model_getter is None:
        raise KeyError("Must give model_getter keyword!")
    c = constants.c.cgs.to(u.km / u.s).value
    vsini *= u.km.to(u.cm)

    # Get the model from the ModelGetter instance
    model = model_getter(temperature, logg, metal, alpha)

    # Next, broaden the model
    broadened = Broaden.RotBroad(model, vsini, linear=True)

    # Finally, interpolate to the x values
    modelfcn = spline(broadened.x, broadened.y)
    return modelfcn(x * (1.0 - rv / c))


if __name__ == "__main__":
    # Define some constants
    c = constants.c
    good_orders = range(41)
    good_orders.pop(33)

    # Set up default values, and then read in command line arguments
    T_min = 9000
    T_max = 9250
    metal_min = -0.8
    metal_max = 0.5
    logg_min = 3.0
    logg_max = 4.5
    alpha_min = 0.0
    alpha_max = 0.4
    modeldir = "models/"
    rv = 0.0 * u.km / u.s
    vsini = 200.0 * u.km / u.s
    R = 80000.0
    output_dir = "{:s}/Dropbox/School/Research/AstarStuff/Parameters/".format(os.environ['HOME'])
    texfile = "Parameters.tex"  # Put in output_dir after command-line arguments are parsed
    file_list = []
    debug = False
    N_iter = 100

    for arg in sys.argv[1:]:
        if "temp" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            T_min = float(values[0])
            T_max = float(values[-1])
        elif "metal" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            metal_min = float(values[0])
            metal_max = float(values[-1])
        elif "logg" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            logg_min = float(values[0])
            logg_max = float(values[-1])
        elif "alpha" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            alpha_min = float(values[0])
            alpha_max = float(values[-1])
        elif "model" in arg.lower():
            modeldir = arg.partition("=")[-1]
            if not modeldir.endswith("/"):
                modeldir += "/"
        elif "rv" in arg.lower():
            rv = float(arg.partition("=")[-1]) * u.km / u.s
        elif "vsini" in arg.lower():
            vsini = float(arg.partition("=")[-1]) * u.km / u.s
        elif "resolution" in arg.lower():
            R = float(arg.partition("=")[-1])
        elif "outdir" in arg.lower():
            output_dir = arg.partition("=")[-1]
            if not output_dir.endswith("/"):
                output_dir += "/"
        elif "texfile" in arg.lower():
            texfile = arg.partition("=")[-1]
        elif "debug" in arg.lower():
            debug = True
            print "Debug mode ON"
        elif "iteration" in arg.lower():
            N_iter = int(arg.partition("=")[-1])
        else:
            file_list.append(arg)

    texfile = "{:s}{:s}".format(output_dir, texfile)

    # Make sure files were give
    if len(file_list) == 0:
        sys.exit("Must give at least one file!")

    # Make guess values for each of the values from the bounds
    temperature = (T_min + T_max) / 2.0
    logg = (logg_min + logg_max) / 2.0
    metal = (metal_min + metal_max) / 2.0
    alpha = (alpha_min + alpha_max) / 2.0

    #Make an instance of the model getter
    mg = StellarModel.KuruczGetter(modeldir,
                                   T_min=T_min,
                                   T_max=T_max,
                                   logg_min=logg_min,
                                   logg_max=logg_max,
                                   metal_min=metal_min,
                                   metal_max=metal_max,
                                   alpha_min=alpha_min,
                                   alpha_max=alpha_max,
                                   wavemin=350.0)

    # Make the appropriate lmfit model
    fitter = HelperFunctions.ListModel(LM_Model, independent_vars=['x'], model_getter=mg)

    #Set default values
    fitter.set_param_hint("rv", value=rv.value, min=-50, max=50)
    fitter.set_param_hint('vsini', value=vsini.value, vary=True, min=0.0, max=500.0)
    fitter.set_param_hint('temperature', value=temperature, min=T_min, max=T_max, vary=True)
    fitter.set_param_hint('logg', value=logg, min=logg_min, max=logg_max, vary=True)
    fitter.set_param_hint('metal', value=metal, min=metal_min, max=metal_max, vary=True)
    fitter.set_param_hint('alpha', value=alpha, min=alpha_min, max=alpha_max, vary=True)

    """
    Here is the main loop over files!
    """
    for filename in file_list:
        print "Fitting parameters for {}".format(filename)
        # Read the data
        all_orders = HelperFunctions.ReadExtensionFits(filename)
        orders = [o[1] for o in enumerate(all_orders) if o[0] in good_orders]

        # Perform the initial fit
        optdict = {"epsfcn": 1e-2}
        params = fitter.make_params()
        result = fitter.fit(orders, fit_kws=optdict, params=params)

        print(result.fit_report())
        if debug:
            for i, order in enumerate(orders):
                m = result.best_fit[i]
                ratio = order.y / m
                order.cont = FittingUtilities.Continuum(order.x, ratio, lowreject=2, highreject=2, fitorder=5)
                plt.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
                plt.plot(order.x, result.best_fit[i], 'r-', alpha=0.5)
            plt.show()


        # Now, re-do the fit several times to get bootstrap error estimates
        fitparams = {"rv": np.zeros(N_iter),
                     "vsini": np.zeros(N_iter),
                     "temperature": np.zeros(N_iter),
                     "logg": np.zeros(N_iter),
                     "metal": np.zeros(N_iter),
                     "alpha": np.zeros(N_iter)}
        params = result.params
        orders_original = [o.copy() for o in orders]
        chainfile = open("chain_temp.dat", "w")
        for n in range(N_iter):
            print "Fitting iteration {:d}/{:d}".format(n + 1, N_iter)
            orders = []
            for order in orders_original:
                o = order.copy()
                o.y += np.random.normal(loc=0, scale=o.err)
                orders.append(o.copy())
            result = fitter.fit(orders, fit_kws=optdict, params=params)
            if debug:
                print "\n**********     Best values      ************"
            for key in fitparams.keys():
                fitparams[key][n] = result.best_values[key]
                chainfile.write("{:g}\t".format(result.best_values[key]))
                if debug:
                    print key, ': ', result.best_values[key]
            print "\n\n"
            chainfile.write("\n")
        chainfile.close()

        # Correct the velocity for barycentric motion
        header = fits.getheader(filename)
        vbary = GenericSearch.HelCorr(header, observatory="CTIO")
        fitparams['rv'] += vbary

        # Save the fitted parameters
        texlog = open(texfile, "a")
        header = fits.getheader(filename)
        date = header['date'].split("T")[0]
        star = header['object']
        texlog.write("{:s} & {:s}".format(star, date))
        print "\n\nBest-fit parameters:\n================================="
        for key in fitparams.keys():
            low, med, up = np.percentile(fitparams[key], [16, 50, 84])
            up_err = up - med
            low_err = med - low
            # Get significant digits
            dist = max(-int(floor(np.log10(up_err))), -int(floor(np.log10(low_err))))
            med = round(med, dist)
            up_err = round(up_err, dist + 1)
            low_err = round(low_err, dist + 1)
            print "{:s} = {:g} + {:g} / - {:g}".format(key, med, up_err, low_err)
            texlog.write(" & $%g^{+ %g}_{- %g$}" % (med, up_err, low_err))
        texlog.write(" \\\\ \n")

        # Save the full results in a directory labeled by the star name and date
        stardir = "{:s}{:s}/".format(output_dir, star.replace(" ", "_"))
        HelperFunctions.ensure_dir(stardir)
        datedir = "{:s}{:s}/".format(stardir, date)
        HelperFunctions.ensure_dir(datedir)
        chain_filename = "{:s}chain.dat".format(datedir)
        chain = np.vstack([fitparams[key] for key in fitparams.keys()]).T
        print "Outputting chain to {:s}".format(chain_filename)
        np.savetxt(chain_filename, chain)
        fig, axes = plt.subplots(len(fitparams), len(fitparams), figsize=(10, 10))
        labeldict = {'rv': '$ \\rm rv$ $ \\rm (km \\cdot s^{-1}$)',
                     'vsini': '$ \\rm v \sin{i}$ $ \\rm (km s^{-1}$)',
                     'temperature': '$ \\rm T_{\\rm eff}$ $\\rm (K)$',
                     'logg': '$\log{g}$',
                     'metal': '$\\rm [Fe/H]$',
                     'alpha': '$\\rm [\\alpha/Fe]$'}
        names = [labeldict[key] for key in fitparams.keys()]
        triangle.corner(chain, labels=names, fig=fig)
        plt.savefig("{:s}corner_plot.pdf".format(datedir))
        # plt.show()

        print "Done with file {:s}\n\n\n".format(filename)






