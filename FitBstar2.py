import sys
import FittingUtilities
import numpy as np
import os
from math import floor

from scipy.interpolate import InterpolatedUnivariateSpline as spline
import matplotlib.pyplot as plt
from astropy import units as u, constants
import HelperFunctions
import Broaden
import StellarModel
import triangle
from astropy.io import fits
import GenericSearch


def LM_Model(x, vsini, rv, temperature, logg, metal, alpha, model_getter=None, **mgargs):
    if model_getter is None:
        raise KeyError("Must give model_getter keyword!")
    c = constants.c.cgs.to(u.km / u.s).value
    # vsini *= u.km.to(u.cm)

    # Get the model from the ModelGetter instance
    mgargs['vsini'] = vsini
    broadened = model_getter(temperature, logg, metal, alpha, **mgargs)

    # Next, broaden the model
    #broadened = Broaden.RotBroad(model, vsini, linear=True)

    # Finally, interpolate to the x values
    modelfcn = spline(broadened.x, broadened.y)
    return modelfcn(x * (1.0 - rv / c))


def Fit(arguments, mg=None):
    # Define some constants
    c = constants.c
    good_orders = range(41)
    good_orders.pop(33)

    # Set up default values, and then read in command line arguments
    T_min = 9000
    T_max = 9250
    metal_min = -0.8
    metal_max = 0.2
    logg_min = 3.0
    logg_max = 5.0
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

    for arg in arguments:
        if "-temp" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            T_min = float(values[0])
            T_max = float(values[-1])
        elif "-metal" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            metal_min = float(values[0])
            metal_max = float(values[-1])
        elif "-logg" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            logg_min = float(values[0])
            logg_max = float(values[-1])
        elif "-alpha" in arg.lower():
            r = arg.partition("=")[-1]
            values = r.partition(",")
            alpha_min = float(values[0])
            alpha_max = float(values[-1])
        elif "-model" in arg.lower():
            modeldir = arg.partition("=")[-1]
            if "," in modeldir:
                #More than one directory is given
                modeldir = modeldir.split(",")
                for m in modeldir:
                    if not m.endswith("/"):
                        m += "/"
            elif not modeldir.endswith("/"):
                modeldir += "/"
        elif "-rv" in arg.lower():
            rv = float(arg.partition("=")[-1]) * u.km / u.s
        elif "-vsini" in arg.lower():
            vsini = float(arg.partition("=")[-1]) * u.km / u.s
        elif "-resolution" in arg.lower():
            R = float(arg.partition("=")[-1])
        elif "-outdir" in arg.lower():
            output_dir = arg.partition("=")[-1]
            if not output_dir.endswith("/"):
                output_dir += "/"
        elif "-texfile" in arg.lower():
            texfile = arg.partition("=")[-1]
        elif "-debug" in arg.lower():
            debug = True
            print "Debug mode ON"
        elif "-iteration" in arg.lower():
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
    if mg is None:
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
    fitter.set_param_hint('alpha', value=0.0, min=alpha_min, max=alpha_max, vary=mg.alpha_varies)

    """
    Here is the main loop over files!
    """
    for filename in file_list:
        # Make output directories
        header = fits.getheader(filename)
        date = header['date'].split("T")[0]
        star = header['object']
        stardir = "{:s}{:s}/".format(output_dir, star.replace(" ", "_"))
        HelperFunctions.ensure_dir(stardir)
        datedir = "{:s}{:s}/".format(stardir, date)
        HelperFunctions.ensure_dir(datedir)
        chain_filename = "{:s}chain.dat".format(datedir)

        # Read the data
        print "Fitting parameters for {}".format(filename)
        all_orders = HelperFunctions.ReadExtensionFits(filename)
        orders = [o[1] for o in enumerate(all_orders) if o[0] in good_orders]

        # Perform the fit
        optdict = {"epsfcn": 1e-2}
        params = fitter.make_params()
        fitparams = {"rv": np.zeros(N_iter),
                     "vsini": np.zeros(N_iter),
                     "temperature": np.zeros(N_iter),
                     "logg": np.zeros(N_iter),
                     "metal": np.zeros(N_iter),
                     "alpha": np.zeros(N_iter)}
        orders_original = [o.copy() for o in orders]
        chainfile = open(chain_filename, "w")
        vbary = GenericSearch.HelCorr_IRAF(header, observatory="CTIO")
        for n in range(N_iter):
            print "Fitting iteration {:d}/{:d}".format(n + 1, N_iter)
            orders = []
            for order in orders_original:
                o = order.copy()
                o.y += np.random.normal(loc=0, scale=o.err)
                orders.append(o.copy())

            # Make a fast interpolator instance if not the first loop
            result = fitter.fit(orders, fit_kws=optdict, params=params)
            result.best_values['rv'] += vbary
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

        # Save the fitted parameters
        texlog = open(texfile, "a")
        texlog.write("{:s} & {:s}".format(star, date))
        print "\n\nBest-fit parameters:\n================================="
        if mg.alpha_varies:
            keys = ['rv', 'temperature', 'metal', 'vsini', 'logg', 'alpha']
        else:
            keys = ['rv', 'temperature', 'metal', 'vsini', 'logg']
        for key in keys:
            low, med, up = np.percentile(fitparams[key], [16, 50, 84])
            up_err = up - med
            low_err = med - low
            # Get significant digits
            dist = max(-int(floor(np.log10(up_err))), -int(floor(np.log10(low_err))))
            med = round(med, dist)
            up_err = round(up_err, dist + 1)
            low_err = round(low_err, dist + 1)
            print "{:s} = {:g} + {:g} / - {:g}".format(key, med, up_err, low_err)
            texlog.write(" & $%g^{+ %g}_{- %g}$" % (med, up_err, low_err))
        if not mg.alpha_varies:
            texlog.write(" & $0.0^{+0.0}_{-0.0}$")
        texlog.write(" \\\\ \n")

        # Save a corner plot of the fitted results
        chain = np.vstack([fitparams[key] for key in fitparams.keys()]).T
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

        print "Done with file {:s}\n\n\n".format(filename)


if __name__ == "__main__":
    Fit(sys.argv[1:])
