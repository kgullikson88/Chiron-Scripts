import sys
import FittingUtilities

from scipy.interpolate import InterpolatedUnivariateSpline as spline
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u, constants
import emcee

import HelperFunctions
import Broaden
import StellarModel





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
    print vsini, rv, temperature, logg, metal, alpha, c
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
        else:
            filename = arg





    #Read the data
    orders = HelperFunctions.ReadExtensionFits(filename)

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
                                   wavemin=orders[0].x[0] - 1)

    #Make guess values for each of the values from the bounds
    temperature = (T_min + T_max) / 2.0
    logg = (logg_min + logg_max) / 2.0
    metal = (metal_min + metal_max) / 2.0
    alpha = (alpha_min + alpha_max) / 2.0

    # Make the appropriate lmfit model
    fitter = HelperFunctions.ListModel(LM_Model, independent_vars=['x'], model_getter=mg)

    #Set default values
    fitter.set_param_hint("rv", value=rv.value, min=-50, max=50)
    fitter.set_param_hint('vsini', value=vsini.value, vary=True, min=0.0, max=500.0)
    fitter.set_param_hint('temperature', value=temperature, min=T_min, max=T_max, vary=True)
    fitter.set_param_hint('logg', value=logg, min=logg_min, max=logg_max, vary=True)
    fitter.set_param_hint('metal', value=metal, min=metal_min, max=metal_max, vary=True)
    fitter.set_param_hint('alpha', value=alpha, min=alpha_min, max=alpha_max, vary=True)
    params = fitter.make_params()

    #Perform the fit
    optdict = {"epsfcn": 1e-2}
    result = fitter.fit(orders, fit_kws=optdict, params=params)
    print(result.fit_report())
    for i, order in enumerate(orders):
        m = result.best_fit[i]
        ratio = order.y / m
        order.cont = FittingUtilities.Continuum(order.x, ratio, lowreject=2, highreject=2, fitorder=5)
        plt.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
        plt.plot(order.x, result.best_fit[i], 'r-', alpha=0.5)
    plt.show()


    # Now, re-do the fit using emcee to get realistic error bars


    #Now, fit the fundamental parameters
    params['rv'].vary = False
    params['vsini'].vary = True
    params['temperature'].vary = True
    params['logg'].vary = True
    params['metal'].vary = True
    params['alpha'].vary = True
    result2 = fitter.fit(orders, params=params)
    print(result2.fit_report())
    for i, order in enumerate(orders):
        m = result.best_fit[i]
        ratio = order.y / m
        order.cont = FittingUtilities.Continuum(order.x, ratio, lowreject=2, highreject=2, fitorder=5)
        plt.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
        plt.plot(order.x, result.best_fit[i], 'r-', alpha=0.5)
    plt.show()

    """
    #Perform the fit - Fit the RV and vsini to an initial guess model first
    params = Parameters()
    params.add('rv', value=rv.value, min=-50, max=50)
    params.add('vsini', value=vsini.value, vary=True, min=0.0, max=500.0)
    params.add('temperature', value=temperature, min=T_min, max=T_max, vary=False)
    params.add('logg', value=logg, min=logg_min, max=logg_max, vary=False)
    params.add('metal', value=metal, min=metal_min, max=metal_max, vary=False)
    params.add('alpha', value=alpha, min=alpha_min, max=alpha_max, vary=False)
    #result = minimize(ErrorFunction, params, args=(data, model_fcn))
    result = minimize(ErrorFunction2, params, args=(orders, mg))
    report_fit(params)

    #plt.plot(orders[7].x, orders[7].cont, 'r--')
    #plt.show()

    #Now, fit the fundamental parameters
    params['rv'].vary = False
    params['vsini'].vary = True
    params['temperature'].vary = True
    params['logg'].vary = True
    params['metal'].vary = True
    params['alpha'].vary = True
    fitter.make_params()
    result = minimize(ErrorFunction, params, args=(orders, mg))

    report_fit(params)
    """

    #Plot
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212, sharex=ax1)
    model_orders = MakeModel(params, orders, mg)
    for order, model in zip(orders, model_orders):
        ratio = order.y / model
        order.cont = FittingUtilities.Continuum(order.x, ratio, fitorder=5, lowreject=2, highreject=2)
        ax1.plot(order.x, order.y / order.cont, 'k-')
        ax1.plot(order.x, model, 'r-')
        ax2.plot(order.x, ratio / order.cont)
    plt.show()

