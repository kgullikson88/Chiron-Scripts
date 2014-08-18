import os
import warnings
import sys
import FittingUtilities

from scipy.interpolate import InterpolatedUnivariateSpline as spline, LinearNDInterpolator, NearestNDInterpolator
import numpy as np
import matplotlib.pyplot as plt
import DataStructures
from astropy import units as u, constants
import HelperFunctions
import Broaden

from lmfit import minimize, Parameters, report_fit


class ModelGetter():
    def __init__(self, modeldir, rebin=True, T_min=7000, T_max=9000, logg_min=3.5, logg_max=4.5, metal_min=-0.5,
                 metal_max=0.5, alpha_min=0.0, alpha_max=0.4, wavemin=0, wavemax=np.inf):
        """
        This class will read in a directory with Kurucz models

        The associated methods can be used to interpolate a model at any
        temperature, gravity, metallicity, and [alpha/Fe] value that
        falls within the grid

        modeldir: The directory where the models are stored
        rebin: If True, it will rebin the models to a constant x-spacing
        other args: The minimum and maximum values for the parameters to search.
                    You need to keep this as small as possible to avoid memory issues!
                    The whole grid would take about 36 GB of RAM!
        """

        # First, read in the grid
        Tvals = []
        loggvals = []
        metalvals = []
        alphavals = []
        spectra = []
        firstkeeper = True
        modelfiles = [f for f in os.listdir(modeldir) if f.startswith("t") and f.endswith(".dat.bin.asc")]
        for i, fname in enumerate(modelfiles):
            T = float(fname[1:6])
            logg = float(fname[8:12])
            metal = float(fname[14:16]) / 10.0
            alpha = float(fname[18:20]) / 10.0
            if fname[13] == "m":
                metal *= -1
            if fname[17] == "m":
                alpha *= -1

            # Read in and save file if it falls in the correct parameter range
            if (T_min <= T <= T_max and
                            logg_min <= logg <= logg_max and
                            metal_min <= metal <= metal_max and
                            alpha_min <= alpha <= alpha_max):

                print "Reading in file {:s}".format(fname)
                x, y = np.loadtxt("{:s}/{:s}".format(modeldir, fname), usecols=(0, 3), unpack=True)
                x *= u.angstrom.to(u.nm)

                left = np.searchsorted(x, wavemin)
                right = np.searchsorted(x, wavemax)
                x = x[left:right]
                y = y[left:right]

                if rebin:
                    xgrid = np.linspace(x[0], x[-1], x.size)
                    fcn = spline(x, y)
                    x = xgrid
                    y = fcn(xgrid)

                if firstkeeper:
                    self.xaxis = x
                    firstkeeper = False
                elif np.max(np.abs(self.xaxis - x) > 1e-4):
                    warnings.warn("x-axis for file {:s} is different from the master one! Not saving!".format(fname))
                    continue

                Tvals.append(T)
                loggvals.append(logg)
                metalvals.append(metal)
                alphavals.append(alpha)
                spectra.append(y)

        # Save the grid as a class variable
        self.grid = np.array((Tvals, loggvals, metalvals, alphavals)).T
        self.spectra = np.array(spectra)

        #Make the interpolator instance
        self.interpolator = LinearNDInterpolator(self.grid, self.spectra, rescale=True)
        self.NN_interpolator = NearestNDInterpolator(self.grid, self.spectra, rescale=True)


    def __call__(self, T, logg, metal, alpha, return_xypoint=True):
        """
        Given parameters, return an interpolated spectrum

        If return_xypoint is False, then it will only return
          a numpy.ndarray with the spectrum

        Before interpolating, we will do some error checking to make
        sure the requested values fall within the grid
        """

        # Get the minimum and maximum values in the grid
        T_min = min(self.grid[:, 0])
        T_max = max(self.grid[:, 0])
        logg_min = min(self.grid[:, 1])
        logg_max = max(self.grid[:, 1])
        metal_min = min(self.grid[:, 2])
        metal_max = max(self.grid[:, 2])
        alpha_min = min(self.grid[:, 3])
        alpha_max = max(self.grid[:, 3])

        # Check to make sure the requested values fall within the grid
        if (T_min <= T <= T_max and
                        logg_min <= logg <= logg_max and
                        metal_min <= metal <= metal_max and
                        alpha_min <= alpha <= alpha_max):

            y = self.interpolator((T, logg, metal, alpha))
        else:
            warnings.warn("The requested parameters fall outside the model grid. Results may be unreliable!")
            print T, T_min, T_max
            print logg, logg_min, logg_max
            print metal, metal_min, metal_max
            print alpha, alpha_min, alpha_max
            y = self.interpolator((T, logg, metal, alpha))

        #Test to make sure the result is valid. If the requested point is
        #outside the Delaunay triangulation, it will return NaN's
        if np.any(np.isnan(y)):
            warnings.warn("Found NaNs in the interpolated spectrum! Falling back to Nearest Neighbor")
            y = self.NN_interpolator((T, logg, metal, alpha))

        #Return the appropriate object
        if return_xypoint:
            return DataStructures.xypoint(x=self.xaxis, y=y)
        else:
            return y


# ########################################################################
# ########################################################################


def ErrorFunction2(params, data, model_getter):
    model_orders = MakeModel(params, data, model_getter)
    loglikelihood = []
    N = 0.0
    for i, (o, m) in enumerate(zip(data, model_orders)):
        ratio = o.y / m
        cont = FittingUtilities.Continuum(o.x, ratio, fitorder=5, lowreject=2, highreject=2)
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

    #First, broaden the model
    broadened = Broaden.RotBroad(model, vsini, linear=True)
    modelfcn = spline(broadened.x, broadened.y)
    model_orders = []
    for order in data:
        model_orders.append(modelfcn(order.x * (1 - rv / c)))
    return model_orders


if __name__ == "__main__":
    # Define some constants
    c = constants.c

    # Set up default values, and then read in command line arguments
    T_min = 9000
    T_max = 9250
    metal_min = -0.8
    metal_max = 0.5
    logg_min = 3.0
    logg_max = 3.5
    alpha_min = 0.0
    alpha_max = 0.4
    modeldir = "models/"
    rv = 0.0 * u.km / u.s
    vsini = 0.0 * u.km / u.x
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
                modeldir = modeldir + "/"
        elif "rv" in arg.lower():
            rv = float(arg.partition("=")[-1]) * u.km / u.s
        elif "vsini" in arg.lower():
            vsini = float(arg.partition("=")[-1]) * u.km / u.s
        elif "resolution" in arg.lower():
            R = float(arg.partition("=")[-1])
        else:
            filename = arg





    #Find the models
    # model_list = sorted(["{:s}{:s}".format(modeldir, f) for f in os.listdir(modeldir)])
    #model_list = sorted(["models/%s" % f for f in os.listdir("models")])

    #Read the data
    orders = HelperFunctions.ReadExtensionFits(filename)

    #plt.plot(orders[7].x, orders[7].y)
    #plt.plot(orders[7].x, orders[7].cont)

    #Make an instance of the model getter
    mg = ModelGetter(modeldir,
                     T_min=T_min,
                     T_max=T_max,
                     logg_min=logg_min,
                     logg_max=logg_max,
                     metal_min=metal_min,
                     metal_max=metal_max,
                     alpha_min=alpha_min,
                     alpha_max=alpha_max,
                     wavemin=orders[0].x[0] - 1)
    # mg = ModelGetter("models", T_min=9000, T_max=9250, metal_min=-1.0, metal_max=1.0, wavemin=orders[0].x[0] - 1.0)

    #Make guess values for each of the values from the bounds
    temperature = (T_min + T_max) / 2.0
    logg = (logg_min + logg_max) / 2.0
    metal = (metal_min + metal_max) / 2.0
    alpha = (alpha_min + alpha_max) / 2.0

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
    result = minimize(ErrorFunction, params, args=(orders, mg))

    report_fit(params)

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

