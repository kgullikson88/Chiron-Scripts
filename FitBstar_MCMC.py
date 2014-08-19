__author__ = 'kgulliks'
import os
import warnings
import sys

from scipy.interpolate import InterpolatedUnivariateSpline as spline, LinearNDInterpolator, NearestNDInterpolator
import numpy as np
import DataStructures
from astropy import units as u, constants
import HelperFunctions
import Broaden
import emcee


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

        # Make the interpolator instance
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

        # Test to make sure the result is valid. If the requested point is
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

def MakeModel(pars, data, model_getter):
    vsini = pars[0]
    rv = pars[1]
    T = pars[2]
    logg = pars[3]
    metal = pars[4]
    alpha = pars[5]
    lnf = pars[6]

    c = constants.c.cgs.to(u.km / u.s).value

    # Get the model from the ModelGetter instance
    model = model_getter(T, logg, metal, alpha)

    #First, broaden the model
    broadened = Broaden.RotBroad(model, vsini, linear=True)
    modelfcn = spline(broadened.x, broadened.y)
    model_orders = []
    for order in data:
        model_orders.append(modelfcn(order.x * (1 - rv / c)))
    return model_orders


def lnprior(theta, bounds):
    vsini = theta[0]
    rv = theta[1]
    T = theta[2]
    logg = theta[3]
    metal = theta[4]
    alpha = theta[5]
    lnf = theta[6]
    if (bounds[0]['lower'] < vsini < bounds[0]['upper'] and
                    bounds[1]['lower'] < rv < bounds[1]['upper'] and
                    bounds[2]['lower'] < T < bounds[2]['upper'] and
                    bounds[3]['lower'] < logg < bounds[3]['upper'] and
                    bounds[4]['lower'] < metal < bounds[4]['upper'] and
                    bounds[5]['lower'] < alpha < bounds[5]['upper'] and
                lnf < 1.0):
        return 0.0

    return -np.inf


def lnlike(theta, data, model_getter):
    model_orders = MakeModel(theta, data, model_getter)
    lnf = theta[6]
    loglikelihood = 0.0
    N = 0.0
    for i, (o, m) in enumerate(zip(data, model_orders)):
        ratio = o.y / m
        # cont = FittingUtilities.Continuum(o.x, ratio, fitorder=5, lowreject=2, highreject=2)
        mean = o.x.mean()
        o.cont = np.poly1d(np.polyfit(o.x - mean, ratio, 5))(o.x - mean)

        inv_sigma2 = 1.0 / (o.err ** 2 + m ** 2 * np.exp(2 * lnf))
        loglikelihood += -0.5*(np.sum((o.y - o.cont * m) ** 2 * inv_sigma2 - np.log(inv_sigma2)))
        N += o.size()

    return loglikelihood


def lnprob(theta, data, bounds, model_getter):
    lp = lnprior(theta, bounds)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, data, model_getter)


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
    vsini = 0.0 * u.km / u.s
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





    # Find the models
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


    #Make guess values for each of the values from the bounds
    T = (T_min + T_max) / 2.0
    logg = (logg_min + logg_max) / 2.0
    metal = (metal_min + metal_max) / 2.0
    alpha = (alpha_min + alpha_max) / 2.0
    vsini = 100.0
    rv = 0.0
    lnf = 0.5
    bounds = [dict(lower=0, upper=500),
              dict(lower=-50, upper=50),
              dict(lower=T_min, upper=T_max),
              dict(lower=logg_min, upper=logg_max),
              dict(lower=metal_min, upper=metal_max),
              dict(lower=alpha_min, upper=alpha_max)]

    pars = np.array([vsini, rv, T, logg, metal, alpha, lnf])

    #Set up the MCMC sampler
    ndim, nwalkers = 3, 100
    pos = [pars + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(orders, bounds, mg))

    #Run the MCMC sampler
    sampler.run_mcmc(pos, 1000)

    #Get the samples after burn-in
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    #Save the output
    import triangle

    fig = triangle.corner(samples, labels=["$v\sin{i}$", "$RV$", "Temperature", "$\log{g}$", "[Fe/H]", "[$\alpha$/Fe]",
                                           "$\ln\,f$"])
    fig.savefig("triangle.pdf")
