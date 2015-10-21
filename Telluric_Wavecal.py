"""
Adjust the wavelength calibration of CHIRON spectra
to match the telluric lines
"""

from astropy.io import fits
import HelperFunctions
import lmfit
import os
from functools import partial
import logging
from collections import defaultdict
import numpy as np
import telfit
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy import constants, units as u 
import FittingUtilities

# Make the speed of light a constant
C_LIGHT = constants.c.cgs.to(u.km/u.s).value

class VelocityFitter(object):
    def __init__(self, filename, tell_orders=[690., 700., 715., 725., 735.]):
        """
        Fit the velocity shift, with uncertainties, to make the spectrum line up with the telluric lines

        Parameters:
        ============
         - filename:         string
                             The filename of a combined spectrum. It should have header keywords 
                             such as 'FILE1', 'FILE2', etc that give the filenames of the original
                             spectra that still have telluric lines in them.

         - tell_orders:      iterable object
                             Should contain either the order numbers to use for fitting, or wavelengths
                             that map to a unique order. The default is the O2 and H2O band from 690-740 nm.

        """
        # First thing: Find all of the original files.
        original_files = self._find_original_files(filename)

        # Compile all of the relevant orders for each file
        self.data = defaultdict(list)
        min_wave = np.inf
        max_wave = -np.inf
        for fname in original_files:
            orders = HelperFunctions.ReadExtensionFits(fname)
            basename = os.path.basename(fname)
            for n in tell_orders:
                if isinstance(n, int) and n < len(orders):
                    self.data[basename].append(orders[n])
                else:
                    ordernum = HelperFunctions.FindOrderNums(orders, [n])[0]
                    self.data[basename].append(orders[ordernum])
                if self.data[basename][-1].x[0] < min_wave:
                    min_wave = self.data[basename][-1].x[0]
                if self.data[basename][-1].x[-1] > max_wave:
                    max_wave = self.data[basename][-1].x[-1]

        # Generate a generic telluric model
        modeler = telfit.Modeler()
        lowfreq = 1e7/(max_wave+2)
        highfreq = 1e7/(min_wave-2)
        original_model = modeler.MakeModel(lowfreq=lowfreq, highfreq=highfreq)

        # Reduce the resolution to match CHIRON
        new_x = np.logspace(np.log(original_model.x[0]), np.log(original_model.x[-1]), original_model.size(), base=np.e)
        full_model = FittingUtilities.RebinData(original_model, new_x)
        full_model = FittingUtilities.ReduceResolutionFFT(full_model, 110000)

        # Interpolate the telluric model
        self.telluric_model = spline(full_model.x, full_model.y)

        # Make a function to adjust the telluric model to better match the data.
        self.new_telluric = lambda tell, alpha, x: tell(x)**alpha


    def lnlike(self, pars, data_orders, telluric_model):
        # Extract parameters
        parvals = pars.valuesdict()
        rv, alpha = pars['rv'].value, pars['alpha'].value

        # Adjust telluric model
        telluric = partial(self.new_telluric, telluric_model, alpha)

        # Calculate log-likelihood
        ll = []
        for order in data_orders:
            x = order.x*(1+rv/C_LIGHT)
            y = order.y / order.cont
            e = order.err / order.cont
            #ll.append(-0.5 * (y - telluric(x))**2 / e**2 + np.log(2*np.pi*(e**2)))
            ll.append((y-telluric(x))/e)
        ll = np.hstack(ll)
        logging.debug('RV: {}\nalpha = {}\nlogL = {}\n'.format(rv, alpha, ll.sum() ))
        return ll

    def _errfcn(self, pars, data_orders, telluric_model):
        return -self.lnlike(pars, data_orders, telluric_model)


    def _fit_one_file(self, data_orders, rv_guess=0.0, alpha_guess=1.0):
        # Create parameters
        params = lmfit.Parameters()
        params.add('rv', value=rv_guess, min=-5, max=5)
        params.add('alpha', value=alpha_guess, min=0, max=3)

        # Perform the fit
        result = lmfit.minimize(self._errfcn, params, args=(data_orders, self.telluric_model))

        logger = logging.getLogger()
        if logger.level <= logging.DEBUG:
            lmfit.report_fit(params)

        return params


    def fit(self, rv_guess=0.0, alpha_guess=1.0):
        """
        Perform the fit on all original files, and then combine into a mean rv.
        """
        rv = np.zeros(len(self.data))
        rv_err = np.zeros(len(self.data))
        for i, key in enumerate(self.data.keys()):
            logging.debug('Fitting file {}'.format(os.path.basename(key)))
            best_pars = self._fit_one_file(self.data[key], rv_guess=rv_guess, alpha_guess=alpha_guess)
            rv[i] = best_pars['rv'].value
            rv_err[i] = best_pars['rv'].stderr

        # Make sure no errors are impossibly small.
        rv_err[np.abs(rv_err) < 1e-6] = 10*rv[np.abs(rv_err) < 1e-6]

        weights = 1.0 / rv_err**2
        best_rv = np.average(rv, weights=weights)
        rv_variance = np.average((rv-best_rv)**2, weights=weights)

        # Get the un-biased variance using the formula from https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
        V1 = np.sum(weights)
        V2 = np.sum(weights**2)
        rv_variance /= (1-V2/V1**2)

        logging.info('RV = {:.2f} +/- {:.3f}'.format(best_rv, np.sqrt(rv_variance)))
        return best_rv, np.sqrt(rv_variance)




    def _find_original_files(self, filename):
        """
        Find the original files that went into the co-added file given. We search first in the
        same directory as the given file, then in the 'backup' subdirectory, and finally ask 
        the user if they are not found.
        """
        header = fits.getheader(filename)
        names = ['{}.fits'.format(header[k].split('_telluric_corrected')[0]) for k in header.keys() if k.startswith('FILE')]
        basedir = os.path.dirname(filename)
        if len(basedir) > 0 and not basedir.endswith('/'):
            basedir += '/'

        # Check in the same directory as the given file
        if all([os.path.exists('{}{}'.format(basedir, f)) for f in names]):
            return ['{}{}'.format(basedir, f) for f in names]

        # Check the backup subdirectory
        if all([os.path.exists('{}backup/{}'.format(basedir, f)) for f in names]):
            return ['{}backup/{}'.format(basedir, f) for f in names]

        # Ask the user
        done = False
        print('Where are the following files located?')
        for n in names:
            print('\t{}'.format(n))
        while not done:
            basedir = raw_input('Enter directory where the files can be found: ')
            if all([os.path.exists('{}{}'.format(basedir, f)) for f in names]):
                done = True
        return ['{}{}'.format(basedir, f) for f in names]




