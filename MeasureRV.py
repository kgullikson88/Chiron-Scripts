"""
Measure the radial velocity and vsini of flattened spectra
"""
from __future__ import print_function, division, absolute_import

import logging
import os
import glob

import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import seaborn as sns

import HelperFunctions
import Fitters

# Set up plotting
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5)

# Set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Get the HDF5 filename. Might want to change this eventually.
#HDF5_FILENAME = '/Volumes/DATADRIVE/Kurucz_Grid/TS23_grid_full.hdf5'
HDF5_FILENAME = '/Users/kevingullikson/StellarLibrary/Kurucz_Grid/CHIRON_grid_air.hdf5'
PAR_LOGFILE = 'Flatten.log'


def fit(filename, model_library, teff, logg, feh=0.0, output_basename='RVFitter'):
    # Read in the (assumed flattened) spectra
    all_orders = HelperFunctions.ReadExtensionFits(filename)
    orders = [o.copy() for o in all_orders if o.x[0] < 475 or o.x[-1] > 495]

    # Set up the fitter
    fitter = Fitters.RVFitter(orders, model_library=model_library,
                              T=teff, logg=logg, feh=feh)
    header = fits.getheader(filename)
    starname = header['OBJECT']
    date = header['DATE-OBS'].split('T')[0]
    stardata_str = '{}_{}-'.format(starname.replace(' ', ''), date.replace('-', ''))
    basename = os.path.join(output_basename, stardata_str)

    # Fit
    fitter.fit(backend='multinest', n_live_points=1000, basename=basename, overwrite=False, init_MPI=False)

    # Make a triangle plot and save it
    fitter.triangle()
    plt.savefig('{}triangle.pdf'.format(basename))

    return fitter

if __name__ == '__main__':
    file_list = glob.glob('201*/*renormalized.fits')
    fitted_df = pd.read_csv(PAR_LOGFILE, header=None, names=['fname', 'star', 'date', 'teff', 'logg', 'rv'])
    print(fitted_df.tail())
    
    for filename in file_list:
        logging.info('Fitting RV for {}'.format(filename))

        # Find this filename in the fitted dataframe (generated while flattening the spectra)
        header = fits.getheader(filename)
        starname = header['OBJECT']
        date = header['DATE-OBS']
        print(starname, date)
        subset = fitted_df.loc[(fitted_df.star==starname) & (fitted_df.date==date)]
        print(subset)
        teff = float(subset.teff)
        logg = float(subset.logg)
        logging.info('Teff = {}\nlogg = {}'.format(teff, logg))

        fitter = fit(filename, HDF5_FILENAME, teff=teff, logg=logg, output_basename='RVFitter_nobalmer_veiling')
