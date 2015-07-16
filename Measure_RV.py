import Fitters
import HelperFunctions
import sys
import StarData
import SpectralTypeRelations
import Mamajek_Table
import os
from astropy.io import fits
import numpy as np

home = os.environ['HOME']
sys.path.append('{}/School/Research/CHIRON_data/Chiron-Scripts/'.format(home))
from Search_slow import hdf5_filename

MT = Mamajek_Table.MamajekTable()
MS = SpectralTypeRelations.MainSequence()
sptnum2teff = MT.get_interpolator('SpTNum', 'Teff') 


def make_fitter(fname, logg=4.0, feh=0.0):
    orders = HelperFunctions.ReadExtensionFits(fname)
    header = fits.getheader(fname)
    starname = header['OBJECT']

    SpT = StarData.GetData(starname).spectype
    Teff = sptnum2teff(MS.SpT_To_Number(SpT))

    # Find the nearest grid temperature
    hdf5_int = Fitters.StellarModel.HDF5Interface(hdf5_filename)
    grid_teffs = np.array([d['temp'] for d in hdf5_int.list_grid_points if d['Z'] == feh and d['logg'] == logg])
    idx = np.argmin((grid_teffs - Teff)**2)
    Teff = grid_teffs[idx]
    print(Teff)


    print(hdf5_filename)
    fitter = Fitters.RVFitter_MultiNest(orders[8:20], model_library=hdf5_filename, T=Teff, logg=logg, feh=feh)

    return fitter



if __name__ == '__main__':
    file_list = []
    for arg in sys.argv[1:]:
        if 1:
            file_list.append(arg)

    for fname in file_list:
        fitter = make_fitter(fname)

