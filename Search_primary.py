import sys

import GenericSearch
import StarData

# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = []
interp_regions = []
trimsize = 10

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
    hdf5_filename = '/Volumes/DATADRIVE/Kurucz_Grid/CHIRON_grid_air.hdf5'
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
    hdf5_filename = '/media/ExtraSpace/Kurucz_FullGrid/CHIRON_grid_air.hdf5'
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    # Get the primary star vsini values
    prim_vsini = [None for _ in fileList]

    # Use this one for the real data search
    GenericSearch.slow_companion_search(fileList, prim_vsini,
                                        hdf5_file=hdf5_filename,
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        metal_values=(0.0),
                                        vsini_values=(10.0, 40, 80, 150),
                                        Tvalues=range(8000, 20000, 1000),
                                        observatory='CTIO',
                                        debug=False,
                                        reject_outliers=False,
                                        vbary_correct=True,
                                        addmode='all',
                                        output_mode='hdf5',
                                        output_file='CCF_primary.hdf5')