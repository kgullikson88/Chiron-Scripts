"""
Sensitivity analysis, using the new search method.
"""
import sys

import Sensitivity
import StarData

import Search_slow

if __name__ == '__main__':
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    """ Use for testing purposes before long runs!
    for fname in fileList:
        header = fits.getheader(fname)
        star = header['object']
        Sensitivity.get_companions(star)
    """

    badregions = Search_slow.badregions
    interp_regions = Search_slow.interp_regions
    trimsize = Search_slow.trimsize
    prim_vsini = StarData.get_vsini(fileList)

    Sensitivity.Analyze(fileList, prim_vsini,
                        hdf5_file='/media/ExtraSpace/PhoenixGrid/CHIRON_Grid.hdf5',
                        extensions=True,
                        resolution=None,
                        trimsize=trimsize,
                        badregions=badregions, interp_regions=interp_regions,
                        metal_values=(0.0,),
                        vsini_values=(5,),
                        Tvalues=range(3000, 6100, 100),
                        debug=False,
                        addmode='simple',
                        output_mode='hdf5')
