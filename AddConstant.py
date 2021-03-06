import sys

import numpy as np

import FitsUtils
import HelperFunctions


"""
  This script will read in a CHIRON reduced spectrum find the lowest point and,
  if that point has a flux < 0, will add a constant such that the lowest flux is 0
"""

if __name__ == "__main__":
    fileList = []
    outfilename = None
    count = False
    for arg in sys.argv[1:]:
        if "-out" in arg:
            outfilename = arg.split("=")[-1]
        else:
            fileList.append(arg)

    if len(fileList) > 1:
        count = True
        counter = 1

    for fname in fileList:
        orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum",
                                        errors="error")
        lowest = 9e30
        for order in orders:
            if np.min(order.y) < lowest:
                lowest = np.min(order.y)

        if lowest < 0:
            column_list = []
            for order in orders:
                order.y -= lowest
                order.cont -= lowest
                order.err -= lowest
                column = {"wavelength": order.x,
                          "flux": order.y,
                          "continuum": order.cont,
                          "error": order.err}
                column_list.append(column)
            if outfilename == None:
                if "-" in fname:
                    idx = int(fname.split("-")[-1].split(".fits")[0])
                    outfilename = "%s-%i.fits" % (fname.split("-")[0], idx + 1)
                else:
                    outfilename = "%s-0.fits" % (fname.split(".fits")[0])
            elif count:
                outfilename = "%s-%i.fits" % (outfilename.split(".fits")[0], counter)
                counter += 1
            print "Outputting new file to %s" % outfilename
            HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
        else:
            print "Lowest point is >0. Not making a new file!"
