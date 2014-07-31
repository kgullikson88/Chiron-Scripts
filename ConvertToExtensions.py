import sys
import os

import FittingUtilities
from astropy.io import fits as pyfits
import numpy as np
import pylab

import FitsUtils
import HelperFunctions


minsnr = 10.0

if __name__ == "__main__":
    fileList = []
    blazecorrect = True
    for arg in sys.argv[1:]:
        if "noblaze" in arg:
            blazecorrect = False
        else:
            fileList.append(arg)
    for fname in fileList:
        outfilename = "%s-0.fits" % (fname.split(".fits")[0])
        header = pyfits.getheader(fname)
        if header["IMAGETYP"] != "object":
            print "Image type is %s. Skipping" % header["IMAGETYP"]
            continue
        orders = HelperFunctions.ReadFits(fname, errors=2)
        orders = orders[::-1]  # Reverse order so the bluest order is first
        if blazecorrect:
            blazefile = [f for f in os.listdir("./") if "eblaze" in f]
            if len(blazefile) != 1:
                blazefile = raw_input("Enter blaze filename: ")
            else:
                blazefile = blazefile[0]

            try:
                blaze = HelperFunctions.ReadFits(blazefile, errors=2)
                blaze = blaze[::-1]
            except IOError:
                print "Error! blaze file %s does not exist!" % blazefile
                print "Not converting file %s" % fname
                continue
        column_list = []
        maxsnr = 0
        for i, order in enumerate(orders):
            # Blaze correction
            if blazecorrect:
                order.y /= blaze[i].y
                order.err /= blaze[i].y

            zeros = np.where(order.y < 0)[0]
            order.y[zeros] = 0.0

            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=4)
            columns = columns = {"wavelength": order.x,
                                 "flux": order.y,
                                 "continuum": order.cont,
                                 "error": order.err}
            column_list.append(columns)
            snr = 1.0 / np.std(order.y / order.cont)
            if snr > maxsnr:
                maxsnr = snr

        print "Maximum S/N ratio = %g" % maxsnr
        if maxsnr > minsnr:
            HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")
        else:
            print "\nSkipping file %s (Very low S/N ratio)\n" % (fname)
      
