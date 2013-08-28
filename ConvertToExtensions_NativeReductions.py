import FitsUtils
import FittingUtilities
import pyfits
import sys
import os
import numpy
import matplotlib.pyplot as plt
from astropy import units


if __name__ == "__main__":
  fileList = []
  blazecorrect = True
  for arg in sys.argv[1:]:
    if "blaze" in arg:
      blazecorrect = False
    else:
      fileList.append(arg)
  for fname in fileList:
    outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    header = pyfits.getheader(fname)
    orders = FitsUtils.MakeXYpoints(fname)
    #orders = orders[::-1]    #Reverse order so the bluest order is first

    column_list = []
    for i, order in enumerate(orders[:-1]):
      #Convert to nanometers
      order.x *= units.Angstrom.to(units.nm)
      plt.plot(order.x, order.y)
      
      #Blaze correction
      if blazecorrect:
        blaze = FittingUtilities.Continuum(order.x-order.x.mean(), order.y, fitorder=9, lowreject=1.5, highreject=5)
        #plt.plot(order.x, blaze)
        #plt.figure(1)
        #plt.plot(order.x, order.y)
        #plt.plot(order.x, blaze)
        blaze /= blaze.max()
        order.y /= blaze
        order.err /= blaze
        
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=5)
      #plt.figure(2)
      #plt.plot(order.x, order.y/order.cont)
      columns = columns = {"wavelength": order.x,
                           "flux": order.y,
                           "continuum": order.cont,
                           "error": order.err}
      column_list.append(columns)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux (counts)")
    plt.show()
    FitsUtils.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')


      
