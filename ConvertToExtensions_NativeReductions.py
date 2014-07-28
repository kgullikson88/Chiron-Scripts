import FitsUtils_NativeReductions as FitsUtils
import FittingUtilities
from astropy.io import fits as pyfits
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import HelperFunctions


def Convert(fileList, blazecorrect=True, plotresult=False):

  #Read blaze orders
  if blazecorrect:
    blazefile = [f for f in os.listdir("./") if "slicerflat" in f][0]
    hdulist = pyfits.open(blazefile)
    data = hdulist[0].data
    blaze_orders = [data[2][i][::-1] / data[2][i].mean() for i in range(data.shape[1])][::-1]
  
  for fname in fileList:
    outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    header = pyfits.getheader(fname)
    if header["IMAGETYP"] != "object":
      print "Image type is %s. Skipping" %header["IMAGETYP"]
      continue
    orders = FitsUtils.MakeXYpoints(fname)
    #orders = orders[::-1]    #Reverse order so the bluest order is first

    column_list = []
    for i, order in enumerate(orders[:-1]):      
      #Convert to nanometers
      order.x *= units.Angstrom.to(units.nm)
      
      if plotresult:
        plt.figure(1)
        plt.plot(order.x, blaze_orders[i] * np.max(order.y) / np.max(blaze_orders[i]))
        plt.plot(order.x, order.y)
      #Blaze correction
      if blazecorrect:
        blaze_orders[i][blaze_orders[i] < 1e-5] = 1.0
        order.y /= blaze_orders[i]

      if plotresult:
        plt.figure(2)
        plt.plot(order.x, order.y)
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=5)
      #plt.figure(2)
      #plt.plot(order.x, order.y/order.cont)
      columns = {"wavelength": order.x,
                           "flux": order.y,
                           "continuum": order.cont,
                           "error": order.err}
      column_list.append(columns)
    if plotresult:
      plt.xlabel("Wavelength (nm)")
      plt.ylabel("Flux (counts)")
      plt.show()
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')


      

if __name__ == "__main__":
  fileList = []
  blazecorrect = True
  plotresult = False
  for arg in sys.argv[1:]:
    if "blaze" in arg:
      blazecorrect = False
    elif "-p" in arg:
      plotresult = True
    else:
      fileList.append(arg)

  Convert(fileList, blazecorrect, plotresult)
