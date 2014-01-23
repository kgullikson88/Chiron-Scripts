import numpy
import FitsUtils
import FittingUtilities
import HelperFunctions
import matplotlib.pyplot as plt
import sys
import os
from astropy import units
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import MakeModel
import HelperFunctions

plot = False

def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10):
  denoised = HelperFunctions.Denoise(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters)
  denoised.y /= denoised.y.max()
  return denoised
  


if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)
  if len(fileList) == 0:
    fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    column_list = []
    for order in orders:
      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = FittingUtilities.RebinData(order, xgrid)
      
      denoised = SmoothData(order, 101, 5, 2, 2, 10)

      column = {"wavelength": denoised.x,
                "flux": order.y / denoised.y,
                "continuum": denoised.cont,
                "error": denoised.err}
      column_list.append(column)
      if plot:
        plt.figure(1)
        plt.plot(order.x, order.y/order.y.mean())
        plt.plot(denoised.x, denoised.y/denoised.y.mean())
        plt.figure(2)
        plt.plot(order.x, order.y/denoised.y)
        plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
    if plot:
      plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
