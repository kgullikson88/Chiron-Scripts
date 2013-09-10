import numpy
import FitsUtils
import FittingUtilities
import matplotlib.pyplot as plt
import sys
from astropy import units
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp


def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10):
  denoised = FittingUtilities.Denoise3(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters)
  denoised.y /= denoised.y.max()
  return denoised
  


if __name__ == "__main__":
  for fname in sys.argv[1:]:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    column_list = []
    for order in orders:
      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = MakeModel.RebinData(order, xgrid)
      
      denoised = SmoothData(order, 91, 5, 2, 2, 10)
      #order2 = order.copy()
      #denoised = FittingUtilities.Denoise3(order2) #, snr=400.0, reduction_factor=0.15)
      #denoised.y = FittingUtilities.Iterative_SV(denoised.y, 91, 5, lowreject=2, highreject=2, numiters=10)

      column = {"wavelength": denoised.x,
                "flux": (order.y - denoised.y) + numpy.median(order.cont),
                "continuum": denoised.cont,
                "error": denoised.err}
      column_list.append(column)
      plt.figure(1)
      plt.plot(order.x, order.y)
      plt.plot(denoised.x, denoised.y)
      plt.figure(2)
      plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
    plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    FitsUtils.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
