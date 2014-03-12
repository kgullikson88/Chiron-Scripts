import numpy
#import FitsUtils
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
from collections import Counter
from sklearn.gaussian_process import GaussianProcess

plot = False

def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10):
  denoised = HelperFunctions.Denoise(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters)
  denoised.y /= denoised.y.max()
  return denoised




def GPSmooth(data, low=0.1, high=10, debug=False):
  """
  This will smooth the data using Gaussian processes. It will find the best
  smoothing parameter via cross-validation to be between the low and high.

  The low and high keywords are reasonable bounds for  A and B stars with 
  vsini > 100 km/s.
  """
    
  gp = GaussianProcess(corr='squared_exponential',
                       theta0 = numpy.sqrt(low*high),
                       thetaL = low,
                       thetaU = high,
                       normalize = False,
                       nugget = (data.err / data.y)**2,
                       random_start=1)
  try:
    gp.fit(data.x[:,None], data.y)
  except ValueError:
    #On some orders with large telluric residuals, this will fail.
    # Just fall back to the old smoothing method in that case.
    return SmoothData(data), 91
  smoothed = data.copy()
  if debug:
    print "\tSmoothing parameter theta = ", gp.theta_
  smoothed.y, smoothed.err = gp.predict(data.x[:,None], eval_MSE=True)
  return smoothed, gp.theta_[0][0]


if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)
  if len(fileList) == 0:
    fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    column_list = []
    header_list = []
    for i, order in enumerate(orders):
      print "Smoothing order %i/%i" %(i+1, len(orders))
      #Fix errors
      order.err[order.err > 1e8] = numpy.sqrt(order.y[order.err > 1e8])

      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = FittingUtilities.RebinData(order, xgrid)
      
      #denoised = SmoothData(order, 101, 5, 2, 2, 10)
      denoised, theta = GPSmooth(order)

      column = {"wavelength": denoised.x,
                "flux": order.y / denoised.y,
                "continuum": denoised.cont,
                "error": denoised.err}
      header_list.append(("Smoother", theta, "Smoothing Parameter"))
      column_list.append(column)
      if plot:
        plt.figure(1)
        plt.plot(order.x, order.y/order.y.mean())
        plt.plot(denoised.x, denoised.y/denoised.y.mean())
        plt.figure(2)
        plt.plot(order.x, order.y/denoised.y)
        plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
        plt.show()
    if plot:
      plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new', headers_info=header_list)
