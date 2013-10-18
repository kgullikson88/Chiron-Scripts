import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from scipy.optimize import leastsq, brute
from scipy import mat
from scipy.linalg import svd, diagsvd
import sys
import os
import FitsUtils
from collections import defaultdict
import DataStructures
from astropy import units, constants
import MakeModel
import RotBroad
import FittingUtilities
import FindContinuum



def Smooth(fname, window_size=91, numiters=100, lowreject=3, highreject=3, smoothorder=3):
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")
  column_list = []
  for i, order in enumerate(orders):
    print "\n\n"
    done = False
    x = order.x.copy()
    y = order.y.copy()
    indices = numpy.array([True]*x.size)
    iteration = 0
    while not done and iteration < numiters:
      done = True
      iteration += 1
      smoothed = FittingUtilities.savitzky_golay(y, window_size, smoothorder)
      residuals = y/smoothed
      if iteration == 1:
        #Save the first and last points, for interpolation
        first, last = smoothed[0], smoothed[-1]
      mean = residuals.mean()
      std = residuals.std()
      print residuals.size, x.size, y.size
      #plt.plot((residuals - mean)/std)
      #plt.show()
      badindices = numpy.where( numpy.logical_or( (residuals-mean)/std > highreject, (residuals-mean)/std < -lowreject ) )[0]
      if badindices.size > 1 and y.size - badindices.size > 2*window_size and iteration < numiters:
        done = False
        x = numpy.delete(x, badindices)
        y = numpy.delete(y, badindices)

    print "iter = %i" %iteration
    if x[0] > order.x[0]:
      x = numpy.append(numpy.array([order.x[0]]), x)
      smoothed = numpy.append(numpy.array([first]), smoothed)
    if x[-1] < order.x[-1]:
      x = numpy.append(x, numpy.array([order.x[-1]]))
      smoothed = numpy.append(smoothed, numpy.array([last]))
    print x.size, y.size, smoothed.size
    smooth_fcn = interp(x, smoothed, k=1)
    smoothed = smooth_fcn(order.x)

    order.cont = FittingUtilities.Continuum(order.x, order.y)
    plt.figure(1)
    plt.plot(order.x, order.y/order.cont, 'k-')
    plt.plot(order.x, smoothed/order.cont, 'r-')
    #plt.figure(2)
    #plt.plot(order.x, order.y/smoothed)
    #plt.plot(order.x, smoothed)
    #orders[i].y /= smoothed
    column = {"wavelength": order.x,
              "flux": order.y/smoothed,
              "continuum": numpy.ones(order.x.size),
              "error": order.err}
    column_list.append(column)
  plt.show()
  outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
  print "Outputting to %s" %outfilename
  HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')



def Filter(fname, vsini=100, numiters=100, lowreject=3, highreject=3):
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")
  column_list = []
  for i, order in enumerate(orders):
    smoothed = FittingUtilities.IterativeLowPass(order, vsini*units.km.to(units.cm), numiter=numiters, lowreject=lowreject, highreject=highreject)
    plt.plot(order.x, order.y)
    plt.plot(order.x, smoothed)
    plt.show()
    continue
  
    done = False
    x = order.x.copy()
    y = order.y.copy()
    indices = numpy.array([True]*x.size)
    iteration = 0
    while not done and iteration < numiters:
      done = True
      iteration += 1
      smoothed = FittingUtilities.savitzky_golay(y, window_size, smoothorder)
      residuals = y/smoothed
      if iteration == 1:
        #Save the first and last points, for interpolation
        first, last = smoothed[0], smoothed[-1]
      mean = residuals.mean()
      std = residuals.std()
      print residuals.size, x.size, y.size
      #plt.plot((residuals - mean)/std)
      #plt.show()
      badindices = numpy.where( numpy.logical_or( (residuals-mean)/std > highreject, (residuals-mean)/std < -lowreject ) )[0]
      if badindices.size > 1 and y.size - badindices.size > 2*window_size and iteration < numiters:
        done = False
        x = numpy.delete(x, badindices)
        y = numpy.delete(y, badindices)

    print "iter = %i" %iteration
    if x[0] > order.x[0]:
      x = numpy.append(numpy.array([order.x[0]]), x)
      smoothed = numpy.append(numpy.array([first]), smoothed)
    if x[-1] < order.x[-1]:
      x = numpy.append(x, numpy.array([order.x[-1]]))
      smoothed = numpy.append(smoothed, numpy.array([last]))
    print x.size, y.size, smoothed.size
    smooth_fcn = interp(x, smoothed, k=1)
    smoothed = smooth_fcn(order.x)

    order.cont = FittingUtilities.Continuum(order.x, order.y)
    plt.figure(1)
    plt.plot(order.x, order.y/order.cont, 'k-')
    plt.plot(order.x, smoothed/order.cont, 'r-')
    #plt.figure(2)
    #plt.plot(order.x, order.y/smoothed)
    #plt.plot(order.x, smoothed)
    #orders[i].y /= smoothed
    column = {"wavelength": order.x,
              "flux": order.y/smoothed,
              "continuum": numpy.ones(order.x.size),
              "error": order.err}
    column_list.append(column)
  plt.show()
  outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
  print "Outputting to %s" %outfilename
  #HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')


if __name__ == "__main__":
  for fname in sys.argv[1:]:
    Smooth(fname, window_size=151, lowreject=5, highreject=5)
    #Filter(fname)
