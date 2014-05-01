import numpy
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

plot = True

def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10, normalize=True):
  denoised = HelperFunctions.Denoise(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters)
  if normalize:
    denoised.y /= denoised.y.max()
  return denoised



def roundodd(num):
  rounded = round(num)
  if rounded%2 != 0:
    return rounded
  else:
    if rounded > num:
      return rounded - 1
    else:
      return rounded + 1



def CrossValidation(order, smoothorder=5, lowreject=3, highreject=3, numiters=10, normalize=True):
  """
    Determine the best window size with cross-validation
  """
  
  #plt.plot(order.x, order.y)
  # First, find outliers by doing a guess smooth
  smoothed = SmoothData(order, windowsize=41, normalize=False)
  temp = smoothed.copy()
  temp.y = order.y/smoothed.y
  temp.cont = FittingUtilities.Continuum(temp.x, temp.y, lowreject=2, highreject=2, fitorder=3)
  outliers = HelperFunctions.FindOutliers(temp, numsiglow=6, numsighigh=6, expand=10)
  data = order.copy()
  if len(outliers) > 0:
    #order.y[outliers] = order.cont[outliers]
    order.y[outliers] = smoothed.y[outliers]
  
  #plt.plot(order.x, order.y)
  #plt.plot(order.x, order.cont)
  #plt.show()


  # First, split the data into a training sample and validation sample
  # Use every 10th point for the validation sample
  cv_indices = range(6, order.size()-1, 6)

  training = DataStructures.xypoint(size=order.size()-len(cv_indices))
  validation = DataStructures.xypoint(size=len(cv_indices))
  cv_idx = 0
  tr_idx = 0
  for i in range(order.size()):
    if i in cv_indices:
      validation.x[cv_idx] = order.x[i]
      validation.y[cv_idx] = order.y[i]
      validation.cont[cv_idx] = order.cont[i]
      validation.err[cv_idx] = order.err[i]
      cv_idx += 1
    else:
      training.x[tr_idx] = order.x[i]
      training.y[tr_idx] = order.y[i]
      training.cont[tr_idx] = order.cont[i]
      training.err[tr_idx] = order.err[i]
      tr_idx += 1

  #Rebin the training set to constant wavelength spacing
  xgrid = numpy.linspace(training.x[0], training.x[-1], training.size())
  training = FittingUtilities.RebinData(training, xgrid)
  dx = training.x[1] - training.x[0]

  #Find the rough location of the best window size
  windowsizes = numpy.logspace(-1.5, 1.0, num=20)
  chisq = []
  for i, windowsize in enumerate(windowsizes):
    npixels = roundodd(windowsize/dx)
    if npixels > training.size:
      windowsizes = windowsizes[:i]
      break
    smoothed = FittingUtilities.Iterative_SV(training.y, npixels, smoothorder, lowreject, highreject, numiters)
    smooth_fcn = interp(training.x, smoothed)
    predict = smooth_fcn(validation.x)
    sig = validation.err
    chisq.append(numpy.sum((predict - validation.y)**2/sig**2)/float(predict.size))
  #plt.loglog(windowsizes, chisq)
  #plt.show()

  chisq = numpy.array(chisq)
  idx = numpy.argmin(abs(chisq-1.0))
  left, right = HelperFunctions.GetSurrounding(chisq, 1, return_index=True)
  if left > right:
    temp = left
    left = right
    right = temp

  #Refine the window size to get more accurate
  windowsizes = numpy.logspace(numpy.log10(windowsizes[left]), numpy.log10(windowsizes[right]), num=10)
  chisq = []
  for i, windowsize in enumerate(windowsizes):
    npixels = roundodd(windowsize/dx)
    if npixels > training.size:
      windowsizes = windowsizes[:i]
      break
    smoothed = FittingUtilities.Iterative_SV(training.y, npixels, smoothorder, lowreject, highreject, numiters)
    smooth_fcn = interp(training.x, smoothed)
    predict = smooth_fcn(validation.x)
    sig = validation.err
    chisq.append(numpy.sum((predict - validation.y)**2/sig**2)/float(predict.size))

  chisq = numpy.array(chisq)
  idx = numpy.argmin(abs(chisq-1.0))

  windowsize = windowsizes[idx]
  npixels = roundodd(windowsize/dx)
  smoothed = order.copy()
  smoothed.y = FittingUtilities.Iterative_SV(order.y, npixels, smoothorder, lowreject, highreject, numiters)

  #plt.plot(data.x, data.y)
  #plt.plot(smoothed.x, smoothed.y)
  #plt.show()

  if normalize:
    smoothed.y /= smoothed.y.max()
  return smoothed, windowsize









def GPSmooth(data, low=0.1, high=10, debug=False):
  """
  This will smooth the data using Gaussian processes. It will find the best
  smoothing parameter via cross-validation to be between the low and high.

  The low and high keywords are reasonable bounds for  A and B stars with 
  vsini > 100 km/s.
  """

  smoothed = data.copy()

  # First, find outliers by doing a guess smooth
  smoothed = SmoothData(data, normalize=False)
  temp = smoothed.copy()
  temp.y = data.y/smoothed.y
  temp.cont = FittingUtilities.Continuum(temp.x, temp.y, lowreject=2, highreject=2, fitorder=3)
  outliers = HelperFunctions.FindOutliers(temp, numsiglow=3, expand=5)
  if len(outliers) > 0:
    data.y[outliers] = smoothed.y[outliers]
    
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
      #denoised, theta = GPSmooth(order.copy())
      denoised, theta = CrossValidation(order.copy(), 5, 2, 2, 10)
      #if i == 0:
      #  denoised, theta = CrossValidation(order, 5, 2, 2, 10)
      #else:
      #  dx = order.x[1] - order.x[0]
      #  npixels = roundodd(theta/dx)
      #  denoised = SmoothData(order, npixels, 5, 2, 2, 10)
      print "Window size = %.4f nm" %theta
      #smoothed, lmbd = scikits.datasmooth.smooth_data(order.x, order.y, weights=order.err)
      #plt.plot(order.x, order.y)
      #plt.plot(order.x, smoothed)
      #plt.show()
      #sys.exit()

      column = {"wavelength": denoised.x,
                "flux": order.y / denoised.y,
                "continuum": denoised.cont,
                "error": denoised.err}
      header_list.append((("Smoother", theta, "Smoothing Parameter"),))
      column_list.append(column)
      if plot:
        plt.figure(1)
        plt.plot(order.x, order.y/order.y.mean())
        plt.plot(denoised.x, denoised.y/denoised.y.mean())
        plt.figure(2)
        plt.plot(order.x, order.y/denoised.y)
        plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
        #plt.show()
    if plot:
      plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new', headers_info=header_list)
