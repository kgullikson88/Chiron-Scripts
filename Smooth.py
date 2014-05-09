import numpy
import FittingUtilities
import HelperFunctions
import matplotlib.pyplot as plt
import sys
import os
from astropy import units
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from scipy.interpolate import UnivariateSpline as smooth
import MakeModel
import HelperFunctions
from collections import Counter
from sklearn.gaussian_process import GaussianProcess
from sklearn import cross_validation
from scipy.stats import gmean
from astropy.io import fits, ascii



def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10, expand=0, normalize=True):
  denoised = HelperFunctions.Denoise(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters, expand=expand)
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


def cost(data, prediction, scale = 1, dx=1):
  retval = numpy.sum((prediction - data)**2/scale**2)/float(prediction.size)
  #retval = gmean(data/prediction) / numpy.mean(data/prediction)
  #idx = numpy.argmax(abs(data - prediction))
  #std = numpy.std(data - prediction)
  #retval = abs(data[idx] - prediction[idx]) / std
  #retval = numpy.std(data/(prediction/prediction.sum())) / scale
  #retval = numpy.std(data - prediction)/numpy.mean(scale)
  return retval# + 1e-10*numpy.mean(numpy.gradient(numpy.gradient(prediction, dx), dx)**2)



def OptimalSmooth(order, normalize=True):
  """
    Determine the best window size with cross-validation
  """

  #Flatten the spectrum
  order.y /= order.cont/order.cont.mean()
  order.err /= order.cont/order.cont.mean()

  #Remove outliers (telluric residuals)
  smoothed = SmoothData(order, windowsize=41, normalize=False)
  temp = smoothed.copy()
  temp.y = order.y/smoothed.y
  temp.cont = FittingUtilities.Continuum(temp.x, temp.y, lowreject=2, highreject=2, fitorder=3)
  outliers = HelperFunctions.FindOutliers(temp, numsiglow=6, numsighigh=6, expand=10)
  data = order.copy()
  if len(outliers) > 0:
    #order.y[outliers] = order.cont[outliers]
    order.y[outliers] = smoothed.y[outliers]
    order.err[outliers] = 9e9

  #Make cross-validation sets
  inp = numpy.transpose((order.x, order.err, order.cont))
  X_train, X_test, y_train, y_test = cross_validation.train_test_split(inp, order.y, test_size=0.2)
  X_train = X_train.transpose()
  X_test = X_test.transpose()
  sorter_train = numpy.argsort(X_train[0])
  sorter_test = numpy.argsort(X_test[0])
  training = DataStructures.xypoint(x=X_train[0][sorter_train], y=y_train[sorter_train], err=X_train[1][sorter_train], cont=X_train[2][sorter_train])
  validation = DataStructures.xypoint(x=X_test[0][sorter_test], y=y_test[sorter_test], err=X_test[1][sorter_test], cont=X_test[2][sorter_test])

  """
  #Try each smoothing parameter
  s_array = numpy.logspace(-3, 1, 100)
  chisq = []
  for s in s_array:
    fcn = smooth(training.x, training.y, w=1.0/training.err, s=s)
    prediction = fcn(validation.x)
    chisq.append(cost(validation.y, prediction, validation.err))
    print s, chisq[-1]


  idx = numpy.argmin(numpy.array(chisq) - 1.0)
  s = s_array[idx]
  """

  s = 0.9*order.size()
  smoothed = order.copy()
  fcn = smooth(smoothed.x, smoothed.y, w=1.0/smoothed.err, s=s)
  smoothed.y = fcn(smoothed.x)
  plt.plot(order.x, order.y)
  plt.plot(smoothed.x, smoothed.y)
  plt.show()
  return smoothed, s







def CrossValidation(order, smoothorder=5, lowreject=3, highreject=3, numiters=10, normalize=True):
  """
    Determine the best window size with cross-validation
  """

  #order = HelperFunctions.Denoise(order.copy())
  order.y /= order.cont/order.cont.mean()
  
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
  #plt.plot(order.x, order.y)
  #plt.plot(denoised.x, denoised.y)
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
  size = 40
  left = xgrid.size/2 - size
  right = left + size*2
  func = numpy.poly1d(numpy.polyfit(training.x[left:right]-training.x[left+size], training.y[left:right], 5))
  sig = numpy.std(training.y[left:right] - func(training.x[left:right]-training.x[left+size]))
  sig = validation.err*0.8
  #print "New std = ", sig
  #plt.figure(3)
  #plt.plot(training.x[left:right], training.y[left:right])
  #plt.plot(training.x[left:right], func(training.x[left:right]))
  #plt.show()
  #plt.figure(1)

  #Find the rough location of the best window size
  windowsizes = numpy.logspace(-1.3, 0.5, num=20)
  chisq = []
  skip = 0
  for i, windowsize in enumerate(windowsizes):
    npixels = roundodd(windowsize/dx)
    if npixels < 6:
      skip += 1
      continue
    if npixels > training.size:
      windowsizes = windowsizes[:i]
      break
    smoothed = FittingUtilities.Iterative_SV(training.y.copy(), npixels, smoothorder, lowreject, highreject, numiters)
    smooth_fcn = interp(training.x, smoothed)
    predict = smooth_fcn(validation.x)
    #sig = validation.err
    #chisq.append(cost(training.y, smoothed, training.err))
    chisq.append(cost(validation.y, predict, sig, validation.x[1] - validation.x[0]))
    #chisq.append(numpy.sum((predict - validation.y)**2/sig**2)/float(predict.size))
    #sig = numpy.std(smoothed / training.y)
    #chisq.append(numpy.std(predict/validation.y) / sig)
    print "\t", windowsize, chisq[-1]
  #plt.loglog(windowsizes, chisq)
  #plt.show()
  
  windowsizes = windowsizes[skip:]
  chisq = numpy.array(chisq)
  idx = numpy.argmin(abs(chisq-1.0))
  sorter = numpy.argsort(chisq)
  chisq = chisq[sorter]
  windowsizes = windowsizes[sorter]
  left, right = HelperFunctions.GetSurrounding(chisq, 1, return_index=True)

  if left > right:
    temp = left
    left = right
    right = temp
  print windowsizes[left], windowsizes[right]

  #Refine the window size to get more accurate
  windowsizes = numpy.logspace(numpy.log10(windowsizes[left]), numpy.log10(windowsizes[right]), num=10)
  chisq = []
  for i, windowsize in enumerate(windowsizes):
    npixels = roundodd(windowsize/dx)
    if npixels > training.size:
      windowsizes = windowsizes[:i]
      break
    smoothed = FittingUtilities.Iterative_SV(training.y.copy(), npixels, smoothorder, lowreject, highreject, numiters)
    smooth_fcn = interp(training.x, smoothed)
    predict = smooth_fcn(validation.x)
    #sig = validation.err
    #chisq.append(cost(training.y, smoothed, training.err))
    chisq.append(cost(validation.y, predict, sig, validation.x[1] - validation.x[0]))
    #chisq.append(numpy.sum((predict - validation.y)**2/sig**2)/float(predict.size))
    #sig = numpy.std(smoothed / training.y)
    #chisq.append(numpy.std(predict/validation.y) / sig)
    print "\t", windowsize, chisq[-1]

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
  plot = False
  vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" %(os.environ["HOME"])
  for arg in sys.argv[1:]:
    if "-p" in arg:
      plot = True
    elif "-vsini" in arg:
      vsini_file = arg.split("=")[-1]
    else:
      fileList.append(arg)

  #Read in the vsini table
  vsini_data = ascii.read(vsini_file)[10:]

  if len(fileList) == 0:
    fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    
    
    #Find the vsini of this star
    header = fits.getheader(fname)
    starname = header["object"]
    found = False
    for data in vsini_data:
      if data[0] == starname:
        vsini = float(data[1])
        found = True
    if not found:
      sys.exit("Cannot find %s in the vsini data: %s" %(starname, vsini_file))
    print starname, vsini
    
    #Begin looping over the orders
    column_list = []
    header_list = []
    for i, order in enumerate(orders):
      print "Smoothing order %i/%i" %(i+1, len(orders))
      #Fix errors
      order.err[order.err > 1e8] = numpy.sqrt(order.y[order.err > 1e8])

      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = FittingUtilities.RebinData(order, xgrid)
      
      dx = order.x[1] - order.x[0]
      smooth_factor = 0.8
      theta = roundodd(vsini/3e5 * order.x.mean()/dx * smooth_factor)
      denoised = SmoothData(order, 
                            windowsize=theta, 
                            smoothorder=3, 
                            lowreject=3, 
                            highreject=3,
                            expand=10, 
                            numiters=10)
      #denoised, theta = GPSmooth(order.copy())
      #denoised, theta = CrossValidation(order.copy(), 5, 2, 2, 10)
      #denoised, theta = OptimalSmooth(order.copy())
      #denoised.y *= order.cont/order.cont.mean()
      print "Window size = %.4f nm" %theta


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
        plt.title(starname)
        plt.figure(2)
        plt.plot(order.x, order.y/denoised.y)
        plt.title(starname)
        #plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
        #plt.show()
    if plot:
      plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new', headers_info=header_list)
