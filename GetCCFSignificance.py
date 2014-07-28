"""
Estimates the significance of a CCF peak
using Gaussian Processes to simulate the 
correlated noise.

Inputs:
  The filename of the correlation functions.
      The CCF is assumed to be a 2-column text
      file with velocity and CCF columns
  A set of regions to ignore. Typically, these
      are velocity ranges that should not be
      modeled as CCF noise (the peak you care
      about, and the peaks from the primary star)
      The regions should have the following form:
      ignore=-210to-150,-10to20 (nm wavelengths
      ranges separated by commas, NO SPACES)
"""
import sys
sys.path.insert(0, '/home/kgullikson/.local/lib/python2.7/site-packages')
import numpy as np
#import HelperFunctions

from sklearn.gaussian_process import GaussianProcess
import matplotlib.pyplot as plt


def rational_quadratic(pars, d):
  """
  Rational quadratic correlation function.
  This allows a range of scale lengths,
  which is useful for CCF noise.
  """
  if len(pars.shape) > 1:
    pars = pars[0]
  print pars.shape, "\t", len(pars.shape), "\t", pars
  l, alpha = pars[0], pars[1]
  r = (1.0 + np.sum(d**2, axis=1) / (2*alpha*l**2))**(-alpha)
  return r



if __name__ == "__main__":
  #Parse command-line arguments
  ignore_regions = []
  for arg in sys.argv[1:]:
    if "ignore" in arg:
      ig = arg.split("=")[1]
      segments = ig.split(",")
      for seg in segments:
        r = seg.split("to")
        left = float(r[0])
        right = float(r[1])
        ignore_regions.append((left, right))
    else:
      filename = arg
  print ignore_regions

  #Read in the correlation function
  vel_original, corr_original = np.loadtxt(filename, usecols=(0,1), unpack=True)
  print "Original: ", vel_original.size
  #corr_original /= corr_original.mean()

  #Remove the requested regions in velocity space 
  vel = vel_original.copy()
  corr = corr_original.copy()
  for region in ignore_regions:
    left = np.searchsorted(vel, region[0])
    right = np.searchsorted(vel, region[1])
    indices = np.r_[np.arange(left), np.arange(right, len(vel))]
    vel = vel[indices]
    corr = corr[indices]
  print "Processed: ", vel.size

  #Plot
  f = 0.03
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(vel_original, corr_original, 'k-')
  ax.plot(vel, corr, 'g-')
  ax.errorbar(vel, corr, yerr=np.std(corr)*f)
  #plt.show(); sys.exit()
  
  #Generate the gaussian process fit
  #f = 1
  gp = GaussianProcess(corr=rational_quadratic,
                       theta0=np.array((0.01, 1e5)),
                       thetaL=np.array((1e-3, 0.1)),
                       thetaU=np.array((1e2, 1e7)),
                       normalize=True,
                       nugget=(np.std(corr)*f/corr)**2)
  #gp = GaussianProcess(corr='squared_exponential',
  #                     theta0=1e3,
  #                     thetaL=10,
  #                     thetaU=1e4,
  #                     nugget=(np.std(corr)*f/corr)**2)
  
  #Fit the theta parameters, with the data excluded
  gp.fit(vel[:,None], corr)

  #Now, use those parameters to fit the whole GP. Can it reproduce the peaks?
  gp2 = GaussianProcess(corr=rational_quadratic,
                        theta0=gp.theta_,
                        nugget=(np.std(corr_original)*f/corr_original)**2)
  gp2.fit(vel_original[:,None], corr_original)

  #prediction, error = gp2.predict(vel_original[:,None], eval_MSE=True)
  prediction, error = gp.predict(vel_original[:,None], eval_MSE=True)

  error = np.sqrt(error)
  print "Best theta = ", gp.theta_
  ax.plot(vel_original, prediction, 'r-')
  ax.fill_between(vel_original, prediction - 3 * error, prediction + 3 * error, color='red', alpha=0.4)

  mu = np.mean(corr)
  std = np.std(corr)
  ax.fill_between(vel_original, mu-3*std, mu+3*std, color='green', alpha=0.3)
  ax.set_xlabel("Velocity (km/s)")
  ax.set_ylabel("CCF")



  plt.show()

