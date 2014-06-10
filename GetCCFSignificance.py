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
import numpy
#import HelperFunctions
import sys
from sklearn.gaussian_process import GaussianProcess


def rational_quadratic(x1, x2, pars):
  """
  Rational quadratic correlation function.
  This allows a range of scale lengths,
  which is useful for CCF noise.
  """
  l, alpha = pars[0], pars[1]
  return (1.0 + (x1 - x2)**2 / (2*alpha*l**2))**(-alpha)


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
	vel, corr = numpy.loadtxt(filename, usecols=(0,1), unpack=True)

	