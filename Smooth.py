import numpy
import FitsUtils
import FittingUtilities
import matplotlib.pyplot as plt
import sys
from astropy import units
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp


if __name__ == "__main__":
  for fname in sys.argv[1:]:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    column_list = []
    for order in orders:
      #Linearize
      yfcn = interp(order.x, order.y)
      errfcn = interp(order.x, order.err)
      contfcn = interp(order.x, order.cont)
      order.x = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order.y = yfcn(order.x)
      order.err = errfcn(order.x)
      order.cont = contfcn(order.x)
      
      order2 = order.copy()
      denoised = FittingUtilities.Denoise3(order2) #, snr=400.0, reduction_factor=0.15)
      denoised.y = FittingUtilities.Iterative_SV(denoised.y, 91, 5, lowreject=1.5, highreject=5)

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
    outfilename = "%s_denoised.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    FitsUtils.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
