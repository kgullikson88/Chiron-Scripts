
import matplotlib.pyplot as plt
import numpy
import FittingUtilities
import os
import sys
import MakeModel
from astropy.io import fits as pyfits
import HelperFunctions
import warnings
import time

if __name__ == "__main__":
  fileList = []
  plot = False
  single_template = False
  find_template=True
  for arg in sys.argv[1:]:
    if "-p" in arg:
      plot = True
    else:
      fileList.append(arg)

  allfiles = os.listdir("./")
  for fname in fileList:
    good = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    hdulist = pyfits.open(fname)
    fname2 = fname.replace("echi", "Corrected_echi")
    if fname2 not in allfiles:
      warnings.warn("Correction file for %s not found! Skipping!" %fname)
      continue
    bad = HelperFunctions.ReadFits(fname2, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    modelorders = HelperFunctions.ReadFits(fname2, extensions=True, x="wavelength", y="model", cont="continuum", errors="error")
   
    
    column_list = []
    header_list = []
    for i, order in enumerate(good):

      #Super-sample the data for more accuracy
      factor = 1.0
      goodorder = order.copy()
      goodorder.x = numpy.arange(1, order.size()+1)
      badorder = bad[i].copy()
      badorder.x = numpy.arange(1, badorder.size()+1)


      # Find the best pixel shift to make the data line up
      bestchisq = 9e9
      bestshift = 0
      size = goodorder.size()
      left = size/2
      right = left+20
      for shift in range(-50, 50):
        chisq = numpy.sum((goodorder.y[left:right] - badorder.y[left+shift:right+shift])**2)
        if chisq < bestchisq:
          bestchisq = chisq
          bestshift = shift

      # Reduce the size of the original and model files, if necessary
      print "Best shift = ", bestshift
      if abs(bestshift) > 45:
        warnings.warn("Very high shift on order %i (%i pixels)" %(i, bestshift))
      if bestshift < 0:
        order = order[abs(bestshift):]
        bad[i] = bad[i][:bestshift]
        model = modelorders[i][:bestshift]
      elif bestshift > 0:
        order = order[:-bestshift]
        bad[i] = bad[i][bestshift:]
        model = modelorders[i][bestshift:]
      else:
        model = modelorders[i]
      bad[i].x = order.x.copy()
      model.x = order.x.copy()

      if plot:
        plt.figure(1)
        plt.plot(order.x, order.y, 'k-')
        plt.plot(bad[i].x, bad[i].y, 'r-')
        #plt.plot(model.x, model.y, 'g-')
        plt.figure(2)
        plt.plot(order.x, order.y/(order.cont*model.y))
      
      #Change the model when it is very low
      order.y[order.y/order.cont < 1e-5] = 1e-5*order.cont[order.y/order.cont < 1e-5]
      badindices = numpy.where(numpy.logical_or(order.y <= 0, model.y < 0.05))[0]
      model.y[badindices] = order.y[badindices]/order.cont[badindices]
      model.y[model.y < 1e-5] = 1e-5

      # Save data
      columns = {"wavelength": order.x,
                 "flux": order.y/model.y,
                 "continuum": order.cont,
                 "error": order.err/model.y}
      column_list.append(columns)
      

    if plot:
      plt.show()
    
    outfilename = "%s_telluric_corrected.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")


