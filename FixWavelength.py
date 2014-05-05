import FitsUtils
import matplotlib.pyplot as plt
import numpy
import FittingUtilities
import os
import sys
import MakeModel
from astropy.io import fits as pyfits
import HelperFunctions
import warnings

if __name__ == "__main__":
  fileList = []
  single_template = False
  find_template=True
  skip = 7
  for arg in sys.argv[1:]:
    if '-template' in arg:
      template = arg.split("=")[-1]
      single_template = True
      find_template = False
    elif "-skip" in arg:
      skip = int(arg.split("=")[-1])
    else:
      fileList.append(arg)

  #Read in the blaze file (we will need it later)
  blazefile = [f for f in os.listdir("./") if "blaze" in f][0]
  try:
    blazeorders = HelperFunctions.ReadFits(blazefile, extensions=False)[::-1]
  except ValueError:
    blazeorders = HelperFunctions.ReadFits(blazefile, extensions=False, errors=2)[::-1]


  if find_template:
    ach_files =[f for f in os.listdir("./") if f.startswith("achi") and f.endswith("-0.fits")]
    objects = {}
    for fname in ach_files:
      header = pyfits.getheader(fname)
      objects[header['object']] = fname

  if single_template:
    native = FitsUtils.MakeXYpoints(template, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

  for fname in fileList:
    hdulist = pyfits.open(fname)
    if all("WaveFixed" in hdu.header.keys() for hdu in hdulist[1:]):
      print "Wavelength calibration already done. Skipping file %s" %fname
      continue


    if not single_template and not find_template:
      if fname.startswith("a"):
        native = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
        fname = "e" + fname[1:]
        mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
      
      elif fname.startswith("e"):
        mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
        fname2 = "a" + fname[1:]
        native = FitsUtils.MakeXYpoints(fname2, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    elif find_template:
      header = pyfits.getheader(fname)
      obj = header['object']
      try:
        fname2 = fname.replace("echi", "achi")
        native = FitsUtils.MakeXYpoints(fname2, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
      except IOError:
        warnings.warn("Warning! Could not find a template file for %s" %fname)
        continue
        #fname2 = objects[obj]
        #native = FitsUtils.MakeXYpoints(fname2, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
      mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    else:
      mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

    column_list = []
    header_list = []
    for i, order in enumerate(mine[skip:]):

      # Get an estimate for the best pixel shift
      shift, corr = FittingUtilities.CCImprove(native[i], order, debug=True, be_safe=False)
      pixelshift = shift*(native[i].x[-1] - native[i].x[0])/float(native[i].x.size)
      pixelshift = int(numpy.searchsorted(order.x, native[i].x[0]) + pixelshift + 0.5)

      #Super-sample the data for more accuracy
      factor = 1.0
      order.x = numpy.arange(1, order.size()+1)
      xgrid = numpy.arange(1, order.size()+1, 1.0/factor)
      order2 = FittingUtilities.RebinData(order, xgrid)
      xgrid = numpy.arange(1, native[i].size()+1, 1.0/factor)
      native2 = native[i].copy()
      native2.x = numpy.arange(1, native2.size()+1)
      native2 = FittingUtilities.RebinData(native2, xgrid)

      # Find the best pixel shift to make the data line up
      bestchisq = 9e9
      bestshift = 0
      order2.y /= order2.y.mean()
      native2.y /= native2.y.mean()
      size = native2.size()
      sizediff = order2.size() - native2.size()
      for shift in range(pixelshift - 10, pixelshift+10):
        chisq = numpy.sum((order2.y[shift:size+shift] - native2.y)**2)
        if chisq < bestchisq:
          bestchisq = chisq
          bestshift = shift
      order = order[bestshift:size+bestshift]
      #plt.plot(order.x, order.y/order.y.mean(), 'k-')
      order.x = native[i].x
      #plt.plot(order.x, order.y/order.y.mean(), 'k-')
      #plt.plot(native[i].x, native[i].y/native[i].y.mean(), 'r-')

      columns = {"wavelength": order.x,
                 "flux": order.y,
                 "continuum": order.cont,
                 "error": order.err}
      column_list.append(columns)
      header_list.append((("WaveFixed", True, "Wavelength calibration taken from native reductions"),))
    
    #plt.show()
    
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, fname, mode='new', headers_info = header_list)


