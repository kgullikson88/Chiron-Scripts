import FitsUtils
import matplotlib.pyplot as plt
import numpy
import FittingUtilities
import os
import sys
import MakeModel
from astropy.io import fits as pyfits
import HelperFunctions

if __name__ == "__main__":
  fileList = []
  single_template = False
  find_template=False
  skip = 7
  for arg in sys.argv[1:]:
    if '-template' in arg:
      template = arg.split("=")[-1]
      single_template = True
    elif "-find" in arg:
      find_template=True
      ach_files =[f for f in os.listdir("./") if f.startswith("achi") and f.endswith("-0.fits")]
      objects = {}
      for fname in ach_files:
        header = pyfits.getheader(fname)
        objects[header['object']] = fname
    elif "-skip" in arg:
      skip = int(arg.split("=")[-1])
    else:
      fileList.append(arg)


  if single_template:
    native = FitsUtils.MakeXYpoints(template, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

  for fname in fileList:
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
        fname2 = objects[obj]
        native = FitsUtils.MakeXYpoints(fname2, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
      mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    else:
      mine = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

    column_list = []
    header_list = []
    for i, order in enumerate(mine[skip:]):
      shift, corr = FittingUtilities.CCImprove(native[i], order, debug=True, be_safe=False)
      pixelshift = shift*(native[i].x[-1] - native[i].x[0])/float(native[i].x.size)
      left = int(numpy.searchsorted(order.x, native[i].x[0]) + pixelshift + 0.5)
      if abs(order.x[left] - native[i].x[0]) > abs(order.x[left-1] - native[i].x[0]):
        left -= 1
      right = left + native[i].size()
      
      order = order[left:right]
      columns = {"wavelength": order.x,
                 "flux": order.y,
                 "continuum": order.cont,
                 "error": order.err}
      column_list.append(columns)
      header_list.append((("WaveFixed", True, "Wavelength calibration taken from native reductions"),))

    
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, fname, mode='new', headers_info = header_list)


