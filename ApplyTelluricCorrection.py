import pyfits
import FitsUtils
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import matplotlib.pyplot as plt
import DataStructures
import os
import FindContinuum
import numpy


def ReadCorrectedFile(fname, yaxis="model"):
  orders = []
  headers = []
  hdulist = pyfits.open(fname)
  numorders = len(hdulist)
  for i in range(1, numorders):
    order = hdulist[i].data
    xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                  y=order.field(yaxis),
                                  cont=order.field("continuum"),
                                  err=order.field("error"))

    orders.append(xypt)
    headers.append(hdulist[i].header)
  return orders, headers


def Correct(original, corrected, offset=None, get_primary=False):
  #Read in the data and model
  original_orders = FitsUtils.MakeXYpoints(original, extensions=True, x="wavelength", y="flux", errors="error", cont="continuum")
  corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
  test_orders, header = ReadCorrectedFile(corrected, yaxis="flux")
  for order, model in zip(test_orders, corrected_orders):
    plt.plot(order.x, order.y/order.cont)
    plt.plot(model.x, model.y)
  plt.show()
  if get_primary:
    primary_orders = ReadCorrectedFile(corrected, yaxis="primary")[0]
  if offset == None:
    offset = len(original_orders) - len(corrected_orders)
  for i in range(offset, len(original_orders)):
    data = original_orders[i]
    data.cont = FindContinuum.Continuum(data.x, data.y)
    try:
      model = corrected_orders[i-offset]
      header = corrected_headers[i-offset]
      if get_primary:
        primary = primary_orders[i-offset]
      print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" %(i, header['h2oval'], header['o2val'])
    except IndexError:
      model = DataStructures.xypoint(x=data.x, y=numpy.ones(data.x.size))
      print "Warning!!! Telluric Model not found for order %i" %i

    #plt.plot(data.x, data.y/data.cont)
    #plt.plot(model.x, model.y)
    #plt.show()
    if model.size() < data.size():
      left = numpy.searchsorted(data.x, model.x[0])
      right = numpy.searchsorted(data.x, model.x[-1])
      if right < data.size():
        right += 1
      data = data[left:right]
    elif model.size() > data.size():
      sys.exit("Error! Model size (%i) is larger than data size (%i)" %(model.size(), data.size()))

    badindices = numpy.where(numpy.logical_or(data.y <= 0, model.y < 0.05))[0]
    model.y[badindices] = data.y[badindices]

    plt.plot(data.x, data.y / model.y)
    data.y /= model.y
    if get_primary:
      data.y /= primary.y
    original_orders[i] = data.copy()
  plt.show()
  return original_orders





def main1():
  primary=False
  if len(sys.argv) > 2:
    original = sys.argv[1]
    corrected = sys.argv[2]
    if len(sys.argv) > 3 and "prim" in sys.argv[3]:
      primary=True
  
    outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
    print "Outputting to %s" %outfilename

    corrected_orders = Correct(original, corrected, offset=None, get_primary=primary)

    column_list = []
    for i, data in enumerate(corrected_orders):
      plt.plot(data.x, data.y/data.cont)
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err}
      column_list.append(columns)
    FitsUtils.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")
    
    plt.show()

  else:
    allfiles = os.listdir("./")
    corrected_files = [f for f in allfiles if "Corrected_" in f and f.endswith(".fits")]
    #original_files = [f for f in allfiles if any(f in cf for cf in corrected_files)]

    #print corrected_files
    #print original_files

    for corrected in corrected_files:
      original = [f for f in allfiles if (f in corrected and f != corrected)]
      if len(original) == 1:
        original = original[0]
      else:
        sys.exit("Error! %i matches found to corrected file %s" %(len(original), corrected))

      print corrected, original
      outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
      print "Outputting to %s" %outfilename

      corrected_orders = Correct(original, corrected, offset=None)

      column_list = []
      for i, data in enumerate(corrected_orders):
        plt.plot(data.x, data.y/data.cont)
        #Set up data structures for OutputFitsFile
        columns = {"wavelength": data.x,
                   "flux": data.y,
                   "continuum": data.cont,
                   "error": data.err}
        column_list.append(columns)
      FitsUtils.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")
        
      plt.title(original)
      plt.xlabel("Wavelength (nm)")
      plt.ylabel("Flux")
      plt.show()




if __name__ == "__main__":
  main1()
