import FitsUtils
import DataStructures
import FindContinuum
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import numpy
import pylab
from astropy.io import fits as pyfits
from collections import defaultdict
import os
import HelperFunctions
import FittingUtilities


def MedianAdd(fileList, outfilename=None):
  all_data = []
  numorders = []
  medians = []
  for fname in fileList:
    observation = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    all_data.append(observation)
    numorders.append(len(observation))
    medians.append([numpy.median(order.y) for order in observation])

  if any(n != numorders[0] for n in numorders):
    print "Error! Some of the files had different numbers of orders!"
    for i in range(len(fileList)):
      print fileList[i], numorders[i]
    sys.exit()

  #If we get this far, all is well. Add each order indidually
  numorders = numorders[0]
  if outfilename == "None":
    outfilename = "Total.fits"
  column_list = []
  for i in range(numorders):
    x = all_data[0][i].x
    total = numpy.zeros((len(all_data), x.size))
    error = numpy.zeros(x.size)
    norm = 0.0
    for j, observation in enumerate(all_data):
      observation[i].y[observation[i].y < 0.0] = 0.0
      flux = interp(observation[i].x, observation[i].y/medians[j][i])
      error += interp(observation[i].x, observation[i].err**2, k=1)(x)
      total[j] = flux(x)
      norm += medians[j][i]
    
    #for j in range(total.shape[0]):
    #  pylab.plot(x, total[j])
    #pylab.show()
    flux = numpy.median(total, axis=0)*norm
    cont = FittingUtilities.Continuum(x, flux, fitorder=3, lowreject=1.5, highreject=5)
    #Set up data structures for OutputFitsFile
    columns = {"wavelength": x,
               "flux": flux,
               "continuum": cont,
               "error": numpy.sqrt(error)}
    column_list.append(columns)

    pylab.plot(x, flux/cont)
    #pylab.plot(total.x, total.cont)

  print "Outputting to %s" %outfilename
  pylab.show()
  FitsUtils.OutputFitsFileExtensions(column_list, fileList[0], outfilename, mode="new")
  
  #Add the files used to the primary header of the new file
  hdulist = pyfits.open(outfilename, mode='update')
  header = hdulist[0].header
  for i in range(len(fileList)):
    header.set("FILE%i" %(i+1), fileList[i], "File %i used in Co-Adding" %(i+1))
  hdulist[0].header = header
  hdulist.flush()
  hdulist.close()
  


def Add(fileList, outfilename=None):
  all_data = []
  numorders = []
  for fname in fileList:
    observation = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    all_data.append(observation)
    numorders.append(len(observation))

  if any(n != numorders[0] for n in numorders):
    print "Error! Some of the files had different numbers of orders!"
    for i in range(len(fileList)):
      print fileList[i], numorders[i]
    sys.exit()

  #If we get this far, all is well. Add each order indidually
  numorders = numorders[0]
  if outfilename == None:
    outfilename = "Total.fits"
  column_list = []
  for i in range(numorders):
    total = all_data[0][i].copy()
    total.y[total.y < 0.0] = 0.0
    total.err = total.err**2
    for observation in all_data[1:]:
      observation[i].y[observation[i].y < 0.0] = 0.0
      flux = interp(observation[i].x, observation[i].y)
      error = interp(observation[i].x, observation[i].err**2, k=1)
      total.y += flux(total.x)
      total.err += error(total.x)
    total.err = numpy.sqrt(total.err)
    total.cont = FindContinuum.Continuum(total.x, total.y, fitorder=3, lowreject=1.5, highreject=5)

     #Set up data structures for OutputFitsFile
    columns = {"wavelength": total.x,
               "flux": total.y,
               "continuum": total.cont,
               "error": total.err}
    column_list.append(columns)

    pylab.plot(total.x, total.y/total.cont)
    #pylab.plot(total.x, total.cont)

  print "Outputting to %s" %outfilename
  pylab.show()
  HelperFunctions.OutputFitsFileExtensions(column_list, fileList[0], outfilename, mode="new")
  
  #Add the files used to the primary header of the new file
  hdulist = pyfits.open(outfilename, mode='update')
  header = hdulist[0].header
  for i in range(len(fileList)):
    header.set("FILE%i" %(i+1), fileList[i], "File %i used in Co-Adding" %(i+1))
  hdulist[0].header = header
  hdulist.flush()
  hdulist.close()




if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  if len(fileList) > 1:
    Add(fileList)
  else:
    allfiles = [f for f in os.listdir("./") if f.startswith("echi") and "-" in f]
    fileDict = defaultdict(list)
    for fname in allfiles:
      header = pyfits.getheader(fname)
      starname = header['OBJECT'].replace(" ", "_")
      fileDict[starname].append(fname)
    for star in fileDict.keys():
      MedianAdd(fileDict[star], outfilename="%s.fits" %star)
    
