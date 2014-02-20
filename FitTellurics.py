import numpy
import sys
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import FitsUtils
import DataStructures
import Units
from astropy import units, constants
import FindContinuum
import HelperFunctions

homedir = os.environ["HOME"]

badregions = [[588.98, 589.037],   #Na D line 1
              [589.567, 589.632],  #Na D line 2
              [627.4, 629.0],  #O2 band
              [686.4, 690.7]]  #O2 band
              
              
namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
                  "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
                  "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
                  "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
                  "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
                  "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
                  "o3": ["O3FIT", "O3VAL", "O3 abundance"],
                  "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
                  "co": ["COFIT", "COVAL", "CO abundance"],
                  "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
                  "o2": ["O2FIT", "O2VAL", "O2 abundance"],
                  "no": ["NOFIT", "NOVAL", "NO abundance"],
                  "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
                  "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
                  "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
                  "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}  


if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetObservatory("CTIO")
  logfile = open("fitlog.txt", "w")
 
  fileList = []
  start = 0
  end = 999
  makenew = True
  edit_atmosphere=False
  for arg in sys.argv[1:]:
    if "-start" in arg:
      makenew = False
      start = int(arg.split("=")[-1])
    elif "-end" in arg:
      end = int(arg.split("=")[-1])
    elif "-atmos" in arg:
      edit_atmosphere = True
      atmosphere_fname = arg.split("=")[-1]
    else:
      fileList.append(arg)


  #START LOOPING OVER INPUT FILES
  for fname in fileList:
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name

    #Read file
    orders = FitsUtils.MakeXYpoints(fname, errors="error", extensions=True, x="wavelength", y="flux")

    header = pyfits.getheader(fname)
    angle = float(header["ZD"])
    resolution = 80000.0
    humidity = header["OUTHUM"]
    pressure = header["OUTPRESS"]
    temperature = header["OUTTEMP"] + 273.15

    if edit_atmosphere:
      #Read in GDAS atmosphere profile information
      Pres,height,Temp,dew = numpy.loadtxt(atmosphere_fname, usecols=(0,1,2,3), unpack=True)
      sorter = numpy.argsort(height)
      height = height[sorter]
      Pres = Pres[sorter]
      Temp = Temp[sorter]
      dew = dew[sorter]
      
      #Convert dew point temperature to ppmv
      Pw = 6.116441 * 10**(7.591386*Temp/(Temp + 240.7263))
      h2o = Pw / (Pres-Pw) * 1e6
      
      height /= 1000.0
      Temp += 273.15
      fitter.EditAtmosphereProfile("Temperature", height, Temp)
      fitter.EditAtmosphereProfile("Pressure", height, Pres)
      fitter.EditAtmosphereProfile("H2O", height, h2o)
      
    
    #Adjust fitter values
    fitter.FitVariable({"h2o": humidity, 
                        "o2": 2.12e5})
    fitter.AdjustValue({"angle": angle,
                        "temperature": temperature,
                        "pressure": pressure,
                        "resolution": resolution})
    fitter.SetBounds({"h2o": [1.0, 98.0],
                      "o2": [5e4, 1e6],
                      "resolution": [resolution/2.0, resolution*2.0]})
    
    #Ignore the interstellar sodium D lines and parts of the O2 bands
    fitter.IgnoreRegions(badregions)
    models = []

    """
    #Make a test model, to determine whether/how to fit each value
    end = min(end, len(orders))
    fitter.AdjustValue({"wavestart": orders[start].x[0]-20,
                        "waveend": orders[end-1].x[-1]+20})
    fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
    fitter.DisplayVariables()
    test_model = fitter.GenerateModel(fitpars, nofit=True)
    numpy.savetxt("Test_Model.dat", numpy.transpose((test_model.x, test_model.y)), fmt="%.8f")
    """

    # Determine the H2O abundance
    resolution = []
    h2o = []
    o2 = []
    waveshifts = []
    wave0 = []
    chisquared = []
    for i in [27, 28, 36, 37]:
      print "\n***************************\nFitting order %i: " %(i)
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_primary=True, 
                                     adjust_wave="model")
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = fitter.parnames.index("h2o")
      h2o.append(fitter.const_pars[idx])
      chisquared.append(1.0/fitter.chisq_vals[-1])

    # Determine the average humidity (weight by chi-squared)
    humidity = numpy.sum(numpy.array(h2o)*numpy.array(chisquared)) / numpy.sum(chisquared)
    fitter.AdjustValue({"h2o": humidity})
    
    # Now, determine the O2 abundance
    for i in [33, 41]:
      order = orders.[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_primary=True, 
                                     adjust_wave="model")
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = fitter.parnames.index("o2")
      o2.append(fitter.const_pars[idx])
      chisquared.append(1.0/fitter.chisq_vals[-1])

    # Determine the average of the other parameter values
    chi2 = numpy.array(chisquared)
    o2 = numpy.array(o2)
    resolution = numpy.array(resolution)
    waveshifts = numpy.array(waveshifts)
    wave0 = numpy.array(wave0)
    o2 = numpy.sum(o2*chi2)/numpy.sum(chi2)
    resolution = numpy.sum(resolution*chi2)/numpy.sum(chi2)
    velshifts = waveshifts/wave0 * constants.c.cgs.value*units.cm.to(units.km)
    vel = numpy.sum(velshifts*chi2) / numpy.sum(chi2)


    # Finally, apply these parameters to all orders in the data
    fitter.shift = 0
    for i, order in enumerate(orders):
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0,
                          "o2": o2,
                          "h2o": h2o,
                          "resolution": resolution})
      fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
      fitter.ImportData(order)
      fitter.resolution_fit_mode = "gauss"
      primary, model = fitter.GenerateModel(fitpars, 
                                            separate_primary=True, 
                                            return_resolution=False)
      data = fitter.data
      data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=3, lowreject=2, highreject=2)

      # Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
                 "model": model.y,
                 "primary": primary.y}
      
      header_info = []
      numpars = len(fitter.const_pars)
      for j in range(numpars):
        try:
          parname = fitter.parnames[j]
          parval = fitter.const_pars[j]
          fitting = fitter.fitting[j]
          header_info.append([namedict[parname][0], fitting, namedict[parname][2] ])
          header_info.append([namedict[parname][1], parval, namedict[parname][2] ])
        except KeyError:
          print "Not saving the following info: %s" %(fitter.parnames[j])
      
      
      if (i == 0 and makenew) or not exists:
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info,], mode="new")
        exists = True
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info,], mode="append")
      
      
   
  logfile.close()
