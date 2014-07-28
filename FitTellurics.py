import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import DataStructures
import Units
from astropy import units, constants
import HelperFunctions
import FittingUtilities
import MakeModel
import GetAtmosphere

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



def FindOrderNums(orders, wavelengths):
  """
    Given a list of xypoint orders and
    another list of wavelengths, this
    finds the order numbers with the
    requested wavelengths
  """
  nums = []
  for wave in wavelengths:
    for i, order in enumerate(orders):
      if order.x[0] < wave and order.x[-1] > wave:
        nums.append(i)
        break
  return nums



if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetObservatory("CTIO")
 
  fileList = []
  start = 0
  end = 999
  makenew = True
  edit_atmosphere=False
  humidity_low = 1.0
  humidity_high = 99.0
  for arg in sys.argv[1:]:
    if "-atmos" in arg:
      edit_atmosphere = True
    elif "-hlow" in arg:
      humidity_low = float(arg.split("=")[1])
    elif "-hhigh" in arg:
      humidity_high = float(arg.split("=")[1])
    else:
      fileList.append(arg)


  #START LOOPING OVER INPUT FILES
  for fname in fileList:
    logfile = open("fitlog_%s.txt" %(fname.split(".fits")[0]), "a")
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name

    #Read file
    orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")

    header = pyfits.getheader(fname)
    angle = float(header["ZD"])
    resolution = 80000.0
    humidity = max(header["OUTHUM"], 5)
    pressure = header["OUTPRESS"]
    temperature = header["OUTTEMP"] + 273.15

    if edit_atmosphere:
      filenames = [f for f in os.listdir("./") if "GDAS" in f]      
      height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'].split("T")[0], header['ut'])

      fitter.EditAtmosphereProfile("Temperature", height, Temp)
      fitter.EditAtmosphereProfile("Pressure", height, Pres)
      fitter.EditAtmosphereProfile("H2O", height, h2o)
      
    
    #Adjust fitter values
    fitter.FitVariable({"h2o": humidity})
    #                    "temperature": temperature})
    #                    "o2": 2.12e5})
    fitter.AdjustValue({"angle": angle,
                        "pressure": pressure,
                        "temperature": temperature,
                        "resolution": resolution,
                        "o2": 2.12e5})
    fitter.SetBounds({"h2o": [humidity_low, humidity_high],
                      "temperature": [temperature-10, temperature+10],
                      "o2": [5e4, 1e6],
                      "resolution": [70000, 90000]})
    
    #Ignore the interstellar sodium D lines and parts of the O2 bands
    fitter.IgnoreRegions(badregions)
    models = []
    
    
    # Determine the H2O abundance
    resolution = []
    h2o = []
    T = []
    o2 = []
    waveshifts = []
    wave0 = []
    chisquared = []
    #for i in [27, 28, 36, 37]:
    #for i in [42, 45, 46, 47]:
    for i in FindOrderNums(orders, [595, 650, 717, 726]):
      print "\n***************************\nFitting order %i: " %(i)
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_source=True, 
                                     return_resolution=True,
				     adjust_wave="model",
				     wavelength_fit_order=3)
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      h2o.append(fitter.GetValue("h2o"))
      T.append(fitter.GetValue("temperature"))
      
      #idx = fitter.parnames.index("h2o")
      #h2o.append(fitter.const_pars[idx])
      #idx = fitter.parnames.index("temperature")
      #T.append(fitter.const_pars[idx])
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average humidity (weight by chi-squared)
    humidity = np.sum(np.array(h2o)*np.array(chisquared)) / np.sum(chisquared)
    temperature = np.sum(np.array(T)*np.array(chisquared)) / np.sum(chisquared)
    logfile.write("Humidity/Temperature values and their chi-squared values:\n")
    for h, t, c in zip(h2o, T, chisquared):
      logfile.write("%g\t%g\t%g\n" %(h, t, 1.0/c))
    logfile.write("Best fit humidity = %g\n" %humidity)
    logfile.write("Best fit temperature = %g\n\n" %temperature)
    fitter.AdjustValue({"h2o": humidity,
                        "temperature": temperature})
    
    # Now, determine the O2 abundance
    fitter.FitVariable({"o2": 2.12e5})
    #for i in [33, 41]:
    for i in FindOrderNums(orders, [630, 690]):
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_source=True,
				     return_resolution=True,
                                     adjust_wave="model",
				     wavelength_fit_order=3)
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      o2.append(fitter.GetValue("o2"))

      #idx = fitter.parnames.index("o2")
      #o2.append(fitter.const_pars[idx])
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average of the other parameter values
    chi2 = np.array(chisquared)
    o2 = np.array(o2)
    resolution = np.array(resolution)
    waveshifts = np.array(waveshifts)
    wave0 = np.array(wave0)
    velshifts = waveshifts/wave0 * constants.c.cgs.value*units.cm.to(units.km)
    vel = np.sum(velshifts*chi2) / np.sum(chi2)
    logfile.write("resolution, velocity shifts and their chi-squared\n")
    for R, v, c in zip(resolution, velshifts, chi2):
      logfile.write("%g\t%g\t%g\n" %(R, v, 1.0/c))
    logfile.write("O2 abundance and their chi-squared:\n")
    for o, c in zip(o2, chi2[-2:]):
      logfile.write("%g\t%g\n" %(o, 1.0/c))
    o2 = np.sum(o2*chi2[-2:])/np.sum(chi2[-2:])
    resolution = np.sum(resolution[:-2]*chi2[:-2])/np.sum(chi2[:-2])
    logfile.write("Best fit o2 mixing ratio = %g ppmv\n" %o2)
    logfile.write("Best fit resolution = %g\n" %resolution)
    logfile.write("Best fit velocity shift = %g km/s\n" %vel)
    """
    
    o2 = 225494
    humidity = 20.9806
    resolution = 78959.6
    vel = 4.4297
    wave0 = orders[0].x.mean()
    fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
    """

    # Finally, apply these parameters to all orders in the data
    for i, order in enumerate(orders):
      print "\n\nGenerating model for order %i of %i\n" %(i, len(orders))
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0,
                          "o2": o2,
                          "h2o": humidity,
                          "resolution": resolution})
      fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      fitter.ImportData(order)
      fitter.resolution_fit_mode = "gauss"
      #fitter.resolution_fit_mode = "svd"
      #wave0 = order.x.mean()
      #fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
      print "fitter.shift = ", fitter.shift
      primary, model = fitter.GenerateModel(fitpars, 
                                            separate_primary=True, 
                                            return_resolution=False)
      
      data = fitter.data
      if min(model.y) > 0.98:
        #The wavelength calibration might be off
        wave0 = order.x.mean()
        fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
        model = fitter.GenerateModel(fitpars, separate_primary=False, nofit=True)
        model.x /= (1.0 + vel/(constants.c.cgs.value*units.cm.to(units.km)))
        model = FittingUtilities.RebinData(model, order.x)
        data = order.copy()
        data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=3, lowreject=2, highreject=5)

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
