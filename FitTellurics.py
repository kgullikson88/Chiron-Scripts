import numpy
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



def GetAtmosphereFile(header):
  fnames = [f for f in os.listdir("./") if "GDAS_atmosphere" in f]
  datestr = header["DATE-OBS"]
  year = int(datestr.split("-")[0])
  month = int(datestr.split("-")[1])
  day = int(datestr.split("-")[2].split("T")[0])
  hour = int(datestr.split("T")[1].split(":")[0])
  minute = int(datestr.split("T")[1].split(":")[1])
  obstime = year + (month*30 + day)/365.0 + (hour*60 + minute)/(365.0*24.0*60.0)
  #obstime = hour + minute/60.0
  filename = ""
  mindiff = 9e12
  for fname in fnames:
    datestr = fname.split("_")[2]
    year = int(datestr.split("-")[0])
    month = int(datestr.split("-")[1])
    day = int(datestr.split("-")[2])
    hour = (float(fname.split("_")[3]) + float(fname.split("_")[4].split(".")[0]))/2.0
    hour = float(fname.split("_")[3])
    atmos_time = year + (month*30 + day)/365.0 + hour/(365.0*24.0)
    #atmos_time = hour
    if abs(atmos_time - obstime) < mindiff:
      mindiff = abs(atmos_time - obstime)
      filename = fname
  return filename





if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetObservatory("CTIO")
  logfile = open("fitlog.txt", "a")
 
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
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name

    #Read file
    orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")

    header = pyfits.getheader(fname)
    angle = float(header["ZD"])
    resolution = 80000.0
    humidity = header["OUTHUM"]
    pressure = header["OUTPRESS"]
    temperature = header["OUTTEMP"] + 273.15

    if edit_atmosphere:
      #Find the appropriate filename
      atmosphere_fname = GetAtmosphereFile(header)
      
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
    humidity = 10.0
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
                      "resolution": [70000, 120000]})
    
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
    for i in [27, 28, 36, 37]:
      print "\n***************************\nFitting order %i: " %(i)
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_primary=True, 
                                     return_resolution=True,
				     adjust_wave="model")
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = fitter.parnames.index("h2o")
      h2o.append(fitter.const_pars[idx])
      idx = fitter.parnames.index("temperature")
      T.append(fitter.const_pars[idx])
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average humidity (weight by chi-squared)
    humidity = numpy.sum(numpy.array(h2o)*numpy.array(chisquared)) / numpy.sum(chisquared)
    temperature = numpy.sum(numpy.array(T)*numpy.array(chisquared)) / numpy.sum(chisquared)
    logfile.write("Humidity/Temperature values and their chi-squared values:\n")
    for h, t, c in zip(h2o, T, chisquared):
      logfile.write("%g\t%g\t%g\n" %(h, t, 1.0/c))
    logfile.write("\n")
    fitter.AdjustValue({"h2o": humidity,
                        "temperature": temperature})
    
    # Now, determine the O2 abundance
    fitter.FitVariable({"o2": 2.12e5})
    for i in [33, 41]:
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_primary=True,
				     return_resolution=True,
                                     adjust_wave="model")
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = fitter.parnames.index("o2")
      o2.append(fitter.const_pars[idx])
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average of the other parameter values
    chi2 = numpy.array(chisquared)
    o2 = numpy.array(o2)
    resolution = numpy.array(resolution)
    waveshifts = numpy.array(waveshifts)
    wave0 = numpy.array(wave0)
    velshifts = waveshifts/wave0 * constants.c.cgs.value*units.cm.to(units.km)
    vel = numpy.sum(velshifts*chi2) / numpy.sum(chi2)
    logfile.write("resolution, velocity shifts and their chi-squared\n")
    for R, v, c in zip(resolution, velshifts, chi2):
      logfile.write("%g\t%g\t%g\n" %(R, v, 1.0/c))
    logfile.write("O2 abundance and their chi-squared:\n")
    for o, c in zip(o2, chi2[-2:]):
      logfile.write("%g\t%g\n" %(o, 1.0/c))
    o2 = numpy.sum(o2*chi2[-2:])/numpy.sum(chi2[-2:])
    resolution = numpy.sum(resolution[:-2]*chi2[:-2])/numpy.sum(chi2[:-2])
    """
    
    o2 = 224773
    humidity = 42.8873
    resolution = 88472.7266
    vel = 4.8
    """

    # Finally, apply these parameters to all orders in the data
    for i, order in enumerate(orders):
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
      wave0 = order.x.mean()
      fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
      print "fitter.shift = ", fitter.shift
      primary, model = fitter.GenerateModel(fitpars, 
                                            separate_primary=True, 
                                            return_resolution=False)
      if min(model.y) > 0.9:
	# The wavelength calibration might be off.
	fitter.shift = 0.0
	model = fitter.GenerateModel(fitpars, separate_primary=False, nofit=True)

	model.x /= (1.0 + vel/(constants.c.cgs.value*units.cm.to(units.km)))
	xgrid = numpy.linspace(model.x[0], model.x[-1], model.size())
	model = FittingUtilities.RebinData(model, xgrid)
	#order2 = order.copy()
	#p = FittingUtilities.Iterative_SV(order2.y/order2.cont, 61, 4)
	#order2.y /= p
	#model2 = FittingUtilities.RebinData(model, order2.x)
	left = numpy.searchsorted(model.x, order.x[0])
	right = numpy.searchsorted(model.x, order.x[-1])
	if min(model.y[left:right]) < 0.99:
	  #Interpolate to finer spacing
	  oversampling = 5.0
          xgrid = numpy.linspace(order.x[0], order.x[-1], order.size()*oversampling)
          order2 = FittingUtilities.RebinData(order, xgrid)
          model2 = FittingUtilities.RebinData(model, xgrid)
          p = FittingUtilities.Iterative_SV(order2.y/order2.cont, 61*oversampling, 4)
	  order2.y /= p
  	  shift, corr = FittingUtilities.CCImprove(order2, 
	                                           model2,
	                                           be_safe=True,
	                                           tol=0.05,
	                                           debug=True)
	  model.x -= shift
	  #plt.figure(1)
	  #plt.plot(order2.x, order2.y/order2.cont)
	  #plt.plot(order2.x, p)
	  #plt.plot(model2.x, model2.y/model2.cont)
	  #plt.figure(2)
	  #plt.plot(corr.x, corr.y)
	  #plt.show()
	  print "Model has low amplitude! Shifting by %.5g nm" %shift
	model = FittingUtilities.ReduceResolution(model, resolution)
	model = FittingUtilities.RebinData(model, order.x)

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
      
      
  #plt.show()
  logfile.close()
