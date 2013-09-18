import numpy
import sys
import os
import matplotlib.pyplot as plt
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import FitsUtils
import DataStructures
import Units
from astropy import units, constants
import FindContinuum

homedir = os.environ["HOME"]
linelist = homedir + "/School/Research/Useful_Datafiles/Linelist_visible.dat"


if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetTelluricLineListFile(linelist)
  fitter.SetObservatory("CTIO")
  LineList = numpy.loadtxt(linelist, usecols=(0,))
  logfile = open("fitlog.txt", "w")
 
  fileList = []
  start = 0
  makenew = True
  for arg in sys.argv[1:]:
    if "-start" in arg:
      makenew = False
      start = int(arg.split("=")[-1])
    else:
      fileList.append(arg)


  #START LOOPING OVER INPUT FILES
  for fname in fileList:
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name

    #Read file
    orders = FitsUtils.MakeXYpoints(fname, errors="error", extensions=True, x="wavelength", y="flux")

    #For some reason, this data goes well below zero (but not just in noise). BIAS error or something?
    #Anyways, just add a constant to each order so that the lowest point is a 0
    lowpoint = 9e30
    for order in orders:
      minimum = order.y.min()
      if minimum < lowpoint:
        lowpoint = minimum

    header = pyfits.getheader(fname)
    angle = float(header["ZD"])
    resolution = 80000.0
    humidity = header["OUTHUM"]
    pressure = header["OUTPRESS"]
    temperature = header["OUTTEMP"] + 273.15
    
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
    models = []
    
    #Make a test model, to determine whether/how to fit each value
    fitter.AdjustValue({"wavestart": orders[0].x[0]-20,
                        "waveend": orders[-1].x[-1]+20})
    fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
    fitter.DisplayVariables()
    test_model = fitter.GenerateModel(fitpars, LineList, nofit=True)
    numpy.savetxt("Test_Model.dat", numpy.transpose((test_model.x, test_model.y)), fmt="%.8f")
    
    print "Starting at order %i" %start
    #START LOOPING OVER ORDERS
    column_list = []
    header_list = []
    for i, order in enumerate(orders[start:]):
      print "\n***************************\nFitting order %i: " %(i+start)
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      fitter.FitVariable({"h2o": humidity, 
                          "o2": 2.12e5})

      if lowpoint < 0:
        order.y -= lowpoint
      
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      
      fitter.ImportData(order)

      #Determine how to fit the data from the initial model guess
      left = numpy.searchsorted(test_model.x, order.x[0])
      right = numpy.searchsorted(test_model.x, order.x[-1])
      
      model = test_model[left:right].copy()
      model_amplitude = 1.0 - min(model.y)
      print "Model amplitude: %g" %model_amplitude
      if model_amplitude < 0.01:
        logfile.write("Skipping order %i\n" %(i+start))
        print "Skipping order %i" %(i+start)
        data = order.copy()
        fitter.resolution_fit_mode = "gauss"
        model = fitter.GenerateModel(fitpars, LineList)
      
      elif model_amplitude >= 0.01 and model_amplitude < 1:
        logfile.write("Fitting order %i with guassian line profiles\n" %(i+start)) 
        print "Fitting line profiles with gaussian profile"
	try:
	  model = fitter.Fit(resolution_fit_mode="gauss", fit_primary=False, adjust_wave="model")
	except ValueError:
	  model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))
	
	models.append(model)
        data = fitter.data
      
      else: 
        logfile.write("Fitting order %i with SVD\n" %(i+start))
        print "Large model amplitude. Using SVD for line profiles"
	try:
          model = fitter.Fit(resolution_fit_mode="SVD", fit_primary=False, adjust_wave="model")
	except ValueError:
	  model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))

	models.append(model)
        data = fitter.data

      logfile.write("Array sizes: wave, flux, cont, error, model, primary\n")
      logfile.write("%i\n%i\n%i\n%i\n%i\n%i\n\n\n" %(data.x.size, data.y.size, data.cont.size, data.err.size, model.y.size, primary.y.size))
      
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
	         "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
		 "model": model.y,
                 "primary": primary.y}
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
      column_list.append(columns)
      header_list.append(header_info)
      
      if i == 0 and makenew:
        FitsUtils.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info,], mode="new")
      else:
        FitsUtils.OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info,], mode="append")
    
    

  logfile.close()
