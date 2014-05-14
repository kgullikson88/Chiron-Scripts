import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import scipy.signal
import os
import sys
import DataStructures
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy.io import ascii
from astropy import units, constants
import StarData
import SpectralTypeRelations
from PlotBlackbodies import Planck
import FittingUtilities
import Smooth
import HelperFunctions
import Broaden
from Search_Fast import Process_Data
import Sensitivity
from re import search
import Correlate

#Ensure a directory exists. Create it if not
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        

homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define some constants to use
vsini = 20*units.km.to(units.cm)
resolution = 80000

#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[567.5, 575.5],
              [588.5, 598.5],
              [627, 632],
              [647,655],
              [686, 706],
              [716, 734],
              [759, 9e9]]

#Set up model list
model_list = [ modeldir + "lte30-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte32-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte34-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte35-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte36-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte37-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte38-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte39-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte40-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte42-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte44-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte46-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte48-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte50-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte51-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte52-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte53-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte54-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte55-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte56-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte57-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte58-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte59-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
""",
               modeldir + "lte61-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte62-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte63-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte64-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte65-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte66-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte72-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte76-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]"""
   
              
"""
star_list = []
temp_list = []
gravity_list = []
metal_list = []
model_data = []
start = 0
end = len(model_list)
fitorders = 8*numpy.ones(len(model_list))
fitorders[0] = 7
fitorders[1] = 6
fitorders[5] = 7
if __name__ == "__main__":
  for i, fname in enumerate(model_list[start:end]):
    if "PHOENIX2004" in fname:
      temp = int(fname.split("lte")[-1][:2])*100
      gravity = float(fname.split("lte")[-1][3:6])
      metallicity = float(fname.split("lte")[-1][6:10])
    elif "PHOENIX-ACES" in fname:
      temp = int(fname.split("lte")[-1][:2])*100
      gravity = float(fname.split("lte")[-1][3:7])
      metallicity = float(fname.split("lte")[-1][7:11])
    print "Reading in file %s" %fname
    x,y = numpy.loadtxt(fname, usecols=(0,1), unpack=True)
    #c = FittingUtilities.Continuum(x, y, fitorder=fitorders[i+start], lowreject=2, highreject=7)
  
    plt.plot(x,y)
    for order in range(6,11):
      print order
      c = FittingUtilities.Continuum(x, y, fitorder=order, lowreject=2, highreject=7, numiter=10)
      plt.plot(x,c, label="%ith order" %order)
    plt.title(fname.split("/")[-1])
    plt.legend(loc='best')
    plt.show()
    continue
    model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y, cont=10**c)
    model = FittingUtilities.RebinData(model, numpy.linspace(model.x[0], model.x[-1], model.size()))
    model = Broaden.RotBroad(model, vsini)
    model = Broaden.ReduceResolution2(model, resolution)
    model_data.append(interp(model.x, model.y/model.cont))
    #model_data.append( model )
    star_list.append(str(temp))
    temp_list.append(temp)
    gravity_list.append(gravity)
    metal_list.append(metallicity)
"""


MS = SpectralTypeRelations.MainSequence()
PMS = SpectralTypeRelations.PreMainSequence(pms_tracks_file="%s/Dropbox/School/Research/Stellar_Evolution/Baraffe_Tracks.dat" %(os.environ["HOME"]), track_source="Baraffe")
PMS2 = SpectralTypeRelations.PreMainSequence()
def GetFluxRatio(sptlist, Tsec, xgrid, age=None):
  """
    Returns the flux ratio between the secondary star of temperature Tsec
    and the (possibly multiple) primary star(s) given in the 
    'sptlist' list (given as spectral types)
    xgrid is a numpy.ndarray containing the x-coordinates to find the 
      flux ratio at (in nm)

    The age of the system is found from the main-sequence age of the 
      earliest spectral type in sptlist, if it is not given
  """
  prim_flux = numpy.zeros(xgrid.size)
  sec_flux = numpy.zeros(xgrid.size)

  #First, get the age of the system
  if age is None:
    age = GetAge(sptlist)

  #Now, determine the flux from the primary star(s)
  for spt in sptlist:
    end = search("[0-9]", spt).end()
    T = MS.Interpolate(MS.Temperature, spt[:end])
    R = PMS2.GetFromTemperature(age, T, key="Radius")
    prim_flux += Planck(T, xgrid*units.nm.to(units.cm)) * R**2

  #Determine the secondary star flux
  R = PMS.GetFromTemperature(age, Tsec, key="Radius")
  sec_flux = Planck(Tsec, xgrid*units.nm.to(units.cm)) * R**2

  return sec_flux / prim_flux


def GetMass(spt, age):
  """
  Returns the mass of the system in solar masses
  spt: Spectral type of the star
  age: age, in years, of the system
  """

  #Get temperature
  end = search("[0-9]", spt).end()
  T = MS.Interpolate(MS.Temperature, spt[:end])

  # Determine which tracks to use
  if spt[0] == "O" or spt[0] == "A" or spt[0] == "F":
    return PMS2.GetFromTemperature(age, T, key="Mass")
  else:
    return PMS.GetFromTemperature(age, T, key="Mass")


def GetAge(sptlist):
  """
  Returns the age of the system, in years, given a list
  of spectral types. It determines the age as the 
  main sequence lifetime of the earliest-type star
  """
  lowidx = 999
  for spt in sptlist:
    if "I" in spt:
      #Pre-main sequence. Just use the age of an early O star, 
      # which is ~ 1Myr
      lowidx = 1
      break
    end = search("[0-9]", spt).end()
    idx = MS.SpT_To_Number(spt[:end])
    if idx < lowidx:
      lowidx = idx
  spt = MS.Number_To_SpT(lowidx)
  Tprim = MS.Interpolate(MS.Temperature, spt)
  age = PMS2.GetMainSequenceAge(Tprim, key="Temperature")
  return age

  





  


if __name__ == "__main__":
  lightspeed = constants.c.cgs.value * units.cm.to(units.km)
  smooth_factor = 0.8
  vel_list = range(-400, 400, 50)
  companion_file = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Multiplicity.csv" %(os.environ["HOME"])
  vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" %(os.environ["HOME"])
  fileList = []
  tolerance = 5.0
  for arg in sys.argv[1:]:
    if "-m" in arg:
      companion_file = arg.split("=")[1]
    elif "-tol" in arg:
      tolerance = float(arg.split("=")[1])
    else:
      fileList.append(arg)

  # Make sure each output file exists:
  logfilenames = {}
  for fname in fileList:
    output_dir = "Sensitivity/"
    outfilebase = fname.split(".fits")[0]
    if "/" in fname:
      dirs = fname.split("/")
      output_dir = ""
      outfilebase = dirs[-1].split(".fits")[0]
      for directory in dirs[:-1]:
        output_dir = output_dir + directory + "/"
      output_dir = output_dir + "Sensitivity/"
    HelperFunctions.ensure_dir(output_dir)
    
    #Make the summary file
    logfile = open(output_dir + "logfile.dat", "w")
    logfile.write("Sensitivity Analysis:\n*****************************\n\n")
    logfile.write("Filename\t\t\tPrimary Temperature\tSecondary Temperature\tMass (Msun)\tMass Ratio\tVelocity\tPeak Correct?\tSignificance\n")
    logfile.close()
    logfilenames[fname] = output_dir + "logfile.dat"

  # Read in the companion file
  companions = ascii.read(companion_file)[20:]

  # Read in the vsini file
  vsini_data = ascii.read(vsini_file)[10:]

  # Now, start loop over the models:
  start = 0
  end = len(model_list)
  fitorders = 8*numpy.ones(len(model_list))
  fitorders = [7, 6, 8, 8, 8, 7, 7, 10, 9, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 8, 8, 8, 8]
  """
  fitorders[0] = 7
  fitorders[1] = 6
  fitorders[5] = 7
  fitorders[6] = 7
  fitorders[7] = 10
  fitorders[8] = 9
  fitorders[9] = 7
  fitorders[10] = 7
  fitorders[11] = 7
  fitorders[12] = 7
  fitorders[13] = 7
  fitorders[14] = 7
  fitorders[15] = 7
  fitorders[16] = 6
  fitorders[17] = 6
  fitorders[18] = 6
  fitorders[19] = 6
  """
  #start = 15
  #end = 16
  for i, modelfile in enumerate(model_list[start:end]):
    if "PHOENIX2004" in modelfile:
      temp = int(modelfile.split("lte")[-1][:2])*100
      gravity = float(modelfile.split("lte")[-1][3:6])
      metallicity = float(modelfile.split("lte")[-1][6:10])
    elif "PHOENIX-ACES" in modelfile:
      temp = int(modelfile.split("lte")[-1][:2])*100
      gravity = float(modelfile.split("lte")[-1][3:7])
      metallicity = float(modelfile.split("lte")[-1][7:11])
    print "Reading in file %s" %modelfile
    x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
    c = FittingUtilities.Continuum(x, y, fitorder=fitorders[i+start], lowreject=2, highreject=7)
    plt.plot(x,y)
    plt.plot(x,c)
    #for order in range(6,11):
    #  print order
    #  c = FittingUtilities.Continuum(x, y, fitorder=order, lowreject=2, highreject=7, numiter=10)
    #  plt.plot(x,c, label="%ith order" %order)
    plt.title(modelfile.split("/")[-1])
    #plt.legend(loc='best')
    plt.title(str(temp))
    plt.show()
    continue
    model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y, cont=10**c)
    model = FittingUtilities.RebinData(model, numpy.linspace(model.x[0], model.x[-1], model.size()))
    model = Broaden.RotBroad(model, vsini)
    model = Broaden.ReduceResolution2(model, resolution)
    modelfcn = interp(model.x, model.y/model.cont)
    

    # Now that we have a spline function for the broadened data,
    # begin looping over the files
    for fname in fileList:
      outfile = open(logfilenames[fname], "a")

      # Read in and process the data like I am about to look for a companion
      orders_original = Process_Data(fname)

      #Find the vsini of the primary star with my spreadsheet
      starname = pyfits.getheader(fname)["object"]
      found = False
      for data in vsini_data:
        if data[0] == starname:
          vsini = float(data[1])
          found = True
      if not found:
        sys.exit("Cannot find %s in the vsini data: %s" %(starname, vsini_file))
      print starname, vsini

      #Check for companions in my master spreadsheet
      known_stars =  []
      if starname in companions.field(0):
        row = companions[companions.field(0) == starname]
        known_stars.append(row['col1'].item())
        ncompanions = int(row['col4'].item())
        for comp in range(ncompanions):
          spt = row["col%i" %(6 + 4*comp)].item()
          if not "?" in spt and (spt[0] == "O" or spt[0] == "B" or spt[0] == "A" or spt[0] == "F"):
            sep = row["col%i" %(7 + 4*comp)].item()
            if (not "?" in sep) and float(sep) < 4.0:
              known_stars.append(spt)
      else:
        sys.exit("Star not found in multiplicity library!")

      #Determine the age of the system and properties of the primary and secondary star
      age = GetAge(known_stars)
      primary_spt = known_stars[0]
      print primary_spt
      end = search("[0-9]", primary_spt).end()
      primary_temp = MS.Interpolate(MS.Temperature, primary_spt[:end])
      primary_mass = GetMass(primary_spt, age)
      secondary_spt = MS.GetSpectralType(MS.Temperature, temp)
      secondary_mass = GetMass(secondary_spt, age)
      massratio = secondary_mass / primary_mass


      for rv in vel_list:
        orders = [order.copy() for order in orders_original]   #Make a copy of orders
        model_orders = []
        for ordernum, order in enumerate(orders):
          #Ensure the x-axis
          #xgrid = numpy.linspace(order.x[0], order.x[-1], order.size())
          #order = FittingUtilities.RebinData(order, xgrid)
          scale = GetFluxRatio(known_stars, temp, order.x)
          print ordernum, scale.mean()

          #Add the model to the data
          model = (modelfcn(order.x*(1.0+rv/lightspeed)) - 1.0) * scale
          order.y += model*order.cont

          #Smooth data using the vsini of the primary star
          dx = order.x[1] - order.x[0]
          npixels = Smooth.roundodd(vsini/lightspeed * order.x.mean()/dx * smooth_factor)
          smoothed = Smooth.SmoothData(order, 
                                       windowsize=npixels, 
                                       smoothorder=3, 
                                       lowreject=3, 
                                       highreject=3,
                                       expand=10, 
                                       numiters=10)
          order.y /= smoothed.y

          # log-space the data
          start = numpy.log(order.x[0])
          end = numpy.log(order.x[-1])
          neworder = order.copy()
          xgrid = numpy.logspace(start, end, order.size(), base=numpy.e)
          logspacing = numpy.log(xgrid[1]/xgrid[0])
          order = FittingUtilities.RebinData(order, xgrid)

          # Generate a model with the same log-spacing (and no rv shift)
          dlambda = order.x[order.size()/2] * 1000*1.5/lightspeed
          start = numpy.log(order.x[0] - dlambda)
          end = numpy.log(order.x[-1] + dlambda)
          xgrid = numpy.exp(numpy.arange(start, end+logspacing, logspacing))
          model = DataStructures.xypoint(x=xgrid, cont=numpy.ones(xgrid.size))
          model.y = modelfcn(xgrid)

          # Save model order
          model_orders.append(model)

        #Do the actual cross-correlation
        corr = Correlate.Correlate(orders, model_orders)#, debug=True, outputdir="Sensitivity_Testing/")

        plt.figure(int(rv/50+ 10))
        plt.plot(corr.x, corr.y)
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("CCF")
        plt.title("RV = %g" %rv)

        # Check if we found the companion
        idx = numpy.argmax(corr.y)
        vmax = corr.x[idx] - 2.0 #There is a systematic offset for some reason!
        fit = FittingUtilities.Continuum(corr.x, corr.y, fitorder=2, lowreject=3, highreject=2.5)
        corr.y -= fit
        goodindices = numpy.where(numpy.abs(corr.x - rv) > 20)[0]
        mean = corr.y[goodindices].mean()
        std = corr.y[goodindices].std()
        significance = (corr.y[idx] - mean)/std
        if abs(vmax - rv) <= tolerance:
          #Found
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" %(fname, primary_temp, temp, secondary_mass, massratio, rv, significance) )
        else:
          #Not found
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\tN/A\n" %(fname, primary_temp, temp, secondary_mass, massratio, rv) )
        
      outfile.close()



  
  plt.show()
  

"""
if __name__ == "__main__":
  #Parse command line arguments:
  fileList = []
  extensions=True
  tellurics=False
  trimsize = 1
  windowsize = 101
  vsini = 10.0*units.km.to(units.cm)
  MS = SpectralTypeRelations.MainSequence()
  PMS = SpectralTypeRelations.PreMainSequence()
  vel_list = range(-400, 400, 50)
  outdir = "Sensitivity/"
  for arg in sys.argv[1:]:
    if "-e" in arg:
      extensions=False
    if "-t" in arg:
      tellurics=True  #telluric lines modeled but not removed
    else:
      fileList.append(arg)

  ensure_dir(outdir)
  outfile = open(outdir + "logfile.dat", "w")
  outfile.write("Sensitivity Analysis:\n*****************************\n\n")
  outfile.write("Filename\t\t\tPrimary Temperature\tSecondary Temperature\tMass (Msun)\tMass Ratio\tVelocity\tPeak Correct?\tSignificance\n")
  
  for fname in fileList:
    orders_original = ProcessData(fname)
    
    #Read in the name of the star from the fits header
    header = pyfits.getheader(fname)
    starname = header["OBJECT"]
    print starname

    #Get spectral type of the primary from the name and simbad
    stardata = StarData.GetData(starname)
    primary_temp = [ MS.Interpolate(MS.Temperature, stardata.spectype[:2]), ]
    #age = 'MS'   #Play with this later!
    primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])
    age = PMS.GetMainSequenceAge(primary_mass)
    

    #Check for close companions
    companions = HelperFunctions.CheckMultiplicityWDS(starname)
    if companions:
      for configuration in companions:
        component = companions[configuration]
        if component["Separation"] < 3.0 and component["Secondary SpT"] != "Unknown":
          print "Known %s companion with a separation of %g arcseconds!" %(component["Secondary SpT"], component["Separation"])
          primary_temp.append(MS.Interpolate(MS.Temperature, component["Secondary SpT"]))
    companion = HelperFunctions.CheckMultiplicitySB9(starname)
    if companion:
      if companion['K1'] != 'Unknown' and companion['K2'] != 'Unknown':
        mass = primary_mass * companion['K1'] / companion['K2']
        spt = MS.GetSpectralType(MS.Mass, mass, interpolate=True)
        primary_temp.append(MS.Interpolate(MS.Temperature, spt))
          

    #Begin loop over model spectra
    for j, model in enumerate(model_data):
            
      #Get info about the secondary star for this model temperature
      secondary_spt = MS.GetSpectralType(MS.Temperature, temp_list[j])
      secondary_radius = PMS.Interpolate(secondary_spt, age, key='Radius')
      secondary_mass = PMS.Interpolate(secondary_spt, age, key='Mass')
      massratio = secondary_mass / primary_mass

      #Rotationally Broaden model
      left = numpy.searchsorted(model.x, orders_original[0].x[0] - 10.0)
      right = numpy.searchsorted(model.x, orders_original[-1].x[-1] + 10.0)
      model = Broaden.RotBroad(model[left:right], vsini, linear=True)

      #Reduce resolution
      model = FittingUtilities.ReduceResolution2(model, 60000, extend=False)

      #Check sensitivity to this star
      orders = [order.copy() for order in orders_original]   #Make a copy of orders
      output_dir = "%s/Sensitivity/" %(os.path.dirname(fname))
      if output_dir.startswith("/"):
        output_dir = output_dir[1:]
      found, significance = Sensitivity.Analyze(orders, 
                                                model, 
                                                prim_temp=primary_temp, 
                                                sec_temp=temp_list[j],
                                                age=age,
                                                smoothing_windowsize=101,
                                                smoothing_order=5,
                                                outdir=output_dir,
                                                outfilebase=fname.split(".fits")[0],
                                                debug=False)

      #Write to logfile
      vels=range(-400, 450, 50)
      for v,f,s in zip(vels, found, significance):
        if f:
          #Signal found!
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" %(fname, primary_temp[0], temp_list[j], secondary_mass, massratio, v, s) )
        else:
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\tN/A\n" %(fname, primary_temp[0], temp_list[j], secondary_mass, massratio, v) )
      
       

          
  outfile.close()
"""      
    


