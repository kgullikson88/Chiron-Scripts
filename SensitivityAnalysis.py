#import Correlate
import FitsUtils
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import scipy.signal
import os
import sys
import DataStructures
#import FindContinuum
import matplotlib.pyplot as plt
import pyfits
from astropy import units, constants
import StarData
import SpectralTypeRelations
from PlotBlackbodies import Planck
import RotBroad
import FittingUtilities
import MakeModel


#Ensure a directory exists. Create it if not
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        

homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[588.8, 589.9],
              [627.1, 635.4]]
badregions = [[0, 466],
#badregions = [[0, 540],
              [567.5, 575.5],
              [587.5, 593],
              [627, 634.5],
              [686, 706],
              [716, 742],
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
   
               
star_list = []
temp_list = []
gravity_list = []
metal_list = []
model_data = []
for fname in model_list:
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
  model_data.append( DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y) )
  star_list.append(str(temp))
  temp_list.append(temp)
  gravity_list.append(gravity)
  metal_list.append(metallicity)


  

if __name__ == "__main__":
  #Parse command line arguments:
  fileList = []
  extensions=True
  tellurics=False
  trimsize = 100
  windowsize = 91
  MS = SpectralTypeRelations.MainSequence()
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
    if extensions:
      orders_original = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="flux", errors="error")
      if tellurics:
        model_orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="model")
        for i, order in enumerate(orders_original):
          orders_original[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
          orders_original[i].y /= model_orders[i].y
          
    else:
      orders_original = FitsUtils.MakeXYpoints(fname, errors=2)

    #Loop over orders, removing bad parts
    numorders = len(orders_original)
    for i, order in enumerate(orders_original[::-1]):
      #Only use the middle half of each order (lots of noise on the edges)
      DATA = interp(order.x, order.y)
      CONT = interp(order.x, order.cont)
      ERROR = interp(order.x, order.err)
      left = int(order.size()/4.0)
      right = int(order.size()*3.0/4.0 + 0.5)
      order.x = numpy.linspace(order.x[left], order.x[right],right - left + 1)
      order.y = DATA(order.x)
      order.cont = CONT(order.x)
      order.err = ERROR(order.x)

      
      #Remove bad regions from the data
      for region in badregions:
        left = numpy.searchsorted(order.x, region[0])
        right = numpy.searchsorted(order.x, region[1])
        order.x = numpy.delete(order.x, numpy.arange(left, right))
        order.y = numpy.delete(order.y, numpy.arange(left, right))
        order.cont = numpy.delete(order.cont, numpy.arange(left, right))
        order.err = numpy.delete(order.err, numpy.arange(left, right))

      #Remove whole order if it is too small
      if order.x.size > windowsize:
        order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
        orders_original[numorders -1 -i] = order.copy()
      else:
        print "Removing order %i" %(numorders - 1 - i)
        orders_original.pop(numorders - 1 - i)
        

    #Read in the name of the star from the fits header
    header = pyfits.getheader(fname)
    starname = header["OBJECT"].split()[0].replace("_", " ")

    #Get spectral type of the primary from the name and simbad
    stardata = StarData.GetData(starname)
    primary_temp = MS.Interpolate(MS.Temperature, stardata.spectype[:2])
    primary_radius = MS.Interpolate(MS.Radius, stardata.spectype[:2])
    primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])

    #Begin loop over model spectra
    for j, model in enumerate(model_data):
      orders = [order.copy() for order in orders_original]   #Make a copy of orders
      #Get info about the secondary star for this model temperature
      secondary_spt = MS.GetSpectralType(MS.Temperature, temp_list[j])
      secondary_radius = MS.Interpolate(MS.Radius, secondary_spt)
      secondary_mass = MS.Interpolate(MS.Mass, secondary_spt)
      massratio = secondary_mass / primary_mass

      left = numpy.searchsorted(model.x, orders[0].x[0] - 10.0)
      right = numpy.searchsorted(model.x, orders[-1].x[-1] + 10.0)
      model = RotBroad.Broaden2(model[left:right], 20*units.km.to(units.cm))
      MODEL = interp(model.x, model.y)

      #Loop over velocities
      for vel in vel_list:
        corrlist = []
        normalization = 0.0
        for i, order in enumerate(orders[::-1]):
          order2 = order.copy()
          #Process the model
          #a: make a segment of the total model to work with
          left = max(0, numpy.searchsorted(model.x, order2.x[0] - 10)-1 )
          right = min(model.size()-1, numpy.searchsorted(model.x, order2.x[-1] + 10))
      
          model2 = DataStructures.xypoint(right - left + 1)
          model2.x = numpy.linspace(model.x[left], model.x[right], right - left + 1)
          model2.y = MODEL(model2.x*(1.0+vel/3e5))
          model3 = model2.copy()
          model3.y = MODEL(model2.x)
          model2.cont = FittingUtilities.Continuum(model2.x, model2.y, fitorder=3, lowreject=1.5, highreject=10.0)
          model3.cont = FittingUtilities.Continuum(model3.x, model3.y, fitorder=3, lowreject=1.5, highreject=10.0)

          #b: Convolve to detector resolution
          model2 = MakeModel.ReduceResolution(model2.copy(), 60000, extend=False)
          model3 = MakeModel.ReduceResolution(model3.copy(), 60000, extend=False)

          #c: rebin to the same spacing as the data
          xgrid = numpy.arange(model2.x[0], model2.x[-1], order2.x[1] - order2.x[0])
          model2 = MakeModel.RebinData(model2.copy(), xgrid)
          model3 = MakeModel.RebinData(model3.copy(), xgrid)

          #d: scale to be at the appropriate flux ratio
          primary_flux = Planck(order2.x.mean()*units.nm.to(units.cm), primary_temp)
          secondary_flux = Planck(order2.x.mean()*units.nm.to(units.cm), temp_list[j])
          scale = secondary_flux / primary_flux * (secondary_radius/primary_radius)**2
          model2.y = (model2.y/model2.cont - 1.0)*scale + 1.0
          model_fcn = interp(model2.x, model2.y)
          order2.y = (order2.y/order2.cont + model_fcn(order2.x))

          #Smooth data in the same way I would normally
          smoothed =  FittingUtilities.savitzky_golay(order2.y, windowsize, 5)
          reduceddata = order2.y/smoothed
          #vsini = 60.0
          #order2.x, order2.y = FittingUtilities.HighPassFilter(order2, vsini*units.km.to(units.cm), linearize=True)
          #x, reduceddata = FittingUtilities.HighPassFilter(order2, vsini*units.km.to(units.cm), linearize=True)
          #filterfcn = interp(x, reduceddata)
          #reduceddata = filterfcn(order2.x)
          #plt.plot(order2.x, order2.y)
          #plt.plot(order2.x, reduceddata+2)
          #plt.show()
          #model3.x, model3.y = FittingUtilities.HighPassFilter(model3, vsini*units.km.to(units.cm), linearize=True)
          #x, reducedmodel = FittingUtilities.HighPassFilter(model3, vsini*units.km.to(units.cm), linearize=True)
          #plt.plot(model3.x, model3.y/model3.cont)
          #plt.plot(x, reducedmodel/model3.cont+1)
          #plt.show()

          #Do the cross-correlations
          reducedmodel = model3.y
          #reduceddata = order2.y
          reducedmodel = model3.y/model3.cont
          meandata = reduceddata.mean()
          meanmodel = reducedmodel.mean()
          data_rms = numpy.sqrt(numpy.sum((reduceddata - meandata)**2))
          model_rms = numpy.sqrt(numpy.sum((reducedmodel - meanmodel)**2))
          left = numpy.searchsorted(model2.x, order2.x[0])
          right = model2.x.size - numpy.searchsorted(model2.x, order2.x[-1])
          delta = left - right

          #plt.plot(order2.x, reduceddata - meandata)
          #plt.plot(model2.x, reducedmodel - meanmodel)
          #plt.show()
          ycorr = scipy.signal.fftconvolve(reduceddata - meandata, (reducedmodel - meanmodel)[::-1], mode='valid')
          xcorr = numpy.arange(ycorr.size)
          lags = xcorr - (model2.x.size + order2.x.size + delta - 1.0)/2.0
          lags = xcorr - right
          distancePerLag = model2.x[1] - model2.x[0]
          offsets = -lags*distancePerLag
          velocity = offsets*3e5 / numpy.median(order2.x)
          velocity, ycorr = velocity[::-1], ycorr[::-1]
          left = numpy.searchsorted(velocity, -1000)
          right = numpy.searchsorted(velocity, +1000)
          corr = DataStructures.xypoint(right - left + 1)
          corr.x = velocity[left:right]
          corr.y = ycorr[left:right]/(data_rms*model_rms) * scale
          corrlist.append(corr.copy())
          normalization += 1.0 * scale

        #Add up the individual CCFs
        master_corr = corrlist[0]
        for corr in corrlist[1:]:
          correlation = interp(corr.x, corr.y)
          master_corr.y += correlation(master_corr.x)
        master_corr.y /= normalization

        #output
        outfilename = "%s%s_t%i_v%i" %(outdir, fname.split(".fits")[0], temp_list[j], vel)
        print "Outputting CCF to %s" %outfilename
        numpy.savetxt(outfilename, numpy.transpose((master_corr.x, master_corr.y)), fmt="%.10g")

        #Write to logfile
        idx = numpy.argmax(master_corr.y)
        vmax = master_corr.x[idx]
        fit = FittingUtilities.Continuum(master_corr.x, master_corr.y, fitorder=2, lowreject=3, highreject=3)
        master_corr.y -= fit
        mean = master_corr.y.mean()
        std = master_corr.y.std()
        significance = (master_corr.y[idx] - mean)/std
        tolerance = 10.0
        if numpy.abs(vmax - vel) <= tolerance:
          #Signal found!
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" %(fname, primary_temp, temp_list[j], secondary_mass, massratio, vel, significance) )
        else:
          outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\t%.2f\n" %(fname, primary_temp, temp_list[j], secondary_mass, massratio, vel, significance) )


          
  outfile.close()
      
    


