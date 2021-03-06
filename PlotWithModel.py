# import Correlate
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import scipy.signal
import os
import sys

import numpy as np
import DataStructures

import FitsUtils

#import FindContinuum
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
import StarData
import SpectralTypeRelations
from PlotBlackbodies import Planck
import RotBroad_Fast as RotBroad
import FittingUtilities
import MakeModel
import Smooth
import HelperFunctions
import Correlate


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
              [749.1, 749.45],
              [759, 9e9]]
badpixels_by_order = {}  #Empty for now

#Set up model list
model_list = [modeldir + "lte30-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
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
              modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
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
              modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]

"""
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
  x,y = np.loadtxt(fname, usecols=(0,1), unpack=True)
  model_data.append( DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y) )
  star_list.append(str(temp))
  temp_list.append(temp)
  gravity_list.append(gravity)
  metal_list.append(metallicity)
"""

if __name__ == "__main__":
    #Parse command line arguments:
    fileList = []
    vsini = 20.0 * units.km.to(units.cm)
    MS = SpectralTypeRelations.MainSequence()
    vel = 0.0
    T = 6000
    for arg in sys.argv[1:]:
        if "-vsini" in arg:
            vsini = float(arg.split("=")[-1]) * units.km.to(units.cm)
        elif "-rv" in arg:
            vel = float(arg.split("=")[-1])
        elif "-temp" in arg:
            T = float(arg.split("=")[-1])
        else:
            fileList.append(arg)

    for fname in fileList:
        orders_original = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")

        #Loop over orders, removing bad parts
        numorders = len(orders_original)
        for i, order in enumerate(orders_original[::-1]):
            #Linearize
            DATA = interp(order.x, order.y)
            CONT = interp(order.x, order.cont)
            ERROR = interp(order.x, order.err)
            left, right = 20, -20
            order.x = np.linspace(order.x[left], order.x[right], order.size())
            order.y = DATA(order.x)
            order.cont = CONT(order.x)
            order.err = ERROR(order.x)

            #Interpolate over bad pixels
            if i in badpixels_by_order:
                left, right = badpixels_by_order[i][0], badpixels_by_order[i][1]
                m = (order.y[right] - order.y[left]) / (order.x[right] - order.x[left])
                order.y[left:right] = m * (order.x[left:right] - order.x[left]) + order.y[left]
                order.err[left:right] = 9e9

            #Remove bad regions from the data
            for region in badregions:
                left = np.searchsorted(order.x, region[0])
                right = np.searchsorted(order.x, region[1])
                if left == 0 or right == order.size():
                    order.x = np.delete(order.x, np.arange(left, right))
                    order.y = np.delete(order.y, np.arange(left, right))
                    order.cont = np.delete(order.cont, np.arange(left, right))
                    order.err = np.delete(order.err, np.arange(left, right))
                else:
                    print "Warning! Bad region covers the middle of order %i" % i
                    print "Interpolating rather than removing"
                    order.y[left:right] = order.cont[left:right]
                    order.err[left:right] = 9e9

            #Remove whole order if it is too small
            remove = False
            if order.x.size <= 10:
                remove = True
            else:
                velrange = 3e5 * (np.median(order.x) - order.x[0]) / np.median(order.x)
                if velrange <= 1000.0:
                    remove = True
            if remove:
                print "Removing order %i" % (numorders - 1 - i)
                orders_original.pop(numorders - 1 - i)
            else:
                order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
                orders_original[numorders - 1 - i] = order.copy()



        #Read in the name of the star from the fits header
        header = pyfits.getheader(fname)
        starname = header["OBJECT"]

        #Get spectral type of the primary from the name and simbad
        stardata = StarData.GetData(starname)
        primary_temp = MS.Interpolate(MS.Temperature, stardata.spectype[:2])
        primary_radius = MS.Interpolate(MS.Radius, stardata.spectype[:2])
        primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])
        companions = HelperFunctions.CheckMultiplicityWDS(starname)

        #Find the appropriate spectrum
        bestindex = 0
        besttemp = 0.0
        for i in range(len(model_list)):
            fname = model_list[i]
            if "PHOENIX2004" in fname:
                temp = float(fname.split("lte")[-1][:2]) * 100
            elif "PHOENIX-ACES" in fname:
                temp = float(fname.split("lte")[-1][:2]) * 100
            if abs(temp - besttemp) > abs(temp - T):
                besttemp = temp
                bestindex = i

        #Begin processing the model
        print "Best model is %s" % model_list[bestindex]
        x, y = np.loadtxt(model_list[bestindex], usecols=(0, 1), unpack=True)
        model = DataStructures.xypoint(x=x * units.angstrom.to(units.nm) / 1.00026, y=10 ** y)
        #Get info about the secondary star for this model temperature
        secondary_spt = MS.GetSpectralType(MS.Temperature, besttemp)
        secondary_radius = MS.Interpolate(MS.Radius, secondary_spt)
        secondary_mass = MS.Interpolate(MS.Mass, secondary_spt)
        massratio = secondary_mass / primary_mass

        #Rotationally Broaden model
        left = np.searchsorted(model.x, orders_original[0].x[0] - 10.0)
        right = np.searchsorted(model.x, orders_original[-1].x[-1] + 10.0)
        model = RotBroad.Broaden(model[left:right], vsini, linear=False)
        MODEL = interp(model.x, model.y)

        #Loop over orders, making a model for each
        orders = [order.copy() for order in orders_original]  #Make a copy of orders
        for i, order in enumerate(orders):
            #a: make a segment of the total model to work with
            left = max(0, np.searchsorted(model.x, order.x[0] - 10) - 1)
            right = min(model.size() - 1, np.searchsorted(model.x, order.x[-1] + 10))

            model2 = DataStructures.xypoint(right - left + 1)
            model2.x = np.linspace(model.x[left], model.x[right], right - left + 1)
            model2.y = MODEL(model2.x * (1.0 + vel / 3e5))
            model2.cont = FittingUtilities.Continuum(model2.x, model2.y, fitorder=3, lowreject=1.5, highreject=10.0)

            #b: Convolve to detector resolution
            model2 = FittingUtilities.ReduceResolution(model2.copy(), 60000, extend=False)

            #c: rebin to the same spacing as the data
            xgrid = order.x
            model2 = FittingUtilities.RebinData(model2.copy(), xgrid)

            #d: scale to be at the appropriate flux ratio
            primary_flux = Planck(order.x.mean() * units.nm.to(units.cm), primary_temp)
            #Check for known secondaries
            if companions:
                for configuration in companions:
                    component = companions[configuration]
                    if component["Separation"] < 3.0 and component["Secondary SpT"] != "Unknown":
                        if i == 0:
                            print "Known %s companion with a separation of %g arcseconds!" % (
                            component["Secondary SpT"], component["Separation"])
                        temperature = MS.Interpolate(MS.Temperature, component["Secondary SpT"])
                        primary_flux += Planck(order.x.mean() * units.nm.to(units.cm), temperature)
            secondary_flux = Planck(order.x.mean() * units.nm.to(units.cm), besttemp)
            scale = secondary_flux / primary_flux * (secondary_radius / primary_radius) ** 2
            model2.y = (model2.y / model2.cont - 1.0) * scale + 1.0

            #Plot
            plt.plot(order.x, order.y / order.cont, 'k-')
            plt.plot(model2.x, model2.y, 'r-')
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Normalized Flux")
            plt.show()

            """


            #Loop over velocities
            for vel in vel_list:
              corrlist = []
              normalization = 0.0
              orders = [order.copy() for order in orders_original]   #Make a copy of orders
              for i, order in enumerate(orders):
                order2 = order.copy()
                #Process the model
                #a: make a segment of the total model to work with
                left = max(0, np.searchsorted(model.x, order2.x[0] - 10)-1 )
                right = min(model.size()-1, np.searchsorted(model.x, order2.x[-1] + 10))

                model2 = DataStructures.xypoint(right - left + 1)
                model2.x = np.linspace(model.x[left], model.x[right], right - left + 1)
                model2.y = MODEL(model2.x*(1.0+vel/3e5))
                model2.cont = FittingUtilities.Continuum(model2.x, model2.y, fitorder=3, lowreject=1.5, highreject=10.0)

                #b: Convolve to detector resolution
                model2 = FittingUtilities.ReduceResolution(model2.copy(), 60000, extend=False)

                #c: rebin to the same spacing as the data
                xgrid = np.arange(model2.x[0], model2.x[-1], order2.x[1] - order2.x[0])
                model2 = FittingUtilities.RebinData(model2.copy(), xgrid)

                #d: scale to be at the appropriate flux ratio
                primary_flux = Planck(order2.x.mean()*units.nm.to(units.cm), primary_temp)
                #Check for known secondaries
                if companions:
                  for configuration in companions:
                    component = companions[configuration]
                    if component["Separation"] < 3.0 and component["Secondary SpT"] != "Unknown":
                      if i == 0:
                        print "Known %s companion with a separation of %g arcseconds!" %(component["Secondary SpT"], component["Separation"])
                      temperature = MS.Interpolate(MS.Temperature, component["Secondary SpT"])
                      primary_flux += Planck(order2.x.mean()*units.nm.to(units.cm), temperature)
                secondary_flux = Planck(order2.x.mean()*units.nm.to(units.cm), temp_list[j])
                scale = secondary_flux / primary_flux * (secondary_radius/primary_radius)**2
                model2.y = (model2.y/model2.cont - 1.0)*scale
                model2.cont = np.ones(model2.size())
                model_fcn = interp(model2.x, model2.y)
                order2.y = (order2.y/order2.cont + model_fcn(order2.x))*order2.cont

                #Smooth data in the same way I would normally
                smoothed = Smooth.SmoothData(order2, windowsize, 5)
                order2.y /= smoothed.y
                order2.cont = FittingUtilities.Continuum(order2.x, order2.y, fitorder=2)
                orders[i] = order2.copy()


              #Do the actual cross-correlation using PyCorr2 (order by order with appropriate weighting)
              corr = Correlate.PyCorr2(orders, resolution=80000, models=[model_data[j],], stars=[star_list[j],], temps=[temp_list[j],], gravities=[gravity_list[j],], metallicities=[metal_list[j],], vsini=0.0, debug=False, save_output=False)[0]

              #output
              outfilename = "%s%s_t%i_v%i" %(outdir, fname.split(".fits")[0], temp_list[j], vel)
              print "Outputting CCF to %s" %outfilename
              np.savetxt(outfilename, np.transpose((corr.x, corr.y)), fmt="%.10g")

              #Write to logfile
              idx = np.argmax(corr.y)
              vmax = corr.x[idx]
              fit = FittingUtilities.Continuum(corr.x, corr.y, fitorder=2, lowreject=3, highreject=3)
              corr.y -= fit
              mean = corr.y.mean()
              std = corr.y.std()
              significance = (corr.y[idx] - mean)/std
              tolerance = 10.0
              if np.abs(vmax - vel) <= tolerance:
                #Signal found!
                outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" %(fname, primary_temp, temp_list[j], secondary_mass, massratio, vel, significance) )
              else:
                outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\t%.2f\n" %(fname, primary_temp, temp_list[j], secondary_mass, massratio, vel, significance) )



        outfile.close()
            """
    


