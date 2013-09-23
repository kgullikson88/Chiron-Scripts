"""
  This script goes through the stars observed, and searches for both known and new
  companions to each target. Right now, it only does known companions automatically
"""

import pySIMBAD as sim
import sys
import os
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
import matplotlib.pyplot as plt
import numpy
import HelperFunctions
import SpectralTypeRelations
import pyfits
import StarData
from astropy import units, constants


"""
  The following dictionary shows the new detections from my data.
  It includes both brand-new detections, and detections of the secondary
    star lines in known SB1s.
  The key is the star name, and the value is the estimated temperature
"""
NewDetections = {"HIP 92855": 0.1} #TODO: actually make this a temperature! Is a mass ratio now for testing!



"""
  This function will search the WDS catalog for known companions within 'sep' arcseconds
"""
def GetWDSCompanions(starname, sep=5.0, MS=None):
  if MS == None:
    MS = SpectralTypeRelations.MainSequence()
  companions = HelperFunctions.CheckMultiplicityWDS(starname)
  companion_info = []
  if companions:
    for configuration in companions:
      component = companions[configuration]
      if component["Separation"] < sep:
        s_spt = component["Secondary SpT"]
        if s_spt == "Unknown":
          print "Warning! Star %s has a companion with unknown magnitude/spectral type in WDS!" %starname
        else:
          mass = MS.Interpolate(MS.Mass, s_spt)
          companion_info.append((component["Separation"], mass))
  return companion_info


"""
   This function searches the SB9 catalog for spectroscopic companions
   The return type is determined by what information is given in the database,
     but always consists of an integer followed by a float

     -If no companion exists, the return values are 0,0
     -If both K1 and K2 are known (double-lined binary):
        -integer returned is 1, float returned is the mass ratio
     -If only K1 is known (the case for most):
        -integer returned is 2, and the float is the mass function f(M2)=M2 (sin(i))^2 / (1+1/q)^2
"""
def GetSB9Companions(starname, MS=None):
  if MS == None:
    MS = SpectralTypeRelations.MainSequence()
  companion = HelperFunctions.CheckMultiplicitySB9(starname)
  if not companion:
    return 0, 0
  if companion["K1"] != "Unknown" and companion["K2"] != "Unknown":
    q = companion["K1"] / companion["K2"]
    return 1, q
  elif companion["K1"] != "Unknown":
    K1 = companion["K1"]
    P = companion["Period"]
    mass_fcn = (K1*units.km.to(units.cm))**3 * (P*units.day.to(units.second)) / (2*numpy.pi*constants.G.cgs.value)
    return 2, mass_fcn * units.gram.to(units.solMass)


if __name__ == "__main__":
  dirlist = []
  for arg in sys.argv[1:]:
    dirlist.append(arg)
  if len(dirlist) == 0:
    dirlist = [d for d in os.listdir("./") if d.startswith("2013")]

  MS = SpectralTypeRelations.MainSequence()

  multiplicity = 0.0
  numstars = 0.0
  mass_ratios = []
  new_massratios = []
  for directory in dirlist:
    starlist = [f for f in os.listdir(directory) if f.startswith("H") and f.endswith("-0.fits")]
    for star in starlist:
      multiple = False
      #First, get the known companions
      multiple = False
      header = pyfits.getheader("%s/%s" %(directory, star))
      starname = header['OBJECT']
      print starname
      stardata = StarData.GetData(starname)
      primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])
      known_companions = GetWDSCompanions(starname, MS=MS, sep=100.0)
      code, value = GetSB9Companions(starname)
      if len(known_companions) > 0:
        multiple = True
        for comp in known_companions:
          print "\tq = %g" %(comp[1]/(primary_mass))
          mass_ratios.append(comp[1]/primary_mass)
      if code == 1:
        multiple = True
        q = value
        wds = False
        for comp in known_companions:
          if abs(q-comp[1]) < 0.1 and comp[0] < 4.0:
            wds = True
        if not wds:
          mass_ratios.append(q)
        else:
          print "Spectroscopic binary found which may match a WDS companion."
          usr = raw_input("Use both (y or n)? ")
          if "y" in usr:
            mass_ratios.append(q)
        print "Spectroscopic companion with q = %g" %q
      elif code == 2:
        print "Single-lined spectroscopic companion to %s found! Double-lined in my data?" %starname
        multiple = True


      #Now, put in my data
      if starname in NewDetections:
        new_massratios.append(NewDetections[starname])
        multiple = True

      #Keep track of total binary fraction
      if multiple:
        multiplicity += 1
      numstars += 1.0


  #Make some plots
  mass_ratios = [min(q, 1.0) for q in mass_ratios]
  print "Multiplicity fraction = %g" %(multiplicity/numstars)
  bins = numpy.arange(0.0, 1.1, 0.1)
  print bins.size, '\t', bins
  plt.figure(1)
  if len(new_massratios) > 0:
    print "Found new entries!"
    mass_ratios = [mass_ratios, new_massratios]
  print len(mass_ratios)
  plt.hist(mass_ratios, bins=bins, color=['gray','green'], stacked=True)
  """
  if len(new_massratios) > 0:
    y,edges = numpy.histogram(new_massratios, bins=bins)
    print y
    print edges
    plt.bar(bins[:-1], y, bottom=numpy.array(height), color='green', align='edge')
    #plt.hist(new_massratios, bins=bins, bottom=height, color='green')
  """
  plt.xlabel(r"$\rm M_s/M_p$")
  plt.ylabel("Frequency")
  plt.title("Mass Ratio Distribution for Companions within 100\"")
  
  plt.figure(2)
  plt.hist(mass_ratios, bins=bins, color=['gray','green'], cumulative=True, normed=True, histtype='step', linewidth=2, stacked=True)
  plt.plot(bins, bins, 'k--', linewidth=2)
  plt.xlabel(r"$\rm M_s/M_p$")
  plt.ylabel("Cumulative Frequency")
  plt.show()
      
