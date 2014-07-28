import numpy as np
import HelperFunctions
import s4
import sys
import matplotlib.pyplot as plt
from astropy import units


def main2():
  vsini=250
  fname = "/Users/kgulliks/t15000g40a00p00_w300_700_r1.dat"
  x,y,c = np.loadtxt(fname, usecols=(0,1,2), unpack=True)
  import DataStructures
  import Broaden
  order = DataStructures.xypoint(x=x, y=y, cont=c)
  order = Broaden.RotBroad(order, vsini*units.km.to(units.cm))
  left = np.searchsorted(order.x, 4780)
  right = np.searchsorted(order.x, 4960)
  arr = order[left:right].toarray(norm=True)
  syn = s4.synthesis.Synplot(15000, 
                             4.0, 
                             idl="/Applications/itt/idl/bin/idl",
                             wstart=order.x[left],
                             wend=order.x[right],
                             observ=arr,
                             relative=1,
                             vrot=vsini,
                             rv=83)
  syn.run()
  syn.plot()
  plt.figure(2)
  #spec = syn.spectrum
  #plt.plot(spec[:,0], spec[:,1])
  #plt.plot(arr[:,0], arr[:,1])
  plt.show()


def main1():
  fname = "HIP_77635.fits"
  orders = HelperFunctions.ReadExtensionFits(fname)
  
  for order in orders:
    if order.x[0] < 492 and order.x[-1] > 492:
      break
  order.x *= units.nm.to(units.angstrom)
  syn = s4.synthesis.Synplot(15000, 
                             4.0, 
                             idl="/Applications/itt/idl/bin/idl",
                             wstart=order.x[0],
                             wend=order.x[-1],
                             observ=order.toarray(norm=True),
                             relative=1,
                             vrot=225,
                             rv=-31)
  syn.run()
  syn.plot()
  plt.show()

  s4.fitting.Fit("teff", 18000, 21000, 500,
                 "logg", 3.5, 4.5, 0.2,
                 "rv", -60, -10, 10,
                 idl="/Applications/itt/idl/bin/idl",
                 wstart=order.x[0],
                 wend=order.x[-1],
                 observ=order.toarray(norm=True),
                 relative=1,
                 vrot=225)



if __name__ == "__main__":
  main2()
