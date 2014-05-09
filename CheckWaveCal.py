"""
  This script loops through all directories and finds the co-added spectra files.
  It then pulls the x-axis of each order, and plots the x-axis in a different 
  figure for each order, and with different lines for each file. It is used
  as a check to see how stable the CHIRON instrument wavelength calibration is,
  and if I can safely use the faster search method on all data at once.
"""

import HelperFunctions
import os
import sys
import matplotlib.pyplot as plt
import astropy.time



if __name__ == "__main__":
  orders = []
  dirs = [d for d in os.listdir("./") if d.startswith("201") and len(d) == 8]
  for d in dirs:
    object_files = [f for f in os.listdir(d) if f.startswith("H") 
                                             and len(f.split("_")) == 2 
                                             and f.endswith(".fits") 
                                             and "-" not in f]
    #object_files = [f for f in os.listdir(d) if f.startswith("achi") and f.endswith("-0.fits")]
    print d
    for f in object_files:
      data = HelperFunctions.ReadExtensionFits("%s/%s" %(d,f))
      print "\t", f, len(data)
      if len(orders) < 1:
        for order in data:
          orders.append({d: order.x})
      else:
        for i, order in enumerate(data):
          orders[i][d] = order.x

  """
  start, end = 999, 1001
  for i, order in enumerate(orders):
    for date in sorted(order.keys()):
      plt.plot(order[date][start:end+1], label=date)
      plt.text((end-start)/2, order[date][(start+end)/2], date)
    plt.title("Order %i" %(i+1))
    plt.xlabel("Pixel number")
    plt.ylabel("Wavelength (nm)")
    #plt.legend(loc='best')
    plt.show()
  """

  pixel = 1000
  for i, order in enumerate(orders):
    for date in sorted(order.keys()):
      date2 = "%s-%s-%s" %(date[:4], date[4:6], date[6:])
      jd = astropy.time.Time(date2, scale='utc', format='iso').jd
      dx = order[date][pixel+1] - order[date][pixel]
      plt.plot(jd, order[date][pixel], 'ro')
    plt.xlabel("Julian Date")
    #plt.ylabel("Delta - Wavelength at pixel %i (nm)" %(pixel))
    plt.ylabel("Wavelength at pixel %i (nm)" %(pixel))
    ax = plt.gca()
    ax.ticklabel_format(style = 'sci', useOffset=False)
    plt.show()
