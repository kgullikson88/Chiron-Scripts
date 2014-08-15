import matplotlib.pyplot as plt
import numpy
import os
import HelperFunctions



if __name__ == "__main__":
  filenames = [f for f in os.listdir("./") if f.endswith("smoothed.fits") and f.startswith("H")]
  corrdir = "Cross_correlations/"
  vsini = "20"
  Temperatures = [3300, 3500, 3700, 3900, 4200, 4500, 5000, 5500]
  metals = ["+0.5", "-0.5"]
  logg = "+4.5"
  HelperFunctions.ensure_dir("Figures/")
  for rootfile in filenames:
    for T in Temperatures[:-1]:
      for metal in metals:
        #if T == 5000 and metal == "-0.5":
        #  logg = "+3.5"
        #else:
        #  logg = "+4.0"
        corrfile = "%s%s.%skps_%iK%s%s" %(corrdir, 
                                          rootfile.split(".fits")[0],
                                          vsini,
                                          T,
                                          logg,
                                          metal)
        print rootfile, T, metal
        vel, corr = numpy.loadtxt(corrfile, unpack=True)
        plt.plot(vel, corr, 'k-', lw=2)
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("CCF")
        plt.title(r"%s:  $T_s$=%iK & [Fe/H]=%s" %(rootfile, T, metal))
        plt.savefig("Figures/%s.pdf" %corrfile.split("/")[-1])
        plt.grid()
        plt.show()
        plt.clf()
