import os

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np

import HelperFunctions


def check_all():
    filenames = [f for f in os.listdir("./") if f.endswith("smoothed.fits") and f.startswith("H")]
    corrdir = "Cross_correlations/"
    vsini_values = [1, 10, 20, 30, 40]
    Temperatures = [3300, 3500, 3700, 3900, 4200, 4500, 5000, 5500]
    Temperatures = range(3000, 6800, 100)
    metals = [-0.5, 0.0, 0.5]
    logg = 4.5
    HelperFunctions.ensure_dir("Figures/")

    for rootfile in sorted(filenames):
        Tvals = []
        Zvals = []
        rotvals = []
        significance = []
        corrval = []
        for T in Temperatures:
            for metal in metals:
                for vsini in vsini_values:
                    corrfile = "{0:s}{1:s}.{2:d}kps_{3:.1f}K{4:+.1f}{5:+.1f}".format(corrdir,
                                                                                   rootfile.split(".fits")[0],
                                                                                   vsini,
                                                                                   T,
                                                                                   logg,
                                                                                   metal)
                    print corrfile
                    try:
                        vel, corr = np.loadtxt(corrfile, unpack=True)
                    except IOError:
                        continue

                    # Check the significance of the highest peak within +/- 500 km/s
                    left = np.searchsorted(vel, -500)
                    right = np.searchsorted(vel, 500)
                    idx = np.argmax(corr[left:right]) + left
                    v = vel[idx]
                    goodindices = np.where(np.abs(vel - v) > vsini)[0]
                    std = np.std(corr[goodindices])
                    mean = np.mean(corr[goodindices])
                    mean = np.median(corr)
                    mad = HelperFunctions.mad(corr)
                    std = 1.4826 * mad
                    sigma = (corr[idx] - mean) / std

                    """
                    # Plot if > 3 sigma peak
                    if sigma > 4:
                        fig = plt.figure(10)
                        ax = fig.add_subplot(111)
                        ax.plot(vel, corr, 'k-', lw=2)
                        ax.set_xlabel("Velocity (km/s)")
                        ax.set_ylabel("CCF")
                        ax.set_title(r'{0:s}:  $T_s$={1:d}K & [Fe/H]={2:.1f}'.format(rootfile, T, metal))
                        ax.grid(True)
                        fig.savefig(u"Figures/{0:s}.pdf".format(corrfile.split("/")[-1]))
                        plt.close(fig)
                    """

                    Tvals.append(T)
                    Zvals.append(metal)
                    rotvals.append(vsini)
                    significance.append(sigma)
                    corrval.append(corr[idx] - np.median(corr))

        # Now, make a plot of the significance as a function of Temperature and metallicity for each vsini
        Tvals = np.array(Tvals)
        Zvals = np.array(Zvals)
        rotvals = np.array(rotvals)
        significance = np.array(significance)
        corrval = np.array(corrval)
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Significance Summary for %s" % (rootfile.split(".fits")[0].replace("_", " ")))
        for i, rot in enumerate(vsini_values):
            goodindices = np.where(abs(rotvals - rot) < 1e-5)[0]
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("[Fe/H]")
            ax.set_zlabel("Significance")
            ax.plot(Tvals[goodindices], Zvals[goodindices], significance[goodindices], 'o', label="%i km/s" % rot)
            #ax.plot(Tvals[goodindices], Zvals[goodindices], corrval[goodindices], 'o', label="{0:d} km/s".format(rot))
        leg = ax.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        fig.savefig("Figures/Summary_{0:s}.pdf".format(rootfile.split(".fits")[0]))
        idx = np.argmax(significance)
        #ax.plot(Tvals[idx], Zvals[idx], significance[idx], 'x', markersize=25, label="Most Significant")
        print os.getcwd()
        plt.show()




def check_high_t(T=6000, metal=0.0, vsini=10):
    filenames = [f for f in os.listdir("./") if f.endswith("smoothed.fits") and f.startswith("H")]
    corrdir = "Cross_correlations/"
    logg = 4.5
    HelperFunctions.ensure_dir("Figures/")

    for rootfile in sorted(filenames):
        corrfile = "{0:s}{1:s}.{2:d}kps_{3:.1f}K{4:+.1f}{5:+.1f}".format(corrdir,
                                                                         rootfile.split(".fits")[0],
                                                                         vsini,
                                                                         T,
                                                                         logg,
                                                                         metal)
        print corrfile
        try:
            vel, corr = np.loadtxt(corrfile, unpack=True)
        except IOError:
            continue

        plt.plot(vel, corr, 'k-')
        plt.xlabel("Velocity")
        plt.ylabel("CCF")
        plt.title(rootfile.split(".fits")[0])
        plt.show()



if __name__ == "__main__":
    check_all()
