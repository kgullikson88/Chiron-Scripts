import os
import sys
from collections import defaultdict
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
import FittingUtilities
from astropy import units

import SpectralTypeRelations
import StarData
import PlotBlackbodies


"""
  Program to analyze the output of SensitivityAnalysis, and make some pretty plots!

  Command line arguments:
     -combine: will combine several output (say as generated by xgrid) NOT YET IMPLEMENTED
     -xaxis: specifies the variable to use as the x axis. Choices are as follows
         SecondarySpectralType
         SecondaryMass
         MassRatio
         DetectionRate
         AverageSignificance
         MagnitudeDifference
     -yaxis: specifies the variable to use for the y axis. Choices are the same as for -xaxis
     -infile: specifies the input filename (default is Sensitivity/summary.dat).
         If combine is True, the input filename should be a list of comma-separated 
         filenames
"""


def MakeSummaryFile(directory, prefix, outfilename="Sensitivity/logfile.dat", tolerance=10.0):
    # Read in all the correlation files
    allfiles = [f for f in os.listdir(directory) if f.startswith(prefix)]

    #Figure out the primary mass
    MS = SpectralTypeRelations.MainSequence()
    header = pyfits.getheader(prefix + ".fits")
    starname = header["OBJECT"].split()[0].replace("_", " ")
    stardata = StarData.GetData(starname)
    primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])
    primary_temp = MS.Interpolate(MS.Temperature, stardata.spectype[:2])

    detections = defaultdict(list)
    outfile = open(outfilename, "w")
    outfile.write("Sensitivity Analysis:\n*****************************\n\n")
    outfile.write(
        "Filename\t\t\tPrimary Temperature\tSecondary Temperature\tMass (Msun)\tMass Ratio\tVelocity\tPeak Correct?\tSignificance\n")
    for fname in allfiles:
        #Read in temperature and expected velocity from filename
        T = float(fname.split(prefix)[-1].split("t")[-1].split("_")[0])
        v = float(fname.split("v")[-1])

        #Figure out the secondary mass from the temperature
        spt = MS.GetSpectralType(MS.Temperature, T)
        secondary_mass = MS.Interpolate(MS.Mass, spt)
        q = secondary_mass / primary_mass

        #Find the maximum in the cross-correlation function
        vel, corr = np.loadtxt(directory + fname, unpack=True)
        idx = np.argmax(corr)
        vmax = vel[idx]
        fit = FittingUtilities.Continuum(vel, corr, fitorder=2, lowreject=3, highreject=3)
        corr -= fit
        mean = corr.mean()
        std = corr.std()
        significance = (corr[idx] - mean) / std
        if np.abs(vmax - v) <= tolerance:
            #Signal found!
            outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" % (
            prefix, primary_temp, T, secondary_mass, q, v, significance))
        else:
            outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\t%.2f\n" % (
            prefix, primary_temp, T, secondary_mass, q, v, significance))
    outfile.close()


def MakePlot(infilename):
    # Set up thing to cycle through matplotlib linestyles
    from itertools import cycle

    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)

    #Defaults
    combine = False
    xaxis = "SecondarySpectralType"
    yaxis = "DetectionRate"


    #Command-line overrides
    for arg in sys.argv:
        if "combine" in arg:
            combine = True
        elif "xaxis" in arg:
            xaxis = arg.split("=")[-1]
        elif "yaxis" in arg:
            yaxis = arg.split("=")[-1]
        elif "infile" in arg:
            infilename = arg.split("=")[-1]

    if combine and "," in infilename:
        infiles = infilename.split(",")
    else:
        infiles = [infilename, ]

    #Set up dictionaries/lists
    p_spt = defaultdict(list)  #Primary spectral type
    s_spt = defaultdict(list)  #Secondary spectral type
    s_temp = defaultdict(list)  #Secondary Temperature
    p_mass = defaultdict(list)  #Primary mass
    s_mass = defaultdict(list)  #Secondary mass
    q = defaultdict(list)  #Mass ratio
    det_rate = defaultdict(list)  #Detection rate
    sig = defaultdict(list)  #Average detection significance
    magdiff = defaultdict(list)  #Magnitude difference
    namedict = {"SecondarySpectralType": s_spt,
                "SecondaryTemperature": s_temp,
                "SecondaryMass": s_mass,
                "MassRatio": q,
                "DetectionRate": det_rate,
                "AverageSignificance": sig,
                "MagnitudeDifference": magdiff}
    labeldict = {"SecondarySpectralType": "Secondary Spectral Type",
                 "SecondaryTemperature": "Secondary Temperature (K)",
                 "SecondaryMass": "SecondaryMass (Solar Masses)",
                 "MassRatio": "Mass Ratio",
                 "DetectionRate": "Detection Rate",
                 "AverageSignificance": "Average Significance",
                 "MagnitudeDifference": "Magnitude Difference"}

    if xaxis not in namedict.keys() or yaxis not in namedict:
        print "Error! axis keywords must be one of the following:"
        for key in namedict.keys():
            print key
        print "You chose %s for the x axis and %s for the y axis" % (xaxis, yaxis)
        sys.exit()

    MS = SpectralTypeRelations.MainSequence()
    vband = np.arange(500, 600, 1) * units.nm.to(units.cm)


    #Read in file/files  WARNING! ASSUMES A CERTAIN FORMAT. MUST CHANGE THIS IF THE FORMAT CHANGES!
    for infilename in infiles:
        infile = open(infilename)
        lines = infile.readlines()
        infile.close()
        print "Reading file %s" % infilename
        current_temp = float(lines[4].split()[2])
        starname = lines[4].split()[0].split("/")[-1].split("_smoothed")[0]
        detections = 0.0
        numsamples = 0.0
        significance = []
        for iternum, line in enumerate(lines[4:]):
            segments = line.split()
            #print segments[2], current_temp
            #print starname, segments[0]
            #print iternum, numsamples, '\n'
            if float(segments[2]) != current_temp or starname not in segments[0] or iternum == 0:
                if iternum != 0:
                    #We are on to the next temperature. Save info!
                    s_spt[starname].append(s_spectype)
                    s_temp[starname].append(current_temp)
                    p_spt[starname].append(p_spectype)
                    p_mass[starname].append(sec_mass / massratio)
                    s_mass[starname].append(sec_mass)
                    q[starname].append(massratio)
                    det_rate[starname].append(detections / numsamples)
                    sig[starname].append(np.mean(significance))
                    magdiff[starname].append(2.5 * np.log10(fluxratio))

                    #Reset things
                    current_temp = float(segments[2])
                    starname = segments[0].split("/")[-1].split("_smoothed")[0]
                    numsamples = 0.0
                    detections = 0.0
                    significance = []

                #Figure out the time-consuming SpectralType calls
                T1 = float(segments[1])
                p_spectype = MS.GetSpectralType(MS.Temperature, T1)
                R1 = MS.Interpolate(MS.Radius, p_spectype)
                T2 = float(segments[2])
                s_spectype = MS.GetSpectralType(MS.Temperature, T2)
                R2 = MS.Interpolate(MS.Radius, s_spectype)
                fluxratio = (PlotBlackbodies.Planck(vband, T1) / PlotBlackbodies.Planck(vband, T2)).mean() * (
                                                                                                             R1 / R2) ** 2

            else:
                #starname = segments[0].split("/")[-1].split("_smoothed")[0]
                sec_mass = float(segments[3])
                massratio = float(segments[4])
                if "y" in segments[6]:
                    detections += 1.
                    significance.append(float(segments[7]))
                numsamples += 1.




        #plot
        print "Plotting now"
        spt_sorter = {"O": 1, "B": 2, "A": 3, "F": 4, "G": 5, "K": 6, "M": 7}
        fcn = lambda s: (spt_sorter[itemgetter(0)(s)], itemgetter(1)(s))
        #print sorted(s_spt.keys(), key=fcn)
        #for starname in sorted(s_spt.keys(), key=fcn):
        print sorted(s_spt.keys())
        for starname in sorted(s_spt.keys()):
            p_spectype = p_spt[starname]
            x = namedict[xaxis][starname]
            y = namedict[yaxis][starname]
            if "SpectralType" in xaxis:
                plt.plot(range(len(x)), y[::-1], linestyle=next(linecycler), linewidth=2,
                         label="%s (%s)" % (starname, p_spectype[0]))
                plt.xticks(range(len(x)), x[::-1], size="small")
            elif "SpectralType" in yaxis:
                plt.plot(x[::-1], range(len(y)), linestyle=next(linecycler), linewidth=2,
                         label="%s (%s)" % (starname, p_spectype[0]))
                plt.yticks(range(len(y)), y[::-1], size="small")
            else:
                plt.plot(x, y, linestyle=next(linecycler), linewidth=2, label="%s (%s)" % (starname, p_spectype[0]))
                if "Magnitude" in xaxis:
                    ax = plt.gca()
                    ax.set_xlim(ax.get_xlim()[::-1])
                elif "Magnitude" in yaxis:
                    ax = plt.gca()
                    ax.set_ylim(ax.get_ylim()[::-1])
            plt.legend(loc='best')
            plt.xlabel(labeldict[xaxis])
            plt.ylabel(labeldict[yaxis])
            if "DetectionRate" in yaxis:
                ax = plt.gca()
                ax.set_ylim([-0.05, 1.05])
            plt.title("Sensitivity Analysis")

    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    if any(["new" in f for f in sys.argv[1:]]):
        directory = "Sensitivity/"
        allfiles = [f for f in os.listdir(directory) if (f.startswith("HIP") or f.startswith("HR"))]
        prefixes = []
        for fname in allfiles:
            prefix = fname.split("_v")[0][:-6]
            if prefix not in prefixes:
                print "New prefix: %s" % prefix
                prefixes.append(prefix)
        for i, prefix in enumerate(prefixes):
            MakeSummaryFile(directory, prefix, outfilename="%slogfile%i.txt" % (directory, i + 1))
            MakePlot("%slogfile%i.txt" % (directory, i + 1))
    else:
        MakePlot("Sensitivity/logfile.dat")
