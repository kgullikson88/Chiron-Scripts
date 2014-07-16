import numpy as np
import HelperFunctions
import StarData
import matplotlib.pyplot as plt 
import sys
import Correlate
from astropy import units
import Broaden
import FittingUtilities
import DataStructures
from scipy.optimize import curve_fit
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astrolib import helcorr
import SpectralTypeRelations
from collections import defaultdict
import os

modelfile = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/lte98-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"

H_orders = [6,7,37]  #Don't use H lines in vsini determination
tell_orders = [27,33,41,44,45,49, ]
good_orders = [6,7,13,32,35, 37,50]
good_orders = [6,7,37,50]
good_orders = [6,7,37]
good_orders = [6,7]
hightemp_orders = [2,3,8,10,13,14,16,17,19,26,28,31,34,38,43,50,51]
lowtemp_orders = [0,1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,30,31,32,33,34,35,36,43,44,45,46,47,48,49,50,51,52]
c = 3e5

Grid_Temperatures = [6750,
7000,
7250,
7500,
7750,
8000,
8250,
8500,
8750,
9000,
9250,
9500,
9750,
10000,
10250,
10500,
10750,
11000,
11250,
11500,
11750,
12000,
12250,
12500,
12750,
13000,
14000,
15000,
16000,
17000,
18000,
19000,
20000,
21000]
Grid_Temperatures = np.array(Grid_Temperatures)


def convert(coord, delim=":"):
    segments = coord.split(delim)
    s = np.sign(float(segments[0]))
    return s * (abs(float(segments[0])) + float(segments[1])/60.0 + float(segments[2])/3600.0)


def getJD(header, rootdir="./"):
    """
    Use the header of a co-added file to determine the julian date
    in the middle of the total observation
    """
    # First, get the filenames that were used
    fileList = [k for k in header.keys() if "FILE" in k]
    fileheader = fits.getheader(rootdir + header[fileList[0]])
    firstjd = fileheader['jd']
    fileheader = fits.getheader(rootdir + header[fileList[-1]])
    lastjd = fileheader['jd'] + fileheader['exptime']/(24.0*3600.0)
    return (firstjd + lastjd)/2.0



def GetModel(Temperature, grid):
    filename = "%st%.5i_g+4.0_p00p00_hr.fits" %(grid, Temperature)
    print "Model file name: ", filename
    hdu = fits.open(filename)[0]
    data = hdu.data
    header = hdu.header
    wave_A = np.arange(data.size)*header['CDELT1'] + header['CRVAL1']
    x = wave_A*units.angstrom.to(units.nm)
    y = data
    left = np.searchsorted(x, 450)
    right = np.searchsorted(x, 900)
    model = DataStructures.xypoint(x=x[left:right], y=y[left:right])
    return model



def profile(x, a, b, amp, mu, sig, vsini):
    a2 = 1.0   #Limb darkening parameter
    b2 = -0.1  #Limb darkening parameter
    c = 3e5
    gauss = a + b*x + amp*np.exp(-(x-mu)**2/(2*sig**2))
    broadening_kernel = a2*np.sqrt(1-(x/vsini)**2) + b2*np.pi/4.0 * (1-(x/vsini)**2)
    keep = -np.isnan(broadening_kernel)
    broadening_kernel = broadening_kernel[keep]
    return np.convolve(gauss, broadening_kernel/broadening_kernel.sum(), mode='same')



if __name__ == "__main__":
    MS = SpectralTypeRelations.MainSequence()
    
    multiplicity_file = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Multiplicity.csv" %(os.environ['HOME'])
    infile = open(multiplicity_file)
    multiplicity = infile.readlines()
    infile.close()

    if "darwin" in sys.platform:
        modeldir = "/Volumes/DATADRIVE/Stellar_Models/Coelho14/s_coelho14_highres/"
    elif "linux" in sys.platform:
        modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Stellar/Coelho14/s_coelho14_highres/"
    else:
        modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

    logfile = open("Measured_RVs.csv", "w")
    
    first = True
    models = defaultdict(str)
    for fnum, fname in enumerate(sys.argv[1:]):
        rootdir = "./"
        dirs = fname.split("/")
        for d in dirs[:-1]:
            rootdir = rootdir + d + "/"
        header = fits.getheader(fname)
        starname = header['OBJECT']
        print starname

        #Make sure this is not a multiple star
        skip = False
        for star in multiplicity:
            segments = star.split("|")
            if starname in segments[0]:
                nummult = int(segments[4])
                if nummult > 0:
                    print "%s is in a multiple system! Skipping!" %starname
                    skip = True
        if skip:
            continue

        # Get the appropriate model for this temperature
        data = StarData.GetData(starname)
        spt = data.spectype[:2]
        Teff = MS.Interpolate(MS.Temperature, spt)
        idx = np.argmin(abs(Grid_Temperatures - Teff))
        Tgrid = Grid_Temperatures[idx]
        if models[Tgrid] == "":
            models[Tgrid] = GetModel(Tgrid, modeldir)
        model = models[Tgrid]
        if Teff > 9000:
            good_orders = hightemp_orders
        else:
            good_orders = lowtemp_orders
        

        # Get the heliocentric correction
        ra = convert(header['RA'])
        dec = convert(header['DEC'])
        jd = getJD(header, rootdir=rootdir)
        longitude = 70.8065
        latitude = -30.1697
        altitude = 2200.0
        vbary = helcorr(longitude, latitude, altitude, ra, dec, jd)[0]
        print "Heliocentric velocity correction: ", vbary


        A = 1.85
        sig0 = 0.6
        all_orders = HelperFunctions.ReadExtensionFits(fname)
        orders = []
        for i, o in enumerate(all_orders):
            if i in good_orders and i not in tell_orders:
                orders.append(o)
            """
            if i in good_orders:
                #o.cont = FittingUtilities.Continuum(o.x, o.y, fitorder=1, lowreject=2, highreject=20)
                o.cont = np.mean(o.y) * np.ones(o.size())
                left = np.searchsorted(model.x, o.x[0]-10)
                right = np.searchsorted(model.x, o.x[-1]+10)
                m = model[left:right].copy()
                #m = FittingUtilities.RebinData(model, o.x)
                m.cont = FittingUtilities.Continuum(m.x, m.y, fitorder=2, lowreject=1.5, highreject=10)
                m = FittingUtilities.RebinData(m, o.x)
                o2 = o.copy()
                o2.y /= (m.y/m.cont)
                o.cont = FittingUtilities.Continuum(o2.x, o2.y, fitorder=1, lowreject=2, highreject=2)
                o = HelperFunctions.Denoise(o)
                orders.append(o)
                #plt.figure(1)
                #plt.plot(m.x, m.y)
                #plt.plot(m.x, m.cont)
                #plt.figure(2)
                #plt.plot(o2.x, o2.y)
                #plt.plot(o2.x, o2.cont)
                #plt.figure(3)
                #plt.plot(o.x, o.y/o2.cont)
                #plt.plot(o2.x, o2.cont)
            """
        #plt.show()

        if first:
            retdict = Correlate.GetCCF(orders, model, vsini=0.0, resolution=900000, addmode="simple", debug=True, oversample=1)
            model_orders = retdict['model']
            #first = False
        else:
            retdict = Correlate.GetCCF(orders, model_orders, vsini=0.0, resolution=3000, addmode='simple', process_model=False, oversample=1)
        ccf = retdict['CCF']
        gauss = lambda x, a, b, amp, mu, sig: a + b*x + amp*np.exp(-(x-mu)**2/(2*sig**2))
        left = np.searchsorted(ccf.x, -200)
        right = np.searchsorted(ccf.x, 200)
        idx = np.argmax(ccf.y[left:right])+left
        y1 = np.median(ccf.y[:50])
        y2 = np.median(ccf.y[-50:])
        m = (y2 - y1) / (ccf.x[-25] - ccf.x[25])
        a = y1 - m*ccf.x[25]
        b = m
        amp = ccf.y[idx] - a
        mu = ccf.x[idx]
        sig = 10.0
        vsini = 100.0
        pars = [a,b,amp, mu, sig, vsini]
        try:
            popt, pcov = curve_fit(profile, ccf.x, ccf.y, p0=(pars))
        except RuntimeError:
            continue

        #Get centroid of peak
        a, b = popt[0], popt[1]
        vel = ccf.x[ccf.y > 0.5*ccf.y[idx]]
        corr = ccf.y[ccf.y > 0.5*ccf.y[idx]]# - (a + b*vel)
        centroid = -np.sum(vel*corr) / np.sum(corr)
        """
        gauss2 = lambda x,amp,mu,sig: amp*np.exp(-(x-mu)**2/(2*sig**2))
        try:
            popt2, pcov2 = curve_fit(gauss2, vel, corr, p0=(popt[2:]))
        except RuntimeError:
            continue
        centroid = -popt2[1]
        """
        centroid = -ccf.x[idx]
        



        rv = -popt[3]
        rv_err = np.sqrt(pcov[3][3])
        #vsini = A*np.sqrt(popt[4]**2 - sig0**2)
        #vsini_err = A*np.sqrt(pcov[4][4]) /  np.sqrt(popt[4]**2 - sig0**2)
        vsini = popt[5]
        vsini_err = np.sqrt(pcov[5][5])
        print "rv = %g km/s" %rv
        print "centroid = %g km/s" %centroid
        print "vsini = %g km/s" %vsini

        ccf.output(fname.replace(".fits", "_CCF.txt"))
        plt.plot(ccf.x, ccf.y)
        plt.plot(ccf.x, profile(ccf.x, *popt))
        plt.plot((-centroid, -centroid), (ccf.y[idx]+0.05, ccf.y[idx]+0.15), 'k-')
        plt.savefig("Figures/%s" %(fname.split("/")[-1].replace("fits", "pdf")))
        plt.clf()
        #plt.show()

        logfile.write("%s  |  %g  |  %g  |  %g  |  %g  |  %g\n" %(starname, rv+vbary, rv_err, centroid+vbary, vsini, vsini_err))
        
        """
        model2 = Broaden.RotBroad(model, vsini*units.km.to(units.cm))
        model2 = FittingUtilities.ReduceResolution2(model2, 80000.0)
        model2.x *= (1-rv/c)
        fcn = spline(model2.x, model2.y)
        for order in all_orders:
            segment = order.copy()
            segment.y = fcn(segment.x)
            #segment = FittingUtilities.RebinData(model2, order.x)
            segment.cont = FittingUtilities.Continuum(segment.x, segment.y, fitorder=3, lowreject=1, highreject=10)
            plt.plot(order.x, order.y/order.cont, 'k-')
            plt.plot(segment.x, segment.y/segment.cont, 'r-')
        plt.show()
        """
    logfile.close()
