import sys
import os
import time

from astropy.io import fits as pyfits
from astropy import units, constants

import TelluricFitter
import FittingUtilities
import HelperFunctions
import GetAtmosphere


sleep = True
homedir = os.environ["HOME"]

badregions = [[588.98, 589.037],  # Na D line 1
              [589.567, 589.632],  # Na D line 2
              [627.4, 629.0],  #O2 band
              [686.4, 690.7]]  # O2 band

namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
            "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
            "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
            "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
            "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
            "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
            "o3": ["O3FIT", "O3VAL", "O3 abundance"],
            "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
            "co": ["COFIT", "COVAL", "CO abundance"],
            "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
            "o2": ["O2FIT", "O2VAL", "O2 abundance"],
            "no": ["NOFIT", "NOVAL", "NO abundance"],
            "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
            "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
            "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
            "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}


def FindOrderNums(orders, wavelengths):
    """
      Given a list of xypoint orders and
      another list of wavelengths, this
      finds the order numbers with the
      requested wavelengths
    """
    nums = []
    for wave in wavelengths:
        for i, order in enumerate(orders):
            if order.x[0] < wave and order.x[-1] > wave:
                nums.append(i)
                break
    return nums


if __name__ == "__main__":
    # Initialize fitter
    fitter = TelluricFitter.TelluricFitter()
    fitter.SetObservatory("CTIO")

    fileList = []
    start = 0
    end = 999
    makenew = True
    edit_atmosphere = False
    humidity_low = 1.0
    humidity_high = 99.0
    for arg in sys.argv[1:]:
        if "-atmos" in arg:
            edit_atmosphere = True
        elif "-hlow" in arg:
            humidity_low = float(arg.split("=")[1])
        elif "-hhigh" in arg:
            humidity_high = float(arg.split("=")[1])
        else:
            fileList.append(arg)


    # START LOOPING OVER INPUT FILES
    for fname in fileList:
        name = fname.split(".fits")[0]
        outfilename = "Corrected_%s.fits" % name

        #Read file
        orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")

        header = pyfits.getheader(fname)
        angle = float(header["ZD"])
        resolution = 80000.0
        humidity = max(header["OUTHUM"], 5)
        pressure = header["OUTPRESS"]
        temperature = header["OUTTEMP"] + 273.15

        if edit_atmosphere:
            filenames = [f for f in os.listdir("./") if "GDAS" in f]
            height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'].split("T")[0],
                                                               header['ut'])

            fitter.EditAtmosphereProfile("Temperature", height, Temp)
            fitter.EditAtmosphereProfile("Pressure", height, Pres)
            fitter.EditAtmosphereProfile("H2O", height, h2o)


        #Adjust fitter values
        fitter.FitVariable({"h2o": humidity})
        #                    "temperature": temperature})
        #                    "o2": 2.12e5})
        fitter.AdjustValue({"angle": angle,
                            "pressure": pressure,
                            "temperature": temperature,
                            "resolution": resolution,
                            "o2": 2.12e5})
        fitter.SetBounds({"h2o": [humidity_low, humidity_high],
                          "temperature": [temperature - 10, temperature + 10],
                          "o2": [5e4, 1e6],
                          "resolution": [70000, 90000]})

        #Ignore the interstellar sodium D lines and parts of the O2 bands
        fitter.IgnoreRegions(badregions)
        models = []

        #Read in the header from the corrected_file
        corrected_file = "Corrected_%s" % fname
        header = pyfits.getheader(corrected_file, ext=2)
        humidity = header['H2OVAL']
        o2 = header['O2VAL']
        resolution = header['RESVAL']
        temperature = header['TEMPVAL']
        fitter.AdjustValue({"o2": o2,
                            "h2o": humidity,
                            "resolution": resolution,
                            "temperature": temperature})
        vel = 4.4

        # Finally, apply these parameters to all orders in the data
        for i, order in enumerate(orders):
            print "\n\nGenerating model for order %i of %i\n" % (i, len(orders))
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
            fitter.ImportData(order)
            fitter.resolution_fit_mode = "gauss"
            #fitter.resolution_fit_mode = "svd"
            #wave0 = order.x.mean()
            #fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
            print "fitter.shift = ", fitter.shift
            primary, model = fitter.GenerateModel(fitpars,
                                                  separate_source=True,
                                                  return_resolution=False)

            data = fitter.data
            if min(model.y) > 0.98:
                #The wavelength calibration might be off
                wave0 = order.x.mean()
                fitter.shift = vel / (constants.c.cgs.value * units.cm.to(units.km)) * wave0
                model = fitter.GenerateModel(fitpars, separate_source=False, nofit=True)
                model.x /= (1.0 + vel / (constants.c.cgs.value * units.cm.to(units.km)))
                model = FittingUtilities.RebinData(model, order.x)
                data = order.copy()
                data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=3, lowreject=2, highreject=5)

            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err,
                       "model": model.y,
                       "primary": primary.y}

            header_info = []
            numpars = len(fitter.const_pars)
            for j in range(numpars):
                try:
                    parname = fitter.parnames[j]
                    parval = fitter.const_pars[j]
                    fitting = fitter.fitting[j]
                    header_info.append([namedict[parname][0], fitting, namedict[parname][2]])
                    header_info.append([namedict[parname][1], parval, namedict[parname][2]])
                except KeyError:
                    print "Not saving the following info: %s" % (fitter.parnames[j])

            if (i == 0 and makenew) or not exists:
                HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info, ],
                                                         mode="new")
                exists = True
            else:
                HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename,
                                                         headers_info=[header_info, ], mode="append")

            if sleep:
                time.sleep(30)
