import sys
import FittingUtilities

import numpy as np
from astropy.io import fits as pyfits
import TelluricFitter

import HelperFunctions


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
        # Make output filename
        name = fname.split(".fits")[0]
        outfilename = u"Corrected_{0:s}.fits".format(name)

        # Read file
        orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")

        header = pyfits.getheader(fname)
        angle = float(header["ZD"])
        resolution = 80000.0
        humidity = max(header["OUTHUM"], 5)
        pressure = header["OUTPRESS"]
        temperature = header["OUTTEMP"] + 273.15

        #Adjust fitter values
        fitter.AdjustValue({"angle": angle,
                            "pressure": pressure,
                            "resolution": resolution,
                            "temperature": temperature,
                            "o2": 2.12e5})
        fitter.FitVariable({"h2o": humidity})
        #                    "temperature": temperature})
        fitter.SetBounds({"resolution": [70000, 90000]})

        #Ignore the interstellar sodium D lines and parts of the O2 bands
        fitter.IgnoreRegions(badregions)

        # Determine the H2O scale factor
        h2o_scale = []
        fitter.DisplayVariables()
        for i in HelperFunctions.FindOrderNums(orders, [595, 700, 717, 730]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0,
                                "o2": 0.0,
                                "h2o": humidity,
                                "resolution": resolution})
            fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
            fitter.ImportData(order)
            fitter.resolution_fit_mode = "gauss"
            fitter.fit_source = False
            fitter.fit_primary = False
            model = fitter.GenerateModel(fitpars, separate_primary=False, return_resolution=False)

            # Find the best scale factor
            model.cont = np.ones(model.size())
            lines = FittingUtilities.FindLines(model, tol=0.95).astype(int)
            if len(lines) > 5:
                scale = np.median(np.log(order.y[lines] / order.cont[lines]) / np.log(model.y[lines]))
            else:
                scale = 1.0
            print i, scale
            h2o_scale.append(scale)

        # Now, find the best O2 scale factor
        o2_scale = []
        for i in HelperFunctions.FindOrderNums(orders, [630, 690]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0,
                                "o2": 2.12e5,
                                "h2o": 0.0,
                                "resolution": resolution})
            fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
            fitter.ImportData(order)
            fitter.resolution_fit_mode = "gauss"
            fitter.fit_source = False
            fitter.fit_primary = False
            model = fitter.GenerateModel(fitpars, separate_primary=False, return_resolution=False)

            # Find the best scale factor
            model.cont = np.ones(model.size())
            lines = FittingUtilities.FindLines(model, tol=0.95).astype(int)
            if len(lines) > 5:
                scale = np.median(np.log(order.y[lines] / order.cont[lines]) / np.log(model.y[lines]))
            else:
                scale = 1.0
            print i, scale
            o2_scale.append(scale)

        # Use the median values
        o2_scale = np.median(o2_scale)
        h2o_scale = np.median(h2o_scale)

        # Now, apply the scale to everything
        o2 = 2.12e5 * o2_scale
        h2o_ppmv = TelluricFitter.MakeModel.humidity_to_ppmv(humidity, temperature, pressure)
        humidity = TelluricFitter.MakeModel.ppmv_to_humidity(h2o_ppmv * h2o_scale, temperature, pressure)

        exists = False
        for i, order in enumerate(orders):
            if (order.x[0] < 470 and order.x[-1] > 470) or max(order.y) < 0.01:
                model = order.copy()
                model.y = np.ones(order.size())
                data = order.copy()
                data.cont = np.ones(data.size())
            else:
                print "\n\nGenerating model for order %i of %i\n" % (i, len(orders))
                fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                    "waveend": order.x[-1] + 20.0,
                                    "o2": o2,
                                    "h2o": humidity,
                                    "resolution": resolution})
                fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
                order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
                fitter.ImportData(order)
                fitter.resolution_fit_mode = "gauss"
                #wave0 = order.x.mean()
                #fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
                fitter.DisplayVariables()
                model = fitter.GenerateModel(fitpars, separate_primary=False, return_resolution=False)

                data = fitter.data

            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err,
                       "model": model.y,
                       "primary": np.ones(data.size())}

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

