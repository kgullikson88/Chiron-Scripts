import sys
import os
import warnings
import FittingUtilities

from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

import DataStructures
import HelperFunctions





# plot = True


def ReadCorrectedFile(fname, yaxis="model"):
    orders = []
    headers = []
    hdulist = pyfits.open(fname)
    numorders = len(hdulist)
    for i in range(1, numorders):
        order = hdulist[i].data
        xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                      y=order.field(yaxis),
                                      cont=order.field("continuum"),
                                      err=order.field("error"))

        orders.append(xypt)
        headers.append(hdulist[i].header)
    return orders, headers


def Correct(original, corrected, offset=None, get_primary=False, plot=False):
    # Read in the data and model
    original_orders = HelperFunctions.ReadFits(original, extensions=True, x="wavelength", y="flux", errors="error",
                                               cont="continuum")
    corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
    test_orders, header = ReadCorrectedFile(corrected, yaxis="flux")

    if plot:
        for order, model in zip(test_orders, corrected_orders):
            plt.plot(order.x, order.y / order.cont)
            plt.plot(model.x, model.y)
        plt.title("Correction in corrected file only")
        plt.show()

    if get_primary:
        primary_orders = ReadCorrectedFile(corrected, yaxis="primary")[0]
    if offset == None:
        offset = len(original_orders) - len(corrected_orders)
    print "Offset = ", offset
    for i in range(offset, len(original_orders)):
        data = original_orders[i]
        data.cont = FittingUtilities.Continuum(data.x, data.y)
        try:
            model = corrected_orders[i - offset]
            header = corrected_headers[i - offset]
            if get_primary:
                primary = primary_orders[i - offset]
            try:
                if i == 3:
                    print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" % (i, header['h2oval'], header['o2val'])
            except KeyError:
                warnings.warn("Atmospheric parameters not found in the header!")
        except IndexError:
            model = DataStructures.xypoint(x=data.x, y=np.ones(data.x.size))
            print "Warning!!! Telluric Model not found for order %i" % i

        if plot:
            plt.figure(1)
            plt.plot(data.x, data.y / data.cont)
            plt.plot(model.x, model.y)

        if model.size() < data.size():
            left = np.searchsorted(data.x, model.x[0])
            right = np.searchsorted(data.x, model.x[-1])
            if right < data.size():
                right += 1
            data = data[left:right]
        elif model.size() > data.size():
            sys.exit("Error! Model size (%i) is larger than data size (%i)" % (model.size(), data.size()))

        # if np.sum((model.x-data.x)**2) > 1e-8:
        # model = FittingUtilities.RebinData(model, data.x)

        data.y[data.y / data.cont < 1e-5] = 1e-5 * data.cont[data.y / data.cont < 1e-5]
        badindices = np.where(np.logical_or(data.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = data.y[badindices] / data.cont[badindices]
        model.y[model.y < 1e-5] = 1e-5

        # plt.plot(data.x, data.y / model.y)
        data.y /= model.y
        data.err /= model.y
        if get_primary:
            data.y /= primary.y
        original_orders[i] = data.copy()
    if plot:
        plt.show()
    return original_orders


def main1():
    primary = False
    plot = False
    if len(sys.argv) > 2:
        original = sys.argv[1]
        corrected = sys.argv[2]
        # plot = True
        if len(sys.argv) > 3:
            if "prim" in sys.argv[3]:
                primary = True
            elif "plot" in sys.argv[3]:
                plot = True

        outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
        print "Outputting to %s" % outfilename

        corrected_orders = Correct(original, corrected, offset=None, get_primary=primary, plot=plot)

        column_list = []
        if plot:
            plt.figure(2)
        for i, data in enumerate(corrected_orders):
            if plot:
                plt.plot(data.x, data.y / data.cont)
                # plt.plot(data.x, data.cont)
            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err}
            column_list.append(columns)
        if plot:
            plt.title("Corrected data")
            plt.show()
        HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

    else:
        if len(sys.argv) > 1 and 'plot' in sys.argv[1]:
            plot = True
        allfiles = os.listdir("./")
        corrected_files = [f for f in allfiles if "Corrected_" in f and f.endswith("-0.fits")]
        # original_files = [f for f in allfiles if any(f in cf for cf in corrected_files)]

        # print corrected_files
        # print original_files

        for corrected in corrected_files:
            original = corrected.split("Corrected_")[-1]  # .split("-")[0] + ".fits"
            if original not in allfiles:
                warnings.warn("File %s not found. Skipping!" % original)
                continue
            # original = [f for f in allfiles if (f in corrected and f != corrected)]
            # if len(original) == 1:
            #  original = original[0]
            #else:
            #  sys.exit("Error! %i matches found to corrected file %s" %(len(original), corrected))

            print corrected, original
            outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
            print "Outputting to %s" % outfilename

            corrected_orders = Correct(original, corrected, offset=None, plot=plot)

            column_list = []
            if plot:
                plt.figure(2)
            for i, data in enumerate(corrected_orders):
                if plot:
                    plt.plot(data.x, data.y / data.cont)
                #Set up data structures for OutputFitsFile
                columns = {"wavelength": data.x,
                           "flux": data.y,
                           "continuum": data.cont,
                           "error": data.err}
                column_list.append(columns)
            HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

            if plot:
                plt.title(original)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Flux")
                plt.show()


if __name__ == "__main__":
    main1()
