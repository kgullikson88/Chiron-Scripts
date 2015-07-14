import sys

import matplotlib.pyplot as plt
import FittingUtilities
import numpy as np

import HelperFunctions


if __name__ == "__main__":
    fileList = []
    tellurics = False
    normalize = False
    byorder = False  # Plots one order at a time
    pixelscale = False
    oneplot = False
    errors = False
    for arg in sys.argv[1:]:
        if "tellcorr" in arg:
            tellurics = True
        elif "-norm" in arg:
            normalize = True
        elif "-order" in arg:
            byorder = True
        elif "-pix" in arg:
            pixelscale = True
            # byorder = True
        elif "-one" in arg:
            oneplot = True
        elif "-err" in arg:
            errors = True
        else:
            fileList.append(arg)

    linestyles = ['k-', 'r-', 'b-', 'g-']

    for fnum, fname in enumerate(fileList):
        ls = linestyles[fnum % len(linestyles)]
        orders = HelperFunctions.ReadExtensionFits(fname)

        print fname, len(orders)
        if not oneplot:
            plt.figure(fnum)
            plt.title(fname)

        for i, order in enumerate(orders):

            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=1, lowreject=3, highreject=3)
            # order.cont = np.median(order.y)*np.ones(order.size())
            if pixelscale:
                order.x = np.arange(order.size())

            if normalize:
                plt.plot(order.x, order.y / order.cont, ls, rasterized=True, alpha=0.4)
                plt.text(order.x.mean(), 1.1, str(i + 1))
            else:
                if errors:
                    plt.errorbar(order.x, order.y, yerr=order.err)
                else:
                    plt.plot(order.x, order.y, ls)
                    #plt.plot(order.x, order.cont)
            if byorder:
                plt.title("Order %i" % i)
                plt.show()
        plt.plot([], [], ls, label=fname)
    if not byorder:
        plt.xlabel("Wavelength (nm)", fontsize=15)
        if normalize:
            plt.ylabel("Normalized Flux", fontsize=15)
            plt.ylim((-0.1, 1.2))
        else:
            plt.ylabel("Flux", fontsize=15)
        leg = plt.legend(loc='best', fancybox=True)
        plt.show()
