__author__ = 'Kevin Gullikson'

"""
This script will add all of the available smoothed data,
  and output an extremely high S/N ratio spectrum. Since
  it is adding lots of different stars, it will essentially
  just be a good estimate of the systematic errors in my data
"""

import os
import FittingUtilities

import numpy as np
import matplotlib.pyplot as plt

import HelperFunctions


def Add(fileList):
    for i, fname in enumerate(fileList):
        orders = HelperFunctions.ReadExtensionFits(fname)
        if i == 0:
            master = []
            for order in orders:
                master.append(order.copy())
                master[-1].x = np.linspace(order.x[0], order.x[-1], order.size())
                master[-1].y = np.zeros(master[-1].size())

        for j, order in enumerate(orders):
            master[j].y += FittingUtilities.RebinData(order, master[j].x).y

    for order in master:
        order.y /= float(len(fileList))
        plt.plot(order.x, order.y, 'k-', alpha=0.4)
    plt.show()


if __name__ == "__main__":
    dirList = [d for d in os.listdir("./") if d.startswith("201") and len(d) == 8]
    all_files = []
    for dir in dirList:
        print dir
        fileList = ["{:s}/{:s}".format(dir, f) for f in os.listdir(dir) if
                    f.startswith("H") and f.endswith("smoothed.fits")]
        [all_files.append(f) for f in fileList]

    Add(all_files)

