import sys
import HelperFunctions
import matplotlib.pyplot as plt
import numpy as np



if __name__ == "__main__":
    trimsize = 100
    fileList = []
    for arg in sys.argv[1:]:
        fileList.append(arg)

    for fname in fileList:
        orders = HelperFunctions.ReadExtensionFits(fname)
        maxsnr = -np.inf
        for order in orders:
            snr = order.y.mean() / order.y.std()
            if snr > maxsnr:
                maxsnr = snr
        print fname, maxsnr