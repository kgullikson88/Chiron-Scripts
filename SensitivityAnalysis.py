import os
import sys

from astropy import units
import Sensitivity

import Search_Fast


if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"


if __name__ == "__main__":
    # Parse command-line arguments
    vsini_secondary = 20 * units.km.to(units.cm)
    resolution = 80000
    smooth_factor = 0.8
    vel_list = range(-400, 400, 50)
    companion_file = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Multiplicity.csv" % (os.environ["HOME"])
    vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" % (os.environ["HOME"])
    fileList = []
    tolerance = 5.0
    debug = False
    vsini_skip = 10
    vsini_idx = 1
    trimsize = 100
    for arg in sys.argv[1:]:
        if "-m" in arg:
            companion_file = arg.split("=")[1]
        elif "-tol" in arg:
            tolerance = float(arg.split("=")[1])
        elif "-d" in arg:
            debug = True
        elif "-vsinifile" in arg:
            vsini_file = arg.split("=")[-1]
        elif "-vsiniskip" in arg:
            vsini_skip = int(arg.split("=")[-1])
        elif "-vsiniidx" in arg:
            vsini_idx = int(arg.split("=")[-1])
        else:
            fileList.append(arg)

    Sensitivity.Analyze(fileList,
                        resolution=resolution,
                        debug=True,
                        badregions=Search_Fast.badregions,
                        trimsize=trimsize,
                        modeldir=modeldir,
                        companion_file=companion_file,
                        vsini_file=vsini_file,
                        vsini_idx=vsini_idx,
                        vsini_skip=vsini_skip,
                        vsini_secondary=vsini_secondary)