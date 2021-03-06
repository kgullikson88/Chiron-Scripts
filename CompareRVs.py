import matplotlib.pyplot as plt
import numpy as np


badnums = [17457, 26126, 46283, 52678, 85537, 77635, 68520, 66249, 93805,
           101589, 105140, 109139, 23362, 23916, 28910, 52736, 60009,
           79387, 80815, 92855, 29735, 31362, 37297, 52736, 71865, 72104]
# badnums = []
bad_names = ["HIP %i" % num for num in badnums]

if __name__ == "__main__":
    known_file = open("Known_RVs.csv")
    known = known_file.readlines()
    known_file.close()

    measured_file = open("Measured_RVs.backup.csv")
    #measured_file = open("Known_RVs.backup.csv")
    measured = measured_file.readlines()
    measured_file.close()

    known_vels = []
    known_velerr = []
    measured_vels = []
    measured_velerr = []
    measured_vsini = []
    for i, line in enumerate(measured):
        segments = line.split("|")
        starname = segments[0].strip()
        if starname in bad_names:
            print "Skipping ", starname
            continue
        mrv = float(segments[1])  #measured rv
        mrv_err = float(segments[2])

        vsini = float(segments[4])
        #vsini = 1

        #Get the known rv
        for line2 in known:
            segments = line2.split("|")
            if segments[0].strip() == starname:
                krv = float(segments[1])
                krv_err = float(segments[2])

        #segments = known[i].split("|")
        #krv = float(segments[1])  
        #krv_err = float(segments[2])

        #if abs(mrv - krv) > 1:
        #    print "Large difference for %s: %g km/s" %(starname, mrv-krv)

        if abs(mrv - krv) > 20:
            continue

        known_vels.append(krv)
        known_velerr.append(krv_err)
        measured_vels.append(mrv)
        measured_velerr.append(mrv_err)
        measured_vsini.append(vsini)
        #plt.errorbar(krv, mrv-krv, xerr=krv_err, yerr=mrv_err)
    krv = np.array(known_vels)
    krv_err = np.array(known_velerr)
    mrv = np.array(measured_vels)
    mrv_err = np.array(measured_velerr)
    vsini = np.array(measured_vsini)
    lowv = np.where(vsini < np.median(vsini))[0]
    highv = np.where(vsini >= np.median(vsini))[0]
    #plt.errorbar(krv, mrv-krv, xerr=krv_err, yerr=mrv_err, fmt='ko')
    plt.errorbar(krv[lowv], (mrv - krv)[lowv], xerr=krv_err[lowv], yerr=mrv_err[lowv], fmt='bo')
    plt.errorbar(krv[highv], (mrv - krv)[highv], xerr=krv_err[highv], yerr=mrv_err[highv], fmt='ro')
    #plt.hist((mrv-krv)[lowv], cumulative=True, bins=50)
    #plt.hist((mrv-krv)[highv], cumulative=True, bins=50)
    print "Deviation = %g km/s" % np.std(mrv - krv)
    print "Offset = %g km/s" % np.mean(mrv - krv)
    plt.xlabel("Known RV (km/s)")
    plt.ylabel("Measured - Known RV (km/s)")
    #plt.ylabel("Measured RV (km/s")
    lims = plt.xlim()
    #plt.plot(lims, lims, 'k-')
    plt.show()
