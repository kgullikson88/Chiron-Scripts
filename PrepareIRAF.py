import sys
import os
from collections import defaultdict
import subprocess

from astropy.io import fits as pyfits

import HelperFunctions
import ConvertToExtensions_NativeReductions as convert


if __name__ == "__main__":
    date = os.getcwd().split("/")[-1]
    outdir = "files_%s/" % (date)
    tarfiles = [f for f in os.listdir("./") if ".tgz" in f]
    caldir = [f.split(".tgz")[0] for f in tarfiles if "cals" in f][0] + "/"
    datadir = [f.split(".tgz")[0] for f in tarfiles if "planid" in f][0] + "/"
    for arg in sys.argv[1:]:
        if "-datadir" in arg:
            datadir = arg.split("=")[-1]
            if not datadir.endswith("/"):
                datadir = datadir + "/"
        elif "-caldir" in arg:
            caldir = arg.split("=")[-1]
            if not caldir.endswith("/"):
                caldir = caldir + "/"
        elif "output" in arg:
            outdir = arg.split("=")[-1]
            if not outdir.endswith("/"):
                outdir = outdir + "/"
        else:
            print "Unrecognized argument: %s" % arg

    HelperFunctions.ensure_dir(outdir)

    # Make a list of the science data, checking some information in the header
    allscience = [f for f in os.listdir(datadir) if f.startswith("chi")]
    # allfile = open("inlist", "w")
    #outfile = open("science.list", "w")
    print "Science data:"
    fname = allscience[0]
    header = pyfits.getheader(datadir + fname)
    decker = header['deckpos']
    ccdsum = header['ccdsum'].strip()
    for fname in allscience:
        header = pyfits.getheader(datadir + fname)
        #print fname, header['IMAGETYP'], header['object'], header['deckpos']
        if header['IMAGETYP'].strip().lower == "calibration" and header['object'].lower() != 'thar':
            print "Calibration image found in data directory: %s\n  Moving to calibration directory" % fname
            cmd = "mv %s%s %s" % (datadir, fname, caldir)
            subprocess.check_call(cmd, shell=True)
            continue

        #Copy to output directory
        cmd = "cp %s%s %s" % (datadir, fname, outdir)
        subprocess.check_call(cmd, shell=True)

        if header['object'].lower() == 'thar':
            #Edit the fits header a bit for iraf
            hdulist = pyfits.open(outdir + fname, mode='update')
            hdulist[0].header['IMAGETYP'] = 'comp'
            hdulist.flush()
            hdulist.close()
        else:
            print "\t", fname
        if abs(header['deckpos'] - decker) > 0.001 or header['ccdsum'].strip() != ccdsum:
            print "Error! deckpos and ccdsum are not the same in all science observations!"
            sys.exit()

    #Make lists of the appropriate calibration data
    allcalib = [f for f in os.listdir(caldir) if f.startswith("chi")]
    FileDict = defaultdict(list)

    for fname in allcalib:
        if "slicer" in fname:
            cmd = "cp %s%s ." % (caldir, fname)
            subprocess.check_call(cmd, shell=True)
        header = pyfits.getheader(caldir + "/" + fname)
        if not "OBJECT" in header:
            print "Header of %s had no 'object' keyword!" % fname
            continue
        otype = header["OBJECT"].lower()
        if header["ccdsum"].strip() == ccdsum and abs(header['deckpos'] - decker) < 0.001:  # and otype != "thar":
            #Copy to output directory
            cmd = "cp %s%s %s" % (caldir, fname, outdir)
            subprocess.check_call(cmd, shell=True)

            #Edit fits files as needed
            if otype.lower() == "bias":
                #Change imagetyp keyword from 'bias' to 'zero'
                hdulist = pyfits.open(outdir + fname, mode='update')
                hdulist[0].header['IMAGETYP'] = 'zero'
                hdulist.flush()
                hdulist.close()
            if otype.lower() == 'quartz':
                #Change imagetyp keyword to 'flat'
                hdulist = pyfits.open(outdir + fname, mode='update')
                hdulist[0].header['IMAGETYP'] = 'flat'
                hdulist.flush()
                hdulist.close()
            if otype.lower() == 'thar':
                hdulist = pyfits.open(outdir + fname, mode='update')
                hdulist[0].header['IMAGETYP'] = 'comp'
                hdulist.flush()
                hdulist.close()

            FileDict[otype].append(fname)

    for otype in FileDict.keys():
        print otype
        #outfile = open("%s.list" %otype, "w")
        files = FileDict[otype]
        for fname in files:
            print "\t", fname

    #copy all achi* files from the science directory to this directory
    cmd = "cp %sachi* ." % (datadir)
    subprocess.check_call(cmd, shell=True)

    achi_files = [f for f in os.listdir("./") if f.startswith("achi")]
    convert.Convert(achi_files, True, False)

    #Copy the files to kepler
    outdir = outdir.strip("/")
    subprocess.check_call(["scp", "-r", outdir, "kgulliks@kepler:~/"])

    #tar.gz the useful files
    archive_dir = "/Volumes/DATADRIVE/CHIRON_data/%s/" % date
    subprocess.check_call(["tar", "czvf", "%s%s.tar.gz" % (archive_dir, outdir), outdir])

    #Remove the unnecessary files
    cmd = "rm -r 15* files*"
    subprocess.check_call(cmd, shell=True)

  

    
