import pyfits
import sys
import os
from collections import defaultdict
import subprocess
import HelperFunctions



if __name__ == "__main__":
  outdir = "files_%s/" %(os.getcwd().split("/")[-1])
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
      print "Unrecognized argument: %s" %arg

  HelperFunctions.ensure_dir(outdir)  
  
  #Make a list of the science data, checking some information in the header
  allscience = [f for f in os.listdir(datadir) if f.startswith("chi")]
  #allfile = open("inlist", "w")
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
      print "Calibration image found in data directory: %s\n  Moving to calibration directory" %fname
      cmd = "mv %s%s %s" %(datadir, fname, caldir)
      subprocess.check_call(cmd, shell=True)
      continue
    
    #Copy to output directory
    cmd = "cp %s%s %s" %(datadir, fname, outdir)
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
    header = pyfits.getheader(caldir + "/" + fname)
    otype = header["OBJECT"].lower()
    if header["ccdsum"].strip() == ccdsum and abs(header['deckpos'] - decker) < 0.001 and otype != "thar":
      #Copy to output directory
      cmd = "cp %s%s %s" %(caldir, fname, outdir)
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
      
      FileDict[otype].append(fname)

  for otype in FileDict.keys():
    print otype
    outfile = open("%s.list" %otype, "w")
    files = FileDict[otype]
    for fname in files:
      print "\t", fname
      #outfile.write("%s/%s\n" %(caldir, fname))
      #allfile.write("%s/%s\n" %(caldir, fname))
    #outfile.close()
  #allfile.close()

  

    
