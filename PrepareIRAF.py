import pyfits
import sys
import os
from collections import defaultdict

caldir = "calibrations"
datadir = "data"
caldir = "130201_cals"
datadir = "130201_planid_115"

if __name__ == "__main__":
  #Make a list of the science data, checking some information in the header
  allscience = [f for f in os.listdir(datadir) if f.startswith("chi")]
  allfile = open("inlist", "w")
  outfile = open("science.list", "w")
  print "Science data:"
  fname = allscience[0]
  header = pyfits.getheader(datadir + "/" + fname)
  decker = header['deckpos']
  ccdsum = header['ccdsum'].strip()
  for fname in allscience:
    header = pyfits.getheader(datadir + "/" + fname)
    if header['object'] == 'ThAr':
      outfile2 = open("ThAr.list", "w")
      outfile2.write("%s/%s\n" %(datadir, fname))
      outfile2.close()
      allfile.write("%s/%s\n" %(datadir, fname))
      hdulist = pyfits.open(datadir + "/" + fname, mode='update')
      hdulist[0].header['IMAGETYP'] = 'comp'
      hdulist.flush()
      hdulist.close()
    else:
      print "\t", fname
      outfile.write("%s/%s\n" %(datadir, fname))
      allfile.write("%s/%s\n" %(datadir, fname))
    if header['deckpos'] != decker or header['ccdsum'].strip() != ccdsum:
      print "Error! deckpos and ccdsum are not the same in all science observations!"
      sys.exit()
  outfile.close()
  
  
  #Make lists of the appropriate calibration data
  allcalib = [f for f in os.listdir(caldir) if f.startswith("chi")]
  FileDict = defaultdict(list)

  for fname in allcalib:
    header = pyfits.getheader(caldir + "/" + fname)
    otype = header["OBJECT"]
    if header["ccdsum"].strip() == ccdsum and header['deckpos'] == decker and otype != "ThAr":
      if otype == "bias":
        #Change imagetyp keyword from 'bias' to 'zero'
        hdulist = pyfits.open(caldir + "/" + fname, mode='update')
        hdulist[0].header['IMAGETYP'] = 'zero'
        hdulist.flush()
        hdulist.close()
      if otype == 'quartz':
        #Change imagetyp keyword to 'flat'
        hdulist = pyfits.open(caldir + "/" + fname, mode='update')
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
      outfile.write("%s/%s\n" %(caldir, fname))
      allfile.write("%s/%s\n" %(caldir, fname))
    outfile.close()
  allfile.close()

  

    
