from astropy.io import fits as pyfits
import sys


for fname in sys.argv[1:]:
  header = pyfits.getheader(fname)
  #print header.ascardlist()
  print fname, header['DATE']

