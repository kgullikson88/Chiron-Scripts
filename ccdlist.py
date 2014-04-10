from astropy.io import fits as pyfits
import sys

for fname in sys.argv[1:]:
  header = pyfits.getheader(fname)
  print fname, header['object'], header['date-obs'], header['zd'], header['exptime']
