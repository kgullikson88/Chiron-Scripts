import pyfits
import sys


for fname in sys.argv[1:]:
  header = pyfits.getheader(fname)
  print fname, header['object']

