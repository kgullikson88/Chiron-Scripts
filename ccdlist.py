import sys

from astropy.io import fits as pyfits


for fname in sys.argv[1:]:
    header = pyfits.getheader(fname)
    hdulist = pyfits.open(fname)
    print fname, header['object'], header['ra'], header['dec'], header['IMAGETYP'], header["ccdsum"], header['deckpos'], \
    header['date-obs'], header['zd'], header['exptime'], len(hdulist)
