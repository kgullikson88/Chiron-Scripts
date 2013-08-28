import pyfits
import numpy
from numpy.polynomial import chebyshev
import DataStructures
#import Units
from astropy import units, constants
import FindContinuum
import readmultispec as multispec
import subprocess


#Make a dictionary for converting from standard polynomial to chebyshev
#   For now, only go up to order 5
#standard_to_cheb = {0: [[1,],],
#                    1: [[1,0], [0,1]],
#                    2: [[1,0,0.5], [0,1,0], [0,0,0.5]],
#                    3: [[1,0,0.5,0], [0,1,0,0.75], [0,0,0.5,0], [0,0,0,0.25]], 
#                    4: [[1,0,0.5,0,0.625], [0,1,0,0.75,0], [0,0,0.5,0,-0.5], [0,0,0,0.25,0], [0,0,0,0,0.125]], 
#                    5: [[1,0,0.5,0,0.625, 0], [0,1,0,0.75,0,0.625], [0,0,0.5,0,-0.5,0], [0,0,0,0.25,0,5./16.], [0,0,0,0,0.125,0], [0,0,0,0,0,1./16.]]}

def Chebyshev(pvals, coefficients):
  order = coefficients.size
  pmid = pvals[pvals.size/2]
  prange = pvals[-1] - pvals[0]
  nvals = (pvals - pmid)/(prange/2.0)
  wave = numpy.zeros(pvals.size)
  x = []
  x.append(numpy.ones(pvals.size))
  x.append(nvals)
  for i in range(2,order):
    x.append(2.0*nvals*x[i-1] - x[i-2])
  for i in range(order):
    wave = wave + coefficients[i]*x[i]
  return wave
  
  
  
def GetChebyshevCoeffs(data, pvals, order=5):
  pmid = pvals[pvals.size/2]
  prange = pvals[-1] - pvals[0]
  nvals = (pvals - pmid)/(prange/2.0)
  
  return chebyshev.chebfit(nvals, data.x, order)
  


"""
  -Function to read in a fits file of echelle data, and turn it into xypoint structures
  -if extensions is True, then it will assume that each echelle order is in a different fits extension
    -Assumes the fits extensions are all binary tables produced by pyfits, which are basically record arrays
    -In this case, x, y, cont, and errors are used. x and y MUST be given, and must be strings
       corresponding to the name of the record array field
    -if cont is given (again, must be a string), it specifies the field which holds the continuum information
    -if errors are given (a string), it specifies the field which holds the errors/sigma information
    
  -if extensions is False (the default), it assumes the data is in multispec format
    -In this case, errors should be given as an integer (not a string as above) specifying which band number
      (in C-style numbering, which starts at 0) the errors/sigma array is.
    -If there is more than one band, the user MUST give the errors keyword or the code will crash!
      -I should probably fix this at some point...
"""
def MakeXYpoints(datafile, extensions=False, x=None, y=None, cont=None, errors=None):
  print "Reading in file %s: " %datafile

  orders = []
  if extensions:
    #This means the data is in fits extensions, with one order per extension
    #At least x and y should be given (and should be strings to identify the field in the table record array)
    if type(x) != str:
      x = raw_input("Give name of the field which contains the x array: ")
    if type(y) != str:
      y = raw_input("Give name of the field which contains the y array: ")
    hdulist = pyfits.open(datafile)
    if cont == None:
      if errors == None:
        for i in range(1,len(hdulist)):
          data = hdulist[i].data
          xypt = DataStructures.xypoint(x=data.field(x), y=data.field(y))
          orders.append(xypt)
      else:
        if type(errors) != str:
          errors = raw_input("Give name of the field which contains the errors/sigma array: ")
        for i in range(1,len(hdulist)):
          data = hdulist[i].data
          xypt = DataStructures.xypoint(x=data.field(x), y=data.field(y), err=data.field(errors))
          orders.append(xypt)
    elif type(cont) == str:
      if errors == None:
        for i in range(1,len(hdulist)):
          data = hdulist[i].data
          xypt = DataStructures.xypoint(x=data.field(x), y=data.field(y), cont=data.field(cont))
          orders.append(xypt)
      else:
        if type(errors) != str:
          errors = raw_input("Give name of the field which contains the errors/sigma array: ")
        for i in range(1,len(hdulist)):
          data = hdulist[i].data
          xypt = DataStructures.xypoint(x=data.field(x), y=data.field(y), cont=data.field(cont), err=data.field(errors))
          orders.append(xypt)
    hdulist.close()

  else:
    #Data is in a big array.
    hdulist = pyfits.open(datafile)
    data = hdulist[0].data
    hdulist.close()
    numorders = data.shape[0]
    for i in range(numorders):
      wave = data[...,...,0][i]
      flux = data[...,...,1][i]
      error = numpy.ones(wave.size)*1e9
      error[flux > 0] = numpy.sqrt(flux[flux>0])
      orders.append(DataStructures.xypoint(x=wave, y=flux, err=error))
    
    
  return orders
  
  
#Function to output a fits file with the same format
#as the template function. Only implementing Chebyshev for now...
def OutputFitsFile(template, orders, func_order=None, outfilename=None, errors=False):
  hdulist = pyfits.open(template)
  header = hdulist[0].header
  
  #First, get what the current wavelength calibration is from the header
  try:
    ltv = header['ltv1']
  except KeyError:
    ltv = 0.0
  try:
    ltm = header['ltm1_1']
  except KeyError:
    ltm = 1.0

  #Make list of the wavelength data for each order
  string = ""
  wave_factor = 1.0   #factor to multiply wavelengths by to get them in nanometers
  for key in sorted(header.keys()):
    if "WAT2" in key:
      string = string + header[key]
      #print key, header[key]
      for i in range(len(header[key]), 68):
        string = string + " "
    elif "WAT1" in key:
      #Get wavelength units
      if "label=Wavelength"  in header[key] and "units" in header[key]:
        units = header[key].split("units=")[-1]
        if units == "angstroms" or units == "Angstroms":
          header.update(key,"wtype=multispec label=Wavelength units=nanometers")
        
  waveinfo = string.split(" spec")
  output_string = waveinfo[0]
  for info, i in zip(waveinfo[1:], range(1,len(waveinfo))):
    output_string = output_string + " spec" + info.split('"')[0] + '"'
    segments = info.split('"')[1].split()
    for segment in segments[:10]:
      output_string = output_string + segment + " "
    size = int(float(segments[5]))
    z = float(segments[6])  #Dopplar correction: 1/(1+z)
    
    xypt = DataStructures.xypoint(size)
    pvals = (numpy.arange(size) - ltv)/ltm

    wlen0 = 0.0
    func_type = int(float(segments[11]))
    if func_type == 1:
      #Chebyshev
      if func_order == None:
        func_order = int(float(segments[12]))
      pmin = int(float(segments[13]))
      pmax = int(float(segments[14]))
      coefficients = GetChebyshevCoeffs(orders[i-1], pvals, func_order-1)

      #Put stuff in output string
      output_string = output_string + ("%.10g " %wlen0).replace("e","E")
      output_string = output_string + segments[11] + " "
      output_string = output_string + str(func_order) + " "
      output_string = output_string + segments[13] + " "
      output_string = output_string + segments[14] + " "
      for coeff in coefficients:
        output_string = output_string + ("%.14g " %coeff).replace("e","E")
     
    else:
      print "Sorry! This function type is not currently implemented!"

    output_string = output_string[:-1] + '"'
  
  #Split output string so it isn't too many characters
  size = 68
  leftindex=0
  rightindex = size
  keyindex = 1
  while rightindex < len(output_string):
    key = "WAT2_"
    if len(str(keyindex)) == 1:
      key = key + "00"
    elif len(str(keyindex)) == 2:
      key = key + "0"
    key = key+str(keyindex)
    
    header.update(key, output_string[leftindex:rightindex])
    leftindex = rightindex
    rightindex = min(leftindex+size, len(output_string))
    keyindex += 1
  
  #Do the last one outside of the loop
  key = "WAT2_"
  if len(str(keyindex)) == 1:
    key = key + "00"
  elif len(str(keyindex)) == 2:
    key = key + "0"
  key = key+str(keyindex)
  #print key, output_string[leftindex:rightindex]
  header.update(key, output_string[leftindex:rightindex])
  
  #Delete any WAT2 keys that were not used
  for key in sorted(header.keys()):
    if "WAT2" in key:
      if int(key[-3:]) > keyindex:
        del header[key]
    
  hdulist[0].header = header

  #Now, update the data (much easier)
  if errors == False:
    numorders = hdulist[0].data.shape[0]
  else:
    numorders = hdulist[0].data.shape[1]
  for i in range(numorders):
    if errors == False:
      hdulist[0].data[i] = orders[i].y
    else:
      hdulist[0].data[0][i] = orders[i].y
      hdulist[0].data[errors][i] = orders[i].err

  if outfilename == None:
    if "-" in template:
      i = int(template.split("-")[-1].split(".fits")[0])
      outfilename = template.split("-")[0] + "-" + str(i+1) + ".fits"
    else:
      outfilename = template.split(".fits")[0] + "-0.fits"
  
  print "Outputting to ", outfilename
  try:
    hdulist.writeto(outfilename)
  except IOError:
    cmd = subprocess.check_call("rm "+outfilename, shell=True)
    hdulist.writeto(outfilename)
  hdulist.close()
  
  
def GetOrderWaveCalPars(header):
  string = ""
  for key in sorted(header.keys()):
    if "WAT2" in key:
      string = string + header[key]
      #print key, header[key]
      for i in range(len(header[key]), 68):
        string = string + " "
        
  return string.split(" spec")
  
  
def CopyWaveCal(copyfrom, copyto, order=None, scale=1.0):
  """
    This function will copy the wavelength calibration of one or more orders and put it in another file. 
    copyfrom and copyto should both be filenames
    The order keyword can be:
       1: left as None, which will copy all of the orders
       2: an integer, which will only copy that order (number in fortran style!)
       3: a list of integers, in which the calibration will be copied for each one
     scale is a factor to apply to the wavelengths (useful for converting between units)
  """
  
  #First, read in the files
  hdulistfrom = pyfits.open(copyfrom)
  hdulistto = pyfits.open(copyto)
  headerfrom = hdulistfrom[0].header
  headerto = hdulistto[0].header
  
  numorders = hdulistfrom[0].data.shape[0]
  if numorders != hdulistto[0].data.shape[0]:
    print "Warning in CopyWaveCal! Files have different number of orders!"
    numorders = min(numorders, hdulistto[0].data.shape[0])
  
  if type(order) == int:
    order = [order,]
  elif order == None:
    order = range(1,numorders+1)
    
  #Make sure we don't try to do too many orders
  order = sorted(order)
  while len(order) > 0 and order[-1] > numorders:
    order.pop()
  
  #Now, get all the orders from both files
  ordersfrom = GetOrderWaveCalPars(headerfrom)
  ordersto = GetOrderWaveCalPars(headerto)

  
  #Do the copy
  for i in order:
    print "Copying order ", i
    print numorders
    ordersto[i] = ordersfrom[i]
    
    if scale != 1.0:
      #Now we have to get more complicated... sigh
      #Only tested with chebyshev, should work with legendre. Maybe others?
      segments = ordersto[i].split()
      
      #Coefficients start in segment 16
      for j in range(17,len(segments)):
        coeff = segments[j]
        print coeff
        coeff = "%.14g" %(float(coeff.replace("E","e").strip('"'))*scale)
        print coeff
        segments[j] = coeff
        
      #re-save as a string
      ordersto[i] = segments[0]
      for segment in segments[1:]:
        ordersto[i] = ordersto[i] + " " + segment
      ordersto[i] = ordersto[i] + '"'
  
  #Convert the orders into one long string
  output_string = ordersfrom[0]
  for segment in ordersto[1:]:
    output_string = output_string + " spec" + segment
  
  #Update the header
  size = 68
  leftindex=0
  rightindex = size
  keyindex = 1
  while rightindex < len(output_string):
    key = "WAT2_"
    if len(str(keyindex)) == 1:
      key = key + "00"
    elif len(str(keyindex)) == 2:
      key = key + "0"
    key = key+str(keyindex)
    
    print key, output_string[leftindex:rightindex]
    headerto.update(key, output_string[leftindex:rightindex])
    leftindex = rightindex
    rightindex = min(leftindex+size, len(output_string))
    keyindex += 1
  
  #Do the last one outside of the loop
  key = "WAT2_"
  if len(str(keyindex)) == 1:
    key = key + "00"
  elif len(str(keyindex)) == 2:
    key = key + "0"
  key = key+str(keyindex)
  print key, output_string[leftindex:rightindex]
  headerto.update(key, output_string[leftindex:rightindex])
  
  #Delete any WAT2 keys that were not used
  for key in sorted(headerto.keys()):
    if "WAT2" in key:
      if int(key[-3:]) > keyindex:
        del headerto[key]
    
  hdulistto[0].header = headerto
  
  #Save with new filename
  if "-" in copyto:
    i = int(copyto.split("-")[-1].split(".fits")[0])
    outfilename = copyto.split("-")[0] + "-" + str(i+1) + ".fits"
  else:
    outfilename = copyto.split(".fits")[0] + "-0.fits"
  print "Outputting to ", outfilename
  try:
    hdulistto.writeto(outfilename)
  except IOError:
    cmd = subprocess.check_call("rm "+outfilename, shell=True)
    hdulistto.writeto(outfilename)
  hdulistto.close()
  hdulistfrom.close()
  
  return outfilename
  



"""
  Function to output a fits file
  column_dict is a dictionary where the key is the name of the column
     and the value is a numpy array with the data. Example of a column
     would be the wavelength or flux at each pixel
  template is the filename of the template fits file. The header will
     be taken from this file and used as the main header
  mode determines how the outputted file is made. Append will just add
     a fits extension to the existing file (and then save it as outfilename)
     "new" mode will create a new fits file. 
     header_info takes a list of lists. Each sub-list should have size 2 where the first element is the name of the new keyword, and the second element is the corresponding value. A 3rd element may be added as a comment
"""
def OutputFitsFileExtensions(column_dicts, template, outfilename, mode="append", headers_info=[]):
  #Get header from template. Use this in the new file
  if mode == "new":
    header = pyfits.getheader(template)
    
  if not isinstance(column_dicts, list):
    column_dicts = [column_dicts, ]
  if len(headers_info) < len(column_dicts):
    for i in range(len(column_dicts) - len(headers_info)):
      headers_info.append([])

  if mode == "append":
    hdulist = pyfits.open(template)
  elif mode == "new":
    hdulist = pyfits.open(template)
    for i in range(len(hdulist)-1):
      hdulist.pop()
      
  for i in range(len(column_dicts)):
    column_dict = column_dicts[i]
    header_info = headers_info[i]
    columns = []
    for key in column_dict.keys():
      columns.append(pyfits.Column(name=key, format="D", array=column_dict[key]))
    cols = pyfits.ColDefs(columns)
    tablehdu = pyfits.new_table(cols)
  
    #Add keywords to extension header
    num_keywords = len(header_info)
    header = tablehdu.header
    for i in range(num_keywords):
      info = header_info[i]
      if len(info) > 2:
        header.update(info[0], info[1], info[2])
      elif len(info) == 2:
        header.update(info[0], info[1])

    hdulist.append(tablehdu)

      
  hdulist.writeto(outfilename, clobber=True, output_verify='ignore')
  hdulist.close()

  
  
if __name__ == "__main__":
  hdulist = pyfits.open("output-1.fits")
  orders = MakeXYpoints(hdulist[0].header, hdulist[0].data)
  OutputFitsFile("output-1.fits", orders)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
