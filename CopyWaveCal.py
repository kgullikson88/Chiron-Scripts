"""
This function takes as input two files. It copies the wavelength
calibration from the first one to the second one. It assumes
that the size of each order is the same
"""
import HelperFunctions
import sys
import warnings
import os
import matplotlib.pyplot as plt

def Copy(inputfile, outputfile):
  print "Copying wavelength calibration from %s ===> %s" %(inputfile, outputfile)
  input_orders = HelperFunctions.ReadExtensionFits(inputfile)
  output_orders = HelperFunctions.ReadExtensionFits(outputfile)

  # Check to make sure they have the same size
  if len(input_orders) != len(output_orders):
    print "Different number of orders in files %s (input) and %s (output)" %(inputfile, outputfile)
    sys.exit()

  # Begin loop over orders
  column_list = []
  for inp, out in zip(input_orders, output_orders):
    # Error checking...
    if inp.size() != out.size():
      sys.exit("Different sized orders for files %s (input) and %s (output)" %(inputfile, outputfile))
    if max(abs(inp.x - out.x)) > 0.2:
      warnings.warn("Wavelength has significantly changed for file %s!" %outputfile)

    plt.plot(inp.x, inp.y, 'k-')
    plt.plot(out.x, out.y, 'r-')
    plt.plot(inp.x, out.y, 'g-')

    columns = {"wavelength": inp.x,
               "flux": out.y,
               "continuum": out.cont,
               "error": out.err}
    column_list.append(columns)
  plt.show()
  return




  HelperFunctions.OutputFitsFileExtensions(column_list, outputfile, outputfile, mode='new')


if __name__ == "__main__":
  if len(sys.argv) == 3:
    Copy(sys.argv[1], sys.argv[2])
  else:
    #Find the appropriate files
    output_files = [f for f in os.listdir("./") if f.endswith("corrected.fits")]
    for out in output_files:
      inp = "%s-0.fits" %(out.split("-0")[0])
      Copy(inp, out)