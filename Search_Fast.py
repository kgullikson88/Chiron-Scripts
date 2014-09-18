import sys

import GenericSearch


#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 466],
              [587.5, 593],
              [627, 634.5],
              [655, 657],  # H alpha
              [485, 487],  #H beta
              [686, 706],
              [716, 742],
              [749.1, 749.45],
              [759, 770],
              [780, 9e9]]

if __name__ == "__main__":
    # Parse command line arguments:
    fileList = []
    extensions = True
    tellurics = False
    trimsize = 100
    for arg in sys.argv[1:]:
        if "-e" in arg:
            extensions = False
        if "-t" in arg:
            tellurics = True  #telluric lines modeled but not removed
        else:
            fileList.append(arg)

    GenericSearch.CompanionSearch(fileList, 
	                          extensions=extensions, 
				  resolution=80000.0, 
				  trimsize=trimsize, 
				  modeldir='/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/')

