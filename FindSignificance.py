import numpy as np
import matplotlib.pyplot as plt
import Bootstrap
import sys


if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  for fname in fileList:
