import pyfits
import matplotlib.pyplot as plt
import sys
import FitsUtils
import FindContinuum


if __name__ == "__main__":
  fileList = []
  tellurics = False
  normalize = False
  for arg in sys.argv[1:]:
    if "tellcorr" in arg:
      tellurics = True
    elif "-norm" in arg:
      normalize = True
    else:
      fileList.append(arg)

  for fnum, fname in enumerate(fileList):
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    print fname, len(orders)
    plt.figure(fnum)
    if tellurics:
      model = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="model")
    for i, order in enumerate(orders):
      order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=3, highreject=3)
      if tellurics:
        plt.plot(order.x, order.y/order.cont, 'k-')
        plt.plot(order.x, model[i].y, 'r-')
      else:
        if normalize:
          plt.plot(order.x, order.y/order.cont)
          plt.text(order.x.mean(), 1.1, str(i+1))
        else:
          plt.plot(order.x, order.y)
          plt.plot(order.x, order.cont)
      
  plt.show()
