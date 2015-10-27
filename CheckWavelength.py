"""
Have the user click on the balmer line
"""

import matplotlib.pyplot as plt 
import HelperFunctions
import logging
import numpy as np

class EventHandler(object):
    def __init__(self, x, ax):
        self.current_x = x
        self.ax = ax
        self.line = ax.plot([x, x], ax.get_ylim(), 'r--')


    def onclick(self, event):
        self.current_x = event.xdata
        self.line[0].remove()
        self.line = self.ax.plot([event.xdata, event.xdata], self.ax.get_ylim(), 'r--')
        plt.draw()


def plot(filename):
    # Read in and combine the data
    orders = HelperFunctions.ReadExtensionFits(filename)
    combined = HelperFunctions.CombineXYpoints(orders)
    good = (combined.x > 475) & (combined.x < 495)
    combined = combined[good]

    guess = combined.y.argmin()
    print('Guess x = {}'.format(combined.x[guess]))

    fig, ax = plt.subplots(1, 1, figsize=(13, 10))
    ax.plot(combined.x, combined.y, 'k-', alpha=0.5)
    ax.set_xlim((combined.x[guess]-2, combined.x[guess]+2))
    xmin = np.searchsorted(combined.x, combined.x[guess] - 2)
    xmax = np.searchsorted(combined.x, combined.x[guess] + 2)
    ymax = min(1.0, combined.y[xmin], combined.y[xmax])
    ymin = ax.get_ylim()[0]
    ax.set_ylim((ymin, ymax))
    ax.set_title(filename)

    # Set up event handler
    eh = EventHandler(combined.x[guess], ax)
    fig.canvas.mpl_connect('button_press_event', eh.onclick)

    plt.show()
    return eh.current_x


if __name__ == '__main__':
    import sys
    with open('wave_log.txt', 'w') as outfile:
        for fname in sys.argv[1:]:
            wave = plot(fname)
            print('User-chosen x = {}'.format(wave))
            outfile.write('{},{}\n'.format(fname, wave))