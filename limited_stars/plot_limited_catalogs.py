import os
import numpy as np
import matplotlib.pyplot as plt

import stars_per_obsid

def plot_star_hists():

    guide, acq = stars_per_obsid.stars_per_obsid()

    gf = plt.figure(figsize=(6, 4))
    n, bins, patches = plt.hist([guide[-7], guide[-11], guide[-14]],
                                histtype='bar',
                                align='left',
                                log=True,
                                bins=range(0,8),
                                color=['red', 'green', 'blue'],
                                label=['-7C', '-11C', '-14C'])
    plt.legend(loc='upper left', fontsize=12)
    plt.title('Guide stars per obsid brighter than limiting mag\nat projected temps in 2018')
    plt.xlabel('N guide stars')
    plt.ylabel('N obsids')
    plt.grid()
    plt.tight_layout()

    af = plt.figure(figsize=(6, 4))
    n, bins, patches = plt.hist([acq[-7], acq[-11], acq[-14]],
                                histtype='bar',
                                align='left',
                                log=True,
                                bins=range(0,10),
                                color=['red', 'green', 'blue'],
                                label=['-7C', '-11C', '-14C'])
    plt.legend(loc='upper left', fontsize=12)
    plt.title('Acq stars per obsid brighter than limiting mag\nat projected temps in 2018')
    plt.xlabel('N ACQ stars')
    plt.ylabel('N obsids')
    plt.grid()
    plt.tight_layout()
