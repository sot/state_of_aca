"""
Plot the fraction of time spent not tracking for stars near a certain mag versus
N100 warm fraction.

In [27]: plot_no_trak(10.6, 0.2)
In [26]: plot_no_trak(10.3, 0.2)
In [28]: plot_no_trak(10.0, 0.2)

"""
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import ascii
dat = ascii.read('gui_warm_fit.dat')

if 'gs' not in globals():
    gs = np.load('guide_stats.npy')
    n100 = gs['n100'] / 1024. ** 2
    no_trak = gs['not_tracking_samples'] / gs['n_samples']
    # Clean things up at the low end for a semilog plot
    no_trak_plot = no_trak + 1e-3 + 1e-3 * np.abs(np.random.normal(size=len(no_trak)))


def plot_no_trak(mag=10.6, dmag=0.1):
    ok = np.abs(gs['aoacmag_mean'] - mag) < dmag

    plt.close(1)
    plt.figure(figsize=(6, 4))

    plt.semilogy(n100[ok], no_trak_plot[ok], '.')
    x0, x1 = plt.xlim()
    plt.plot([x0, x1], [0.05, 0.05], '--r', linewidth=2, label='5% no track')
    plt.xlim(0.02, 0.16)
    plt.ylim(0.9e-3, 1.5)
    plt.grid()
    plt.xlabel('N > 100 warm fraction')
    plt.ylabel('Not tracking fraction')
    plt.title('Not tracking vs. N > 100 for {:.1f} +/- {:.1f} mag stars'.format(mag, dmag))
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.savefig('no_trak_{:.1f}.png'.format(mag))
