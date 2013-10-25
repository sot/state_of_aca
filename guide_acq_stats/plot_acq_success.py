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
import os

from astropy.io import ascii

ROOT = os.path.abspath(os.path.dirname(__file__))

if 'acq' not in globals():
    acq = np.load(os.path.join(ROOT, 'acq_stats.npy'))
    n100 = acq['n100'] / 1024. ** 2


def plot_scale_offset_fit():
    """
    Fit the acq_warm_fit values for scale and offset vs. mag.  Use a quadratic fit in log
    space.  The warm fit data are based on sherpa fitting of an empirical model:

      P_acq_success(time_bin) = offset(mag) + scale(mag) * warm_frac(time_bin)

    """
    dat = ascii.read(os.path.join(ROOT, 'acq_warm_fit.dat'))
    dat = dat[1:]   # Skip first point at mag=8.25
    mag = dat['mag_mean']
    mag10 = mag - 10.0
    log_scale = np.log10(dat['scale'])
    log_offset = np.log10(dat['offset'])
    p_scale = np.polyfit(mag10, log_scale, 2)
    p_offset = np.polyfit(mag10, log_offset, 2)
    print('log10(scale) = {:.3f} + {:.3f} * (mag - 10) + {:.3f} * (mag - 10)**2'
          .format(*p_scale[::-1]))
    print('log10(offset) = {:.3f} + {:.3f} * (mag - 10) + {:.3f} * (mag - 10)**2'
          .format(*p_offset[::-1]))

    fit_scale = 10 ** np.polyval(p_scale, mag10)
    fit_offset = 10 ** np.polyval(p_offset, mag10)

    plt.close(1)
    plt.figure(1, figsize=(6, 4))

    plt.errorbar(mag, dat['scale'], yerr=dat['scale_parmin'], fmt='.b')
    plt.plot(mag, fit_scale, '-b', label='scale', linewidth=1.5)

    plt.errorbar(mag, dat['offset'], yerr=dat['offset_parmin'], fmt='.r')
    plt.plot(mag, fit_offset, '-r', label='offset', linewidth=1.5)

    ax = plt.gca()
    ax.set_yscale('log')
    plt.xlabel('ACA mag')
    plt.title('Acquisition warm pixel fit vs. mag')
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.grid()
    plt.show()
    plt.savefig('acq_warm_fit_vs_mag.png')


def plot_acq_success(mag=10.6, dmag=0.1):
    ok = np.abs(acq['mag'] - mag) < dmag
    acqok = acq[ok]
    n100ok = n100[ok]

    plt.close(1)
    plt.figure(figsize=(6, 4))

    bins = np.arange(0.02, 0.1601, 0.02)
    x = (bins[1:] + bins[:-1]) / 2.
    succ = n100ok[acqok['obc_id'] == 'ID']
    fail = n100ok[acqok['obc_id'] != 'ID']
    y_succ, x_bins = np.histogram(succ, bins)
    y_fail, x_bins = np.histogram(fail, bins)
    plt.plot(x, y_fail / (y_succ + y_fail), '-')

    x0, x1 = plt.xlim()
    plt.plot([x0, x1], [0.5, 0.5], '--r', linewidth=2, label='50% success')
    plt.xlim(0.02, 0.16)
    plt.ylim(0., 1.1)
    plt.grid()
    plt.xlabel('N > 100 warm fraction')
    plt.ylabel('Acq success')
    plt.title('Acq success vs. N > 100 for {:.1f} +/- {:.1f} mag stars'.format(mag, dmag))
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()
    plt.savefig('acq_success_{:.1f}.png'.format(mag))
