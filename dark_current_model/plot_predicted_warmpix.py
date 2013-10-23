"""
Plot the warm pixel fraction from the 2013 baseline dark current model vs. observed
fractions from dark current calibrations.

Note that the N=1000 prediction doesn't match observation that well (off by a scale
factor).  This is indicative that the pure powerlaw doesn't hold out that far, but
it's not a big issue here.
"""

from glob import glob
import os
import re
from itertools import izip

import numpy as np

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from Ska.Matplotlib import plot_cxctime
from Chandra.Time import DateTime

import dark_models


cals = ascii.read('../warm_pix_from_dc/dark_cals_with_temp.txt')
cals_map = {re.sub(r':', '', cal['day']): cal for cal in cals}
t_ccd_ref = -19.0  # CCD reference temperature


def get_warm_frac(dateglob='20?????', outroot=None, warm_thresholds=None):
    """
    Get the warm pixel fraction from the 2013 baseline dark current model.
    """
    if warm_thresholds is None:
        warm_thresholds = [50., 100., 200.]

    darkfiles = glob('aca_dark_cal/{}'.format(dateglob))
    pred_warm_fracs = {x: list() for x in warm_thresholds}
    obs_warm_fracs = {x: list() for x in warm_thresholds}
    dates = []

    for darkdir in sorted(darkfiles):
        date = re.search(r'(\d+)', darkdir).group(1)
        if date not in cals_map:
            print '{} not in cals map'.format(date)
            continue

        print 'Processing', darkdir, date
        hdus = fits.open(os.path.join(darkdir, 'imd.fits'))
        dark = hdus[0].data.flatten()
        hdus.close()

        cal = cals_map[date]
        warm_fracs = dark_models.get_warm_fracs(warm_thresholds, date=cal['day'], T_ccd=cal['temp'])

        dates.append(DateTime(cal['day']).secs)
        for warm_threshold, warm_frac in zip(warm_thresholds, warm_fracs):
            pred_warm_fracs[warm_threshold].append(warm_frac)
            obs_warm_fracs[warm_threshold].append(np.sum(dark > warm_threshold) / 1024.0 ** 2)

    return dates, pred_warm_fracs, obs_warm_fracs


def plot_warm_frac():
    """
    Plot the warm pixel fraction from the 2013 baseline dark current model vs. observed
    fractions from dark current calibrations.

    Note that the N=1000 prediction doesn't match observation that well (off by a scale
    factor).  This is indicative that the pure powerlaw doesn't hold out that far, but
    it's not a big issue here.
    """
    warm_thresholds = [50., 100., 200.]
    dates, pred_warm_fracs, obs_warm_fracs = get_warm_frac(warm_thresholds=warm_thresholds)

    plt.close(1)
    plt.figure(1, figsize=(6, 4))
    for color, warm_threshold in zip(['b', 'g', 'r', 'k', 'c', 'm', 'y'], warm_thresholds):
        plot_cxctime(dates, obs_warm_fracs[warm_threshold], color + '.')
        plot_cxctime(dates, pred_warm_fracs[warm_threshold], color + '-',
                     label='N > {:.0f} e-/sec'.format(warm_threshold))

    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.03
    plt.xlim(x0 - dx, x1 + dx)
    plt.grid()
    plt.legend(loc='best', fontsize=12)
    plt.title('Warm pixel fraction: model (line), observed (dots) vs. time', fontsize=12)
    plt.ylabel('Warm pixel fraction')
    plt.tight_layout()
    plt.savefig('warm_frac_model.png')


def plot_warm_frac_future():
    """
    Plot the future warm pixel fraction from the 2013 baseline dark current model
    for the N100 case.  Include observed fractions from dark current calibrations.
    """
    warm_threshold = 100.
    dates, pred_warm_fracs, obs_warm_fracs = get_warm_frac(warm_thresholds=[warm_threshold])

    plt.close(1)
    plt.figure(1, figsize=(6, 4))
    color = 'g'
    plot_cxctime(dates, obs_warm_fracs[warm_threshold], color + '.')
    plot_cxctime(dates, pred_warm_fracs[warm_threshold], color + '-')

    for t_ccd_max, offset, color in izip((-14, -11, -7), (0, 0.005, 0.01), ('r', 'g', 'b')):
        dates, t_ccds = dark_models.ccd_temperature_model(t_ccd_max, start='2014:001',
                                                          stop='2018:001', n=20)
        warm_fracs = []
        for time, t_ccd in izip(dates.secs, t_ccds):
            warm_frac = dark_models.get_warm_fracs(warm_threshold, date=time, T_ccd=t_ccd)
            warm_fracs.append(warm_frac + offset)
        plot_cxctime(dates.secs, warm_fracs, color + '-',
                     label='T_ccd < {} C'.format(t_ccd_max))

    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.03
    plt.xlim(x0 - dx, x1 + dx)
    plt.grid()
    # plt.legend(loc='best', fontsize=12)
    plt.title('Warm pixel fraction: model (line), observed (dots) vs. time', fontsize=12)
    plt.ylabel('N > 100 warm pixel fraction')
    plt.tight_layout()
    plt.savefig('warm_frac_future.png')
