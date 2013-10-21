"""
Compare dark current histograms for the 2006:329 (T=-15 C) and 2007:069 (T=-20 C)
calibrations.

There is a 5-6% mismatch in the cool-warm region that might point to
possible improvements in the way that dark current is scaled to compensate
for temperature changes.
"""

import os
import re

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import darkbins
from astropy.io import ascii

cals = ascii.read('../warm_pix_from_dc/dark_cals_with_temp.txt')
cals_map = {re.sub(r':', '', cal['day']): cal for cal in cals}
t_ccd_ref = -19.0  # CCD reference temperature


def get_peak(dark):
    y, xb = np.histogram(dark, np.arange(0, 20, 0.1))
    imax = np.argmax(y)
    plt.figure(2)
    x = (xb[:-1] + xb[1:]) / 2.
    plt.plot(x, y)
    peak = np.mean(xb[imax:imax + 2])
    ok = np.abs(x - peak) < 2
    p = np.polyfit(x[ok], y[ok], 2)
    print p
    fit_peak = -p[1] / (2 * p[0])
    plt.plot(x[ok], np.polyval(p, x[ok]))

    return fit_peak


def compare_darkhist():
    plt.figure(1)
    plt.clf()
    bins = darkbins.bins

    for date in ('2006329', '2007069'):
        darkdir = 'aca_dark_cal/{}'.format(date)
        cal = cals_map[date]
        print darkdir, date
        hdus = fits.open(os.path.join(darkdir, 'imd.fits'))
        dark = hdus[0].data.flatten()
        hdus.close()

        dark += np.random.normal(0.0, 0.8, size=1024 ** 2)

        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.
        scale = 10 ** ((t_ccd_ref - cal['temp']) / 21.0)
        dark *= scale

        peak = get_peak(dark)
        dark += (8.0 - peak)
        n100 = np.sum(dark > 100)
        print (' t_ccd={:.2f} scale={:.2f} peak={:.2f} zodib={:.2f} n100={}'
               .format(cal['temp'], scale, peak, cal['zodib'], n100))

        y, xb = np.histogram(dark, bins)
        x = (xb[1:] + xb[:-1]) / 2

        plt.figure(1)
        plt.loglog(x, y)
        plt.loglog(x, y, '.')
        plt.draw()
        plt.show()
