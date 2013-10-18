from glob import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.io import ascii

import darkbins

cals = ascii.read('../warm_pix_from_dc/dark_cals_with_temp.txt')
cals_map = {re.sub(r':', '', cal['day']): cal for cal in cals}
t_ccd_ref = -19.0  # CCD reference temperature


def get_peak(dark):
    y, xb = np.histogram(dark, np.arange(0, 20, 0.1))
    imax = np.argmax(y)
    x = (xb[:-1] + xb[1:]) / 2.
    peak = np.mean(xb[imax:imax + 2])
    ok = np.abs(x - peak) < 2
    p = np.polyfit(x[ok], y[ok], 2)
    fit_peak = -p[1] / (2 * p[0])

    return fit_peak


def make_darkhist(dateglob='20?????', plot=True, outroot=None, peak_norm=8.0, zodi_norm=False):
    if plot:
        plt.clf()

    darkfiles = glob('aca_dark_cal/{}'.format(dateglob))
    for darkdir in sorted(darkfiles):
        date = re.search(r'(\d+)', darkdir).group(1)
        if date not in cals_map:
            print '{} not in cals map'.format(date)
            continue

        cal = cals_map[date]
        print 'Processing', darkdir, date
        hdus = fits.open(os.path.join(darkdir, 'imd.fits'))
        dark = hdus[0].data.flatten()
        hdus.close()

        # Put in random gaussian noise to smooth things out at the low end
        dark += np.random.normal(0.0, 0.8, size=1024 ** 2)

        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.
        scale = 10 ** ((t_ccd_ref - cal['temp']) / 21.0)
        dark *= scale

        if peak_norm is not None:
            peak = get_peak(dark)
            dark += peak_norm - peak

        if zodi_norm:
            dark -= cal['zodib']

        print (' t_ccd={:.2f} scale={:.2f} peak={:.2f} zodib={:.2f}'
               .format(cal['temp'], scale, peak, cal['zodib']))

        y, xb = np.histogram(dark, darkbins.bins)
        x = (xb[1:] + xb[:-1]) / 2
        if plot:
            plt.loglog(x, y)
            plt.draw()
            plt.show()

        if outroot:
            out = open('{}/{}.dat'.format(outroot, date), 'w')
            for i in range(len(x)):
                print >>out, x[i], y[i], 1 + y[i] * 0.1
            out.close()
