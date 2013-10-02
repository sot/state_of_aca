from glob import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from darkbins import x0, x1, dx, bins
from astropy.io import ascii

cals = ascii.read('../warm_pix_from_dc/dark_cals_with_temp.txt')
cals_map = {re.sub(r':', '', cal['day']): cal for cal in cals}
t_ccd_ref = -19.0  # CCD reference temperature


def get_peak(dark):
    y, xb = np.histogram(dark, np.arange(0, 30, 0.5))
    imax = np.argmax(y)
    return np.mean(xb[imax:imax + 2])


def make_darkhist(dateglob='201????', plot=True, outroot=None):
    if plot:
        plt.clf()

    for darkdir in glob('aca_dark_cal/{}'.format(dateglob)):
        date = re.search(r'(\d+)', darkdir).group(1)
        cal = cals_map[date]
        print darkdir, date
        hdus = fits.open(os.path.join(darkdir, 'imd.fits'))
        dark = hdus[0].data.flatten()
        hdus.close()

        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.
        scale = 10 ** ((t_ccd_ref - cal['temp']) / 21.0)

        dark *= scale
        dark += np.random.uniform(-0.5, 0.5, 1024 ** 2)
        peak = get_peak(dark)
        dark -= peak
        print ' t_ccd={:.2f} scale={:.2f} peak={:.2f} zodib={:.2f}'.format(cal['temp'], scale, peak, cal['zodib'])
        y, xb = np.histogram(dark, bins)
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
