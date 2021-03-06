"""Define models to represent the observed ACA CCD dark current distribution.
Provide a function to degrade the CCD for a delta time."""
import os
import re

import numpy as np
from numpy import exp, log, arange

import Ska.Numpy
from Chandra.Time import DateTime
from astropy.io import ascii

# Define a common fixed binning of dark current distribution
import darkbins

# Some constants and globals.  Done this way to support sherpa fitting.
# Needs to be re-worked to be nicer.

# Fixed gaussian for smoothing the broken power law
dx = 0.1
sigma = 0.30                            # Gaussian sigma in log space
xg = arange(-2.5 * sigma, 2.5 * sigma, dx, dtype=float)
yg = exp(-0.5 * (xg / sigma) ** 2)
yg /= np.sum(yg)

NPIX = 1024 ** 2

# Fixed
xbins = darkbins.bins
xall = darkbins.bin_centers
imin = 0
imax = len(xall)

basedir = os.path.dirname(__file__)

cals = ascii.read(os.path.join(basedir, '../warm_pix_from_dc/dark_cals_with_temp.txt'))
cals_map = {re.sub(r':', '', cal['day']): cal for cal in cals}


def get_dark_map(date, scale_to_m19=True):
    """
    Read a dark current map and (by default) scale values to the canonical
    CCD temperature of -19 C.

    Returns a 1024**2 float array (1-D).
    """
    from astropy.io import fits
    hdus = fits.open(os.path.join(basedir, 'aca_dark_cal', date, 'imd.fits'))
    dark = hdus[0].data
    hdus.close()

    if scale_to_m19:
        # Scale factor to adjust data to an effective temperature of t_ccd_ref.
        # For t_ccds warmer than t_ccd_ref this scale factor is < 1, i.e. the
        # observed dark current is made smaller to match what it would be at the
        # lower reference temperature.
        scale = 1.0 / temp_scalefac(cals_map[date]['temp'])
        dark *= scale

    return dark


def integrate_from(y, x0, dx=1.0):
    """
    Return the integral of y from x0 to inf (aka max(xall)).
    """
    import scipy.interpolate
    import scipy.integrate
    x_samp = np.arange(x0, xall[-1], dx)
    y_samp = scipy.interpolate.interp1d(xall, y)(x_samp)
    result = scipy.integrate.simps(y_samp, x_samp)
    return result / NPIX


def broken_pow(pars, x):
    """Broken power-law.  Pars are same as bpl1d:
    1: gamma1
    2: gamma2
    3: x_b (break point)
    4: x_r (normalization reference point)
    5: ampl1"""
    (gamma1, gamma2, x_b, x_r, ampl1) = pars
    ampl2 = ampl1 * (x_b / x_r) ** (gamma2 - gamma1)
    ok = x > x_b
    y = ampl1 * (x / x_r) ** (-gamma1)
    y[ok] = ampl2 * (x[ok] / x_r) ** (-gamma2)
    return y


def smooth_broken_pow(pars, x):
    """Smoothed broken power-law.  Pars are same as bpl1d (NOT + gaussian sigma):
    1: gamma1
    2: gamma2
    3: x_b (break point)
    4: x_r (normalization reference point)
    5: ampl1
    #   NOT 6: sigma (bins)"""
    (gamma1, gamma2, x_b, x_r, ampl1) = pars
    ampl2 = ampl1 * (x_b / x_r) ** (gamma2 - gamma1)
    ok = xall > x_b
    y = ampl1 * (xall / x_r) ** (-gamma1)
    y[ok] = ampl2 * (xall[ok] / x_r) ** (-gamma2)
    imin = np.searchsorted(xall, x[0] - 1e-3)
    imax = np.searchsorted(xall, x[-1] + 1e-3)
    return np.convolve(y, yg, mode='same')[imin:imax]


def smooth_twice_broken_pow(pars, x):
    """Smoothed broken power-law.  Pars are same as bpl1d (NOT + gaussian sigma):
    1: gamma1
    2: gamma2
    3: x_b (break point)
    4: x_r (normalization reference point)
    5: ampl1
    #   NOT 6: sigma (bins)"""
    gamma0 = -4.0
    x_b0 = 5

    (gamma1, gamma2, x_b, x_r, ampl1) = pars
    y = ampl1 * (x / x_r) ** (-gamma1)

    ok = x > x_b
    ampl2 = ampl1 * (x_b / x_r) ** (gamma2 - gamma1)
    y[ok] = ampl2 * (x[ok] / x_r) ** (-gamma2)

    ok = x < x_b0
    ampl0 = ampl1 * (x_b0 / x_r) ** (gamma0 - gamma1)
    y[ok] = ampl0 * (x[ok] / x_r) ** (-gamma0)

    return np.convolve(y, yg, mode='same')


def zero_dark_ccd():
    return np.zeros(NPIX)


def pristine_ccd():
    """Generate and return a 'pristine' (zero-year) dark current map.  This was
    empirically derived."""
    nlow = 350000
    dark = np.append(np.random.lognormal(log(4.5), 0.25, NPIX - nlow),
                     np.random.lognormal(log(2.5), 0.45, nlow))
    return dark


def nompars(dyear):
    """Return nominal degradation parameters for dyear interval"""
    gamma1 = 0.14
    gamma2 = 3.15
    x_b = 125.
    x_r = 155.
    ampl1 = 1032. * dyear
    return (gamma1, gamma2, x_b, x_r, ampl1)


def degrade_ccd(dark, dyear, inplace=True):
    """Degrade CCD (dark current map) for dyear years.  The input 'dark' map is
    updated in place if inplace is True."""
    pars = nompars(dyear)
    darkmodel = smooth_twice_broken_pow(pars, xall)

    if not inplace:
        dark_shape = dark.shape
        dark = dark.copy()

    # Work with a flattened reference (does not copy)
    dark = dark.ravel()

    ndegrade = 800000
    ipix = np.random.randint(0, NPIX, ndegrade)
    dark[ipix] += np.random.poisson(dyear * 0.7, ndegrade)

    darkran = np.random.poisson(darkmodel)
    for i, n in enumerate(darkran):
        # Generate n log-uniform variates within bin
        if n > 0:
            logdark = np.random.uniform(log(xbins[i]), log(xbins[i + 1]), n)
            ipix = np.random.randint(0, NPIX, n)
            dark[ipix] += exp(logdark)

    if not inplace:
        return dark.reshape(dark_shape)


def temp_scalefac(T_ccd):
    """Return the multiplicative scale factor to convert a CCD dark map from
    the nominal -19C temperature to the temperature T.  Based on best global fit for
    dark current model in plot_predicted_warmpix.py.  Previous value was 0.62 instead
    of 0.70.
    """
    return exp(log(0.70) / 4.0 * (-19.0 - T_ccd))


def as_array(vals):
    if np.array(vals).ndim == 0:
        is_scalar = True
        vals = np.array([vals])
    else:
        is_scalar = False

    vals = np.array(vals)
    return vals, is_scalar


def get_sbp_pars(dates):
    """
    Return smooth broken powerlaw parameters at ``date``.  This is based on the
    sbp fits for the darkhist_peaknorm histograms, with parameters derived from
    by-hand inspection of fit trending.  See NOTES.
    """
    dates, is_scalar = as_array(dates)
    n_dates = len(dates)

    years = DateTime(dates).frac_year

    ones = np.ones(n_dates)
    g1 = 0.05 * ones
    g2 = 3.15 * ones
    x_r = 50.0 * ones

    ampl = (years - 2000.0) * 1390.2 + 1666.4

    bp_years = np.array([1999.0, 2000.9, 2003.5, 2007.0, 2007.01, 2011.5, 2011.51, 2013.7])
    bp_vals = np.array([125.000, 125.00, 110.00, 117.80, 111.50, 125.00, 108.80, 115.00])
    bp = Ska.Numpy.interpolate(bp_vals, bp_years, years, method='linear')

    if is_scalar:
        g1 = g1[0]
        g2 = g2[0]
        ampl = ampl[0]
        bp = bp[0]
        x_r = x_r[0]

    return g1, g2, bp, x_r, ampl


def get_darkhist(date='2013:001', T_ccd=-19.0):
    """
    Get model dark current histogram at the given ``date`` and ``T_ccd``.

    Returns bins, darkhist
    """
    pars = get_sbp_pars(date)
    x = darkbins.bin_centers
    y = smooth_broken_pow(pars, x)

    scale = temp_scalefac(T_ccd)
    xbins = darkbins.bins * scale
    x = x * scale

    return x, xbins, y


def get_warm_fracs(warm_thresholds, date='2013:001', T_ccd=-19.0):
    x, xbins, y = get_darkhist(date, T_ccd)
    warm_thresholds, is_scalar = as_array(warm_thresholds)

    warmpixes = []
    for warm_threshold in warm_thresholds:
        # First get the full bins to right of warm_threshold
        ii = np.searchsorted(xbins, warm_threshold)
        warmpix = np.sum(y[ii:])
        lx = np.log(warm_threshold)
        lx0 = np.log(xbins[ii - 1])
        lx1 = np.log(xbins[ii])
        ly0 = np.log(y[ii - 1])
        ly1 = np.log(y[ii])
        m = (ly1 - ly0) / (lx1 - lx0)
        partial_bin = y[ii] * (lx1 ** m - lx ** m) / (lx1 ** m - lx0 ** m)
        # print ii, x[ii], xbins[ii - 1], xbins[ii], y[ii], partial_bin
        warmpix += partial_bin
        warmpixes.append(warmpix)

    if is_scalar:
        out = warmpixes[0]
    else:
        out = np.array(warmpixes)

    return out / (1024.0 ** 2)


def ccd_temperature_model(t_ccd_max, start='2014:001', stop='2018:001', n=20):
    """
    Canonical model of future CCD temperature which is increasing linearly
    from current observed maxes until it gets clipped at ``t_ccd_max``.
    """
    start = DateTime(start)
    stop = DateTime(stop)
    x0 = 2014.0
    y0 = -15.0
    x1 = 2016.0
    y1 = -11.0
    x = np.linspace(start.frac_year, stop.frac_year, n)
    y = (x - x0) / (x1 - x0) * (y1 - y0) + y0
    y = y.clip(None, t_ccd_max)

    return DateTime(x, format='frac_year'), y
