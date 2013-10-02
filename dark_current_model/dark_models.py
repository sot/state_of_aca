"""Define models to represent the observed ACA CCD dark current distribution.
Provide a function to degrade the CCD for a delta time."""

import numpy as np
from numpy import exp, log, arange

# Define a common fixed binning of dark current distribution
from darkbins import x0, x1, dx, bins as xbins

# Some constants and globals.  Done this way to support sherpa fitting.
# Needs to be re-worked to be nicer.

# Fixed gaussian for smoothing the broken power law
dx = 0.1
sigma = 0.25                            # Gaussian sigma in log space
xg = arange(-2.5 * sigma, 2.5 * sigma, dx, dtype=float)
yg = exp(-0.5 * (xg / sigma) ** 2)
yg /= np.sum(yg)

NPIX = 1024 ** 2

# Fixed
xall = (xbins[:-1] + xbins[1:]) / 2.0
imin = 0
imax = len(xall)


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
    dark = np.append(np.random.lognormal(log(4.5), 0.25, NPIX-nlow),
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

def degrade_ccd(dark, dyear):
    """Degrade CCD (dark current map) for dyear years.  The input 'dark' map is
    updated in place."""
    pars = nompars(dyear)
    darkmodel = smooth_twice_broken_pow(pars, xall)

    ndegrade = 800000
    ipix = np.random.randint(0, NPIX, ndegrade)
    dark[ipix] += np.random.poisson(dyear * 0.7, ndegrade)

    darkran = np.random.poisson(darkmodel)
    for i, n in enumerate(darkran):
        # Generate n log-uniform variates within bin
        if n > 0:
            logdark = np.random.uniform(log(xbins[i]), log(xbins[i+1]), n)
            ipix = np.random.randint(0, NPIX, n)
            dark[ipix] += exp(logdark)

def temp_scalefac(T_ccd):
    """Return the multiplicative scale factor to convert a CCD dark map from
    the nominal -19C temperature to the temperature T.  Based on an overall reduction
    of 0.62 after changing from -15 to -19.
    """
    return exp(log(0.62) / 4.0 * (-19.0 - T_ccd))
