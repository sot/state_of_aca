from pylab import *
import re
import pickle
import numpy as np
# from Chandra.Time import DateTime


def extract(par, fits):
    ('parnames', 'parvals', 'parmins', 'parmaxes')
    par = 'sbp.' + par
    records = []
    for key in sorted(fits):
        yeardoy = re.search(r'(\d{7})', key).group(1)
        year = float(yeardoy[0:4]) + float(yeardoy[4:7]) / 365.
        fit = fits[key]
        ipar = list(fit['parnames']).index(par)
        records.append((year, fit['parvals'][ipar], fit['parmins'][ipar], fit['parmaxes'][ipar]))
    return np.rec.fromrecords(records, names=('year', 'val', 'lo', 'hi'))


def extract1(par, fits):
    ('parnames', 'parvals', 'parmins', 'parmaxes')
    par = 'sbp.' + par
    records = []
    yeardoy = re.search(r'(\d{7})', key).group(1)
    year = float(yeardoy[0:4]) + float(yeardoy[4:7]) / 365.
    fit = fits[key]
    ipar = list(fit['parnames']).index(par)
    records.append((year, fit['parvals'][ipar], fit['parmins'][ipar], fit['parmaxes'][ipar]))
    return np.rec.fromrecords(records, names=('year', 'val', 'lo', 'hi'))


def make_plots(g1, g2, amp, xb, rootdir=None):
    if g1 is not None:
        plt.figure(1, figsize=(5.5, 4))
        plt.clf()
        plt.errorbar(g1['year'], g1['val'], yerr=-g1['lo'])
        plt.title('Cool/warm powerlaw index')
        plt.xlabel('Year')
        plt.ylabel('Gamma 1')
        plt.grid(True)
        if rootdir:
            plt.savefig('{}/gamma1.png'.format(rootdir))

    if g2 is not None:
        plt.figure(2, figsize=(5.5, 4))
        plt.clf()
        plt.errorbar(g2['year'], g2['val'], yerr=g2['lo'])
        plt.title('Warm/hot powerlaw index')
        plt.xlabel('Year')
        plt.ylabel('Gamma 2')
        plt.grid(True)
        if rootdir:
            plt.savefig('{}/gamma2.png'.format(rootdir))

    if amp is not None:
        plt.figure(3, figsize=(5.5, 4))
        plt.clf()
        plt.errorbar(amp['year'], amp['val'], yerr=amp['lo'])
        plt.title('Powerlaw normalization')
        plt.xlabel('Year')
        plt.ylabel('Amplitude 1')
        plt.grid(True)
        if rootdir:
            plt.savefig('{}/ampl1.png'.format(rootdir))

    if xb is not None:
        plt.figure(4, figsize=(5.5, 4))
        plt.clf()
        plt.errorbar(xb['year'], xb['val'], yerr=xb['lo'])
        plt.title('Broken powerlaw break point')
        plt.xlabel('Year')
        plt.ylabel('Break point (e-/sec)')
        plt.grid(True)
        if rootdir:
            plt.savefig('{}/break.png'.format(rootdir))


def extrap_amp():
    """Extrapolate dark current amplitude to 15, 20 years MET"""
    figure(1, figsize=(5.5, 4))
    year = amp['year'] - 2000
    yrs = np.array([1, 14.5, 19.5])
    r = polyfit(year, amp['val'], 1)
    amps = polyval(r, yrs)
    print r
    print amps

    clf()
    errorbar(year, amp['val'], yerr=amp['lo'])
    plot(yrs, amps, 'o-')
    title('Powerlaw normalization')
    xlabel('Year - 2000')
    ylabel('Amplitude')
    savefig('ampl_extrap.png')


def get_par_fits(filename='darkhist_zodib/fits.pkl'):
    fits = pickle.load(open(filename))
    try:
        g1 = extract('gamma1', fits)
    except Exception as err:
        print 'g1', err
        g1 = None
    try:
        g2 = extract('gamma2', fits)
    except Exception as err:
        print 'g2', err
        g2 = None
    try:
        xb = extract('x_b', fits)
    except Exception as err:
        print 'xb', err
        xb = None
    try:
        amp = extract('ampl1', fits)
    except Exception as err:
        print 'amp', err
        amp = None
    return g1, g2, amp, xb
