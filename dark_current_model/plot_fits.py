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
        year = float(yeardoy[0:4]) + float(yeardoy[4:7])/365.
        fit = fits[key]
        ipar = list(fit['parnames']).index(par)
        records.append((year, fit['parvals'][ipar], fit['parmins'][ipar], fit['parmaxes'][ipar]))
    return np.rec.fromrecords(records, names=('year', 'val', 'lo', 'hi'))

def make_plots():
    figure(figsize=(5.5, 4))
    clf()
    errorbar(g1['year'], g1['val'], yerr=-g1['lo'])
    title('Cool/warm powerlaw index')
    xlabel('Year')
    ylabel('Gamma 1')
    savefig('gamma1.png')

    clf()
    errorbar(g2['year'], g2['val'], yerr=g2['lo'])
    title('Warm/hot powerlaw index')
    xlabel('Year')
    ylabel('Gamma 2')
    savefig('gamma2.png')

    clf()
    errorbar(amp['year'], amp['val'], yerr=amp['lo'])
    title('Powerlaw normalization')
    xlabel('Year')
    ylabel('Amplitude 1')
    savefig('ampl1.png')

    clf()
    errorbar(xb['year'], xb['val'], yerr=xb['lo'])
    title('Broken powerlaw break point')
    xlabel('Year')
    ylabel('Break point (e-/sec)')
    savefig('break.png')
    
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

# fits = pickle.load(open('fits.pickle'))
fits = pickle.load(open('fits_cash_6000.pickle'))
# fits = pickle.load(open('fits_fixed_6000.pickle'))
g1 = extract('gamma1', fits)
g2 = extract('gamma2', fits)
xb = extract('x_b', fits)
amp = extract('ampl1', fits)

