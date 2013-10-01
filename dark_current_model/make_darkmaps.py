"""
Evolve the CCD, generating fake dark current map using the
smoothed-broken-powerlaw spectral model fitting (from fit_evol.py).

First use the dates of actual dark current cals (for direct comparison), then
go once per year.
"""

import dark_models
from Chandra.Time import DateTime
import matplotlib.ticker as ticker
import os
import re
import ParseTable
import numpy as np

def plot_overlays():
    """Plot 1999223 actual, degradation model, 2008273 actual and 2008273 predicted"""
    figure(1, figsize=(7,5))
    clf()
    plot_darkcal('aca_dark_cal/1999223/imd.fits', -10, '-g', label='1999 SSD closed', step=True)
    plot_darkcal('aca_dark_cal/2008273/imd.fits', -19, '-b', label='2008 actual', step=True)
    plot_darkcal('worst/from1999223/2008273.fits', -19, '.b', label='2008 predicted')
    plot_darkcal('worst/from1999223/2019182.fits', -18, '.r', label='2019 predicted')

    pars = dark_models.nompars(dyear=1.0)
    darkmodel = dark_models.smooth_twice_broken_pow(pars, dark_models.xall)
    plot(dark_models.xall, darkmodel, label='Degradation model')

    xlim(1, 1e4)
    ylim(0.5, 3e5)
    xlabel('Dark current (e-/sec)')
    ylabel('Number per bin')
    title('Dark current distributions')
    legend(loc=3)
    savefig('dark_overlays.png')
    
def plot_darkcal(darkfile, Tccd, plotsymcol='-', label=None, step=None):
    """
    Plot a real dark cal, scaled from temp C{Tccd} to -19

    @param darkfile: FITS file containing a dark current map image.
    @param Tccd: CCD temperature (degC)
    @param plotsymcol: Matplotlib symbol/color specifierr
    @param label: Matplotlib plot label
    @param step: Matplotlib plot step specifier

    @return: C{(dark, x, y)}
      - C{dark}: Dark current map (scaled to temperature C{Tccd})
      - C{x, y}: Histogram values plotted
    """
    hdus = pyfits.open(darkfile)
    # Convert to an effective T_ccd = -19
    dark = hdus[0].data.flatten() / dark_models.temp_scalefac(Tccd)
    x, y = plot_dark(dark, label, plotsymcol, step=step)
    return dark, x, y

def plot_dark(dark, label=None, plotsymcol='.', step=False):
    """Plot dark cal data as a histogram.  'dark' is a flat array."""
    dark = dark + np.random.uniform(-0.5, 0.5, len(dark))
    y, x = histogram(dark, dark_models.xbins, new=True)
    if step:
        x, y = make_steps(x, y)
    else:
        x = (x[1:] + x[:-1])/2
    y = np.array(y, dtype=float)
    y[y < 1e-5] = 0.1
    loglog(x, y, plotsymcol, label=label)
    return x, y

def make_steps(x, y):
    """For x as bin edges (n+1 vals) and y as bin vals (n vals), make the corresponding
    arrays that can be plotted as a line."""
    lx = x[:-1]
    rx = x[1:]
    newx = np.array([lx, rx]).flatten(2)
    newy = np.array([y, y]).flatten(2)
    return newx, newy

def plot_warm_frac(thresh, warmpix):
    mxdates = (DateTime(x).mxDateTime for x in alldarkcals['date'])
    yrs = [x.year + x.day_of_year / yeardays for x in mxdates]
    fracs = alldarkcals[thresh] / npix

    figure(1, figsize=(6, 4.5))
    ax = subplot(1, 1, 1)
    plot(yrs, fracs, '.')
    formatter = ticker.FormatStrFormatter('%d')
    ax.xaxis.set_major_formatter(formatter)
    plot(warmpix['year'], np.array(warmpix[thresh]) / npix)
    xlim(2000, 2020)
    ylim(-0.005, 0.251)
    xlabel('Year')
    ylabel('Fraction')
    title('Warm pixel fraction (Nominal and Worst)')
    draw()
    # savefig(thresh + '_' + case + '.png')

npix = 1024.0**2
secs2k = DateTime('2000:001:00:00:00').secs
sec2year = 1 / (86400. * 365.25)
yeardays = 365.25
alldarkcals = ParseTable.parse_table('darkcal_stats.dat')

def make_darkmaps(case='nominal', initial='pristine'):
    """Make dark current maps by degrading the CCD using the best-fit model
    from fit_evol.py"""

    date0 = DateTime('1999-05-23T18:00:00')

    if initial == 'pristine':
        # Start from a synthetic pristine CCD.  Probably no good reason to use this,
        # prefer starting from the pre-SSD-open dark cal on 1999223.
        dark = dark_models.pristine_ccd()
        darkcals = alldarkcals

    elif re.match('from|zero', initial):
        darkmap, year, doy = re.match(r'(\S+)(\d{4})(\d{3})', initial).groups()
        yeardoy = year + ':' + doy
        date0 = DateTime(yeardoy)
        darkcals = alldarkcals[alldarkcals['date'] > yeardoy]

        if darkmap == 'zero':
            dark = dark_models.zero_dark_ccd()
        else:
            ok = alldarkcals['date'] == yeardoy
            tccd = alldarkcals[ok]['tccd'][0]
            scalefac = dark_models.temp_scalefac(tccd)
            print 'Scaling dark cal', yeardoy, 'by', '%.2f' % scalefac
            hdus = pyfits.open(os.path.join('aca_dark_cal', year+doy, 'imd.fits'))
            dark = hdus[0].data.flatten() / scalefac

    warmpix = dict(year=[], n100=[], n200=[], n2000=[])

    print 'Calculating for', case, 'case starting from', initial, 'CCD'

    # First evolve from date0 through last dark cal
    dates = [x['date'] for x in darkcals] + [str(yr) + ':182' for yr in range(2009, 2020)]

    datelast = date0
    for datestr in dates:
        date = DateTime(datestr)
        datemx = date.mxDateTime
        print date.date

        # Evolve (degrade) the CCD dark map.  'dark' is for an effective T=-19 temperature.
        dyear = (datemx - datelast.mxDateTime).days / yeardays
        dark_models.degrade_ccd(dark, dyear)

        # Determine actual or predicted CCD temperature
        try:
            ok = darkcals['date'] == datestr
            tccd = darkcals[ok]['tccd'][0]
        except IndexError:
            if datemx.year > 2014:
                if case == 'nominal':
                    tccd = -18
                else:
                    tccd = -15
            else:
                tccd = -19

        # Calculate dark map at actual temperature and possibly write out to file
        scalefac = dark_models.temp_scalefac(tccd)
        dark_tccd = dark * scalefac         
        outdir = os.path.join(case, initial)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.join(outdir, '%04d%03d.fits' % (datemx.year, datemx.day_of_year))

        if os.path.exists(outfile):
            os.unlink(outfile)
        hdu = pyfits.PrimaryHDU(np.float32(dark_tccd.reshape(1024,1024)))
        hdu.writeto(outfile)

        # Calculate the standard warm/hot pixel stats for prediction
        warmpix['year'].append(datemx.year + datemx.day_of_year / yeardays)
        warmpix['n100'].append(len(where(dark_tccd > 100)[0]))
        warmpix['n200'].append(len(where(dark_tccd > 200)[0]))
        warmpix['n2000'].append(len(where(dark_tccd > 2000)[0]))

        datelast = date

##     out = open('warmpix_' + case + '.dat', 'w')
##     for i, n100 in enumerate(warmpix['n100']):
##         print >>out, warmpix['year'][i], warmpix['n100'][i], warmpix['n200'][i], warmpix['n2000'][i]
##     out.close()

    return warmpix
