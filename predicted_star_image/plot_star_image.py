"""
Do a by-eye comparison of current limiting mag with 25 year mission values.
The commands below are the "keepers" where the by-eye signal to noise of the
star was similar.

>>> draw_star(t_ccd=-15, mag=10.6, year=2014, figure=3, savefig=True)

>>> draw_star(t_ccd=2.5, mag=8.3, year=2024, figure=3, savefig=True)
>>> draw_star(t_ccd=0, mag=8.5, year=2024, figure=3, savefig=True)
>>> draw_star(t_ccd=-5, mag=9.0, year=2024, figure=3, savefig=True)
>>> draw_star(t_ccd=-10, mag=9.4, year=2024, figure=3, savefig=True)
"""

import sys

import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits

sys.path.insert(0, '../dark_current_model')
import dark_models

if not 'DARKS' in globals():
    DARK0 = dark_models.get_dark_map('2013191')
    STAR = fits.open('obs890_adat41.fits')[1].data
    DARKS = {}
    DARKS[2014] = dark_models.degrade_ccd(DARK0, 0.5, inplace=False)
    for year in range(2015, 1999 + 25 + 1):
        print "Creating year {} image".format(year)
        DARKS[year] = dark_models.degrade_ccd(DARKS[year - 1], 1.0, inplace=False)


# 237 == 1
# 279 == 100
# dt = 42 degC => 10*2 => 21 degC per decade

STAR_IMG = STAR[400]['img_corr'] / 1.7  # e-/sec from a 1.7 sec readout
STAR_MAG = 10.3
SL = slice(400, 440)


def draw_star(year=2014, t_ccd=-16.0, mag=10.3, figure=1, savefig=False):
    year = int(year)

    plt.close(figure)
    plt.figure(figure, figsize=(6, 6))

    img = dark_models.temp_scalefac(t_ccd) * DARKS[year][SL, SL]
    nn = img.shape[0] / 2
    img[nn - 3:nn + 3, nn - 3:nn + 3] += STAR_IMG * 10 ** ((STAR_MAG - mag) * 0.4)

    plt.imshow(img, interpolation='nearest', cmap=cm.afmhot)
    plt.colorbar()
    plt.title('Year={} T_ccd={:.1f} Mag={:.1f}'.format(year, t_ccd, mag))
    if savefig:
        plt.savefig('star_{}_t{:.1f}_mag{:.1f}.png'.format(year, t_ccd, mag))


def draw_star_dark_range():
    plt.figure(figsize=(3.5, 3.5))
    for t in [-15, -10, -5, 0]:
        for mag in [10.6, 10.0]:
            draw_star(t, mag)
            plt.title('T = {} degC, Mag_ACA = {}'.format(t, mag))
            ax = plt.gca()
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            plt.tight_layout()
            plt.savefig('star_img_mag{}_T{}.png'.format(mag, t))


def draw_star_dark2000():
    plt.figure(figsize=(3.5, 3.5))
    t = -20.3
    mag = 10.6

    draw_star(t, mag)
    plt.title('T = -10 degC, Mag_ACA = {}'.format(t, mag))
    ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.tight_layout()
    plt.savefig('star_img_2000_mag{}_T{}.png'.format(mag, t))
