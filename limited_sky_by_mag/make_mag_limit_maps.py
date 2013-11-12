from __future__ import division

import cPickle as pickle
import os
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from astropy.table import Table, Column
import healpy

import agasc

NSIDE = 64
limits = [10.13, 10.02, 9.92, 10.38, 10.29, 10.2, 10.0, 10.1, 10.3]


MAX_RADIUS = 0.7836  # same area (1.92 deg**2) as CCD (5000 arcsec on a side).

if 'stars_list' not in globals() and os.path.exists('ras_decs_stars_list.pkl'):
    print 'Reading ras_decs_stars_list.pkl'
    ras, decs, stars_list = pickle.load(open('ras_decs_stars_list.pkl', 'rb'))

    print 'Filtering stars to max radius of {}'.format(MAX_RADIUS)
    # Filter down to MAX_RADIUS
    new_stars_list = []
    for stars in stars_list:
        ok = stars['radius'] < MAX_RADIUS
        new_stars_list.append(stars[ok])
    stars_list = new_stars_list


def make_ras_decs_stars_list():
    """
    Messed up a bit in make_ra_dec_stars list.  Should have added radii and pickled
    a more sensible object.  Fix this here.
    """
    ra_dec_stars_list = pickle.load(open('ra_dec_stars_list.pkl', 'rb'))
    ras = np.array([x[0] for x in ra_dec_stars_list])
    decs = np.array([x[1] for x in ra_dec_stars_list])
    stars_list = []
    for ra, dec, stars in ra_dec_stars_list:
        sys.stdout.write('{:9.3f} {:9.3f}\r'.format(ra, dec))
        radii = agasc.sphere_dist(ra, dec, stars['RA'], stars['DEC'])
        stars = Table(stars)
        stars.add_column(Column(name='radius', data=radii))
        stars_list.append(np.array(stars))
    sys.stdout.write('\n\n')
    pickle.dump([ras, decs, stars_list], open('ras_decs_stars_list.pkl', 'rb'), protocol=-1)


def plot_good_guide_catalogs(mag_limit=10.6, min_06=2, min_03=3, min_00=4):
    """
    Plot the status of guide catalogs which meet specific magnitude distribution
    requirements.  Make a mollweide visualization using healpy.
    """
    plt.close(11)

    good_guide_catalogs = np.zeros(len(stars_list), dtype=int)
    for ii, stars in enumerate(stars_list):
        n_06 = np.sum((stars['MAG_ACA'] <= mag_limit - 0.6))
        n_03 = np.sum((stars['MAG_ACA'] <= mag_limit - 0.3))
        n_00 = np.sum((stars['MAG_ACA'] <= mag_limit))
        good = n_06 >= min_06 and n_03 >= min_03 and n_00 >= min_00
        good_guide_catalogs[ii] = good

    healpy.mollview(good_guide_catalogs, cmap=cm.jet_r, fig=11, min=-1, max=1, cbar=False)
    return np.sum(good_guide_catalogs) / len(good_guide_catalogs)


def plot_series_good_guide_catalogs():
    """
    Make a series of images suitable for animation:
    % convert -delay 100 good_guide_catalogs_*.png  -loop 0 good_guide_catalogs_anim.gif
    """
    for ii, mag_limit in enumerate(np.arange(10.6, 9.5, -0.1)):
        outfile = 'good_guide_catalogs_{:02d}_{:.1f}.png'.format(ii, mag_limit)
        frac_good = plot_good_guide_catalogs(mag_limit)
        plt.title('Good catalogs for mag limit {:5.1f} mag ({:6.2f}%)'
                  .format(mag_limit, frac_good * 100.0))
        plt.savefig(outfile)


def plot_n_stars(mag_limit=10.6):
    """
    Plot the number of stars brighter than ``mag_limit`` within ``radius`` degrees.
    Makes a mollweide visualization using healpy.
    """
    plt.close(10)

    n_map = [np.sum((stars['MAG_ACA'] <= mag_limit)) for stars in stars_list]
    n_map = np.array(n_map)

    healpy.mollview(n_map, cmap=cm.jet_r, fig=10, min=0, max=8)


def plot_series_n_stars():
    """
    Make a series of images suitable for animation:
    % convert -delay 50 n_stars_*.png  -loop 0 n_stars_anim.gif
    """
    for ii, mag_limit in enumerate(np.arange(10.6, 9.5, -0.1)):
        outfile = 'n_stars_{:02d}_{:.1f}.png'.format(ii, mag_limit)
        plot_n_stars(mag_limit)
        plt.title('Magnitude limit {:.1f} mag'.format(mag_limit))
        plt.savefig(outfile)


def make_ra_dec_stars_list():
    """
    Make a list of stars structured arrays.  Each list element has::

      (ra, dec, stars)

    ``stars`` is a structured array from the microagasc that has the brightest 30
    candidate guide/acq stars within a 1.0 degree radius of ``ra``, ``dec``.
    """
    npix = healpy.nside2npix(NSIDE)
    ra_dec_stars_list = []
    for ipix in range(0, npix):
        theta, phi = healpy.pix2ang(NSIDE, ipix)
        dec = 90 - np.degrees(theta)
        ra = np.degrees(phi)
        print ra, dec

        stars = agasc.get_agasc_cone(ra, dec, 1.0, agasc_file='microagasc.h5')
        ok = ((stars['CLASS'] == 0)
              & (stars['MAG_ACA'] < 10.6)
              & (stars['MAG_ACA'] > 5.8)
              & (stars['ASPQ1'] == 0)
              & (stars['COLOR1'] != 0.700))
        stars_ok = stars[ok]
        stars_ok = np.sort(stars_ok, order=['MAG_ACA'])
        ra_dec_stars_list.append((ra, dec, stars_ok[:30]))

    pickle.dump(ra_dec_stars_list, open('ra_dec_stars_list.pkl', 'w'), protocol=-1)
