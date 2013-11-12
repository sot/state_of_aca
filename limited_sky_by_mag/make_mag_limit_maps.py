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


if 'stars_list' not in globals() and os.path.exists('ra_dec_stars_list.pkl'):
    ras, decs, stars_list = pickle.load(open('ras_decs_stars_list.pkl', 'rb'))


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


def plot_n_stars(mag_limit=10.6):
    """
    Plot the number of stars brighter than ``mag_limit`` within ``radius`` degrees.
    Makes a mollweide visualization using healpy.
    """
    plt.close(10)
    max_radius = 0.7836  # same area (1.92 deg**2) as CCD (5000 arcsec on a side).

    n_map = [np.sum((stars['MAG_ACA'] <= mag_limit) & (stars['radius'] < max_radius))
             for stars in stars_list]
    n_map = np.array(n_map)

    healpy.mollview(n_map, cmap=cm.jet_r, fig=10, min=0, max=8)


def plot_series_n_stars():
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
