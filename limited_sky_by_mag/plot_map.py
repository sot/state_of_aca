import os
import matplotlib.pyplot as plt
import agasc
import healpy
import numpy as np

NSIDE = 128
limits = [10.13, 10.02, 9.92, 10.38, 10.29, 10.2, 10.0, 10.1, 10.3]

def plot_map(n, limit, radius='small'):
    """
    Plot the sky area there there are fewer than N stars brighter than limit
    in a radius of either 'small' (0.7 deg) or 'large' (0.9 deg).  For speed, limit
    should be one of the mag limits already evaluated
    limits = [10.13, 10.02, 9.92, 10.38, 10.29, 10.2, 10.0, 10.1, 10.3]
    Makes a mollweide visualization using healpy.
    """
    if limit not in limits:
        raise ValueError("not calculated for limit {}".format(limit))
    if radius == 'small':
        if not os.path.exists('small_{}.npy'.format(limit)):
            make_maps()
        if radius == 'small':
            field = np.load('small_{}.npy'.format(limit))
        else:
            field = np.load('big_{}.npy'.format(limit))
        if len(field) != healpy.nside2npix(NSIDE):
            make_maps()
            if radius == 'small':
                field = np.load('small_{}.npy'.format(limit))
            else:
                field = np.load('big_{}.npy'.format(limit))
        map = (field > n).astype('int')
        healpy.mollview(map, xsize=2000)

def make_maps():
    npix = healpy.nside2npix(NSIDE)
    small_ns = {limit: [] for limit in limits}
    big_ns = {limit: [] for limit in limits}
    ras = []
    decs = []
    for ipix in range(0, npix):
        theta, phi = healpy.pix2ang(NSIDE, ipix)
        dec = 90 - np.degrees(theta)
        ra = np.degrees(phi)
        ras.append(ra)
        decs.append(dec)
        print ra, dec
        star_small = agasc.get_agasc_cone(ra, dec, 0.7)
        star_big = agasc.get_agasc_cone(ra, dec, 0.9)
        for limit in limits:
            stars = star_small
            ok = ((stars['CLASS'] == 0)
                  & (stars['MAG_ACA'] < limit)
                  & (stars['MAG_ACA'] > 5.8)
                  & (stars['ASPQ1'] == 0)
                  & (stars['COLOR1'] != 0.700))
            n = len(stars[ok])
            small_ns[limit].append(n)
            stars = star_big
            ok = ((stars['CLASS'] == 0)
                  & (stars['MAG_ACA'] < limit)
                  & (stars['MAG_ACA'] > 5.8)
                  & (stars['ASPQ1'] == 0)
                  & (stars['COLOR1'] != 0.700))
            n = len(stars[ok])
            big_ns[limit].append(n)
    np.save('ra', np.array(ras))
    np.save('dec', np.array(decs))
    for limit in limits:
        np.save('small_{}.npy'.format(limit), np.array(small_ns[limit]))
        np.save('big_{}.npy'.format(limit), np.array(big_ns[limit]))

