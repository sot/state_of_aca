import os
import numpy as np
from itertools import izip
import sys
sys.path.insert(0, '../dark_current_model')
from scipy import optimize
import dark_models

def acq_success_prob(mag, warm_frac, prob_offset=0):
    """
    Calculate probability of acquisition success for a star with ``mag``
    magnitude and a CCD warm fraction ``warm_frac``.  Uses the empirical relation:

       P_acq_success = offset(mag) + scale(mag) * warm_frac

    In ../guide_acq_success/plot_acq_success.py we find the best fit relation:

      log10(scale) = 0.185 + 0.990 * (mag - 10) + -0.491 * (mag - 10)**2
      log10(offset) = -1.489 + 0.888 * (mag - 10) + 0.280 * (mag - 10)**2
    """
    mag10 = mag - 10.0
    scale = 10. ** (0.185 + 0.990 * mag10 + -0.491 * mag10 ** 2)
    offset = 10. ** (-1.489 + 0.888 * mag10 + 0.280 * mag10 ** 2)

    return offset + scale * warm_frac - prob_offset


def save_stars_per_catalog():
    target_year = 2018
    guide_stars = np.load('/proj/sot/ska/data/gui_stat_reports/guide_stats_with_temp_and_warmpix.npy')
    acq_stars = np.load('/proj/sot/ska/data/acq_stat_reports/acq_stats_with_temp_and_warmpix.npy')

    warm_threshold = 100

    # guide stars
    gui_mag_limit = {}
    for t_ccd_max, offset, color in izip((-14, -11, -7), (0, 0.0025, 0.005), ('r', 'g', 'b')):

        dates, t_ccds = dark_models.ccd_temperature_model(t_ccd_max,
                                                          start='{}:001'.format(target_year),
                                                          stop='{}:001'.format(target_year), n=1)
        warm_frac = (dark_models.get_warm_fracs(warm_threshold, date=dates.secs[0], T_ccd=t_ccds[0])
                     + offset)
        gui_mag_limit[t_ccd_max] = np.log10(0.5 / warm_frac) / 1.18 + 10.09

    per_obsid_guide = { -7: [],
                       -11: [],
                       -14: [] }
    obsids = np.unique(guide_stars['obsid'])
    for obsid in obsids:
        # ignore cal for the time being
        if obsid > 40000:
            continue
        gs = guide_stars[(guide_stars['obsid'] == obsid) & (guide_stars['type'] != 'FID') & (guide_stars['type'] != 'MON')]
        # ignore multi-obi
        if len(np.unique(gs['obi'])) > 1:
            continue
        for mag in gui_mag_limit:
            n_bright_stars = len(np.flatnonzero(gs['mag_exp'] < gui_mag_limit[mag]))
            per_obsid_guide[mag].append(n_bright_stars)
            print obsid, mag, n_bright_stars


    per_obsid_acq = { -7: [],
                          -11: [],
                          -14: [] }

    acq_mag_limit = {}
    for t_ccd_max, offset, color in izip((-14, -11, -7), (0, 0.0025, 0.005), ('r', 'g', 'b')):
        dates, t_ccds = dark_models.ccd_temperature_model(t_ccd_max,
                                                          start='{}:001'.format(target_year),
                                                          stop='{}:001'.format(target_year), n=1)

        warm_frac = (dark_models.get_warm_fracs(warm_threshold, date=dates.secs[0], T_ccd=t_ccds[0])
                     + offset)

        for time, t_ccd in izip(dates.secs, t_ccds):
            warm_frac = dark_models.get_warm_fracs(warm_threshold, date=time, T_ccd=t_ccd)
            mag_limit = optimize.brentq(acq_success_prob, 6.0, 12.0, args=(warm_frac + offset, 0.5))
            acq_mag_limit[t_ccd_max] = mag_limit

    obsids = np.unique(acq_stars['obsid'])
    for obsid in obsids:
        # ignore cal for the time being
        if obsid > 40000:
            continue
        acq = acq_stars[(acq_stars['obsid'] == obsid)]
        # ignore multi-obi
        if len(np.unique(acq['obi'])) > 1:
            continue
        if len(acq) < 8:
            continue
        for mag in acq_mag_limit:
            n_bright_stars = len(np.flatnonzero(acq['mag'] < acq_mag_limit[mag]))
            per_obsid_acq[mag].append(n_bright_stars)
            print obsid, mag, n_bright_stars

    for mag in [-7, -11, -14]:
        np.save('guide_{}_bright'.format(mag), per_obsid_guide[mag])
        np.save('acq_{}_bright'.format(mag), per_obsid_acq[mag])


def stars_per_obsid():
    mags = [-7, -11, -14]
    for mag in mags:
        if not (os.path.exists('guide_{}_bright.npy'.format(mag))
                and os.path.exists('acq_{}_bright.npy'.format(mag))):
            save_stars_per_catalog()
            break
    guide = dict()
    acq = dict()
    for mag in mags:
        guide[mag] = np.load('guide_{}_bright.npy'.format(mag))
        acq[mag] = np.load('acq_{}_bright.npy'.format(mag))
    return guide, acq
