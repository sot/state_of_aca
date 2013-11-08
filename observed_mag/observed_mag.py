"""
Plot observed magnitude offsets for the last 4 years.

>>> run -i observed_mag.py
>>> plot_mag_diff_hist()
"""

import matplotlib.pyplot as plt
import numpy as np

from Chandra.Time import DateTime


if 'gui' not in globals() or 'acq' not in globals():
    gui = np.load('/proj/sot/ska/data/gui_stat_reports/guide_stats_with_temp_and_warmpix.npy')
    acq = np.load('/proj/sot/ska/data/acq_stat_reports/acq_stats_with_temp_and_warmpix.npy')


def plot_mag_diff_hist():
    plt.close('all')

    years = np.arange(-3, 1) + int(DateTime().frac_year)
    for year in years:
        tstart = DateTime(year, format='frac_year').secs
        tstop = DateTime(year + 1, format='frac_year').secs

        range_acq = acq[(acq['mag_obs'] != 13.9375)
                        & (acq['obc_id'] == 'ID')
                        & (acq['tstart'] > tstart)
                        & (acq['tstart'] < tstop)]
        range_gui = gui[(gui['type'] != 'FID')
                        & (gui['kalman_tstart'] > tstart)
                        & (gui['kalman_tstart'] < tstop)]

        acq_deltas = range_acq['mag_obs'] - range_acq['mag']
        gui_deltas = range_gui['aoacmag_mean'] - range_gui['mag_exp']
        deltas = np.hstack([acq_deltas, gui_deltas])

        std = np.std(deltas)
        mean = np.mean(deltas)
        ok = np.abs(deltas - mean) < std
        sigma_clip_mean = np.mean(deltas[ok])
        print 'Mean within 1-sigma for year {} is {:.4f}'.format(year, sigma_clip_mean)

        plt.figure(figsize=(4, 3))
        plt.hist(deltas, bins=np.arange(-2.3, 2.3, .04), log=True, histtype='stepfilled', alpha=0.5)
        plt.xlim(-2.3, 2.3)
        plt.ylabel('N Stars')
        plt.xlabel('Observed Mag - Catalog Mag')
        plt.title('Observed Magnitude Offsets {}'.format(year), fontsize='medium')
        plt.grid()
        plt.text(0.3, 4000., 'Mean = {:6.4f} mag'.format(sigma_clip_mean),
                 bbox=dict(facecolor='red', alpha=0.3), fontsize=8)
        plt.tight_layout()
        plt.savefig('mag_diff_hist_{}.png'.format(year))
