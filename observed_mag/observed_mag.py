import matplotlib.pyplot as plt
import numpy as np
from Chandra.Time import DateTime

def plot_mag_diff_hist():
    gui = np.load('/proj/sot/ska/data/gui_stat_reports/guide_stats_with_temp_and_warmpix.npy')
    acq = np.load('/proj/sot/ska/data/acq_stat_reports/acq_stats_with_temp_and_warmpix.npy')

    for year in [2010, 2011, 2012, 2013]:
        tstart = DateTime('{}:001'.format(year)).day_start().secs
        tstop = DateTime('{}:001'.format(year + 1)).day_start().secs
        range_acq = acq[(acq['mag_obs'] != 13.9375)
                  & (acq['obc_id'] == 'ID')
                  & (acq['tstart'] > tstart)
                  & (acq['tstart'] < tstop)]
        range_gui = gui[(gui['type'] != 'FID')
                  & (gui['kalman_tstart'] > tstart)
                  & (gui['kalman_tstart'] < tstop)]

        acq_deltas = range_acq['mag'] - range_acq['mag_obs']
        gui_deltas = range_gui['mag_exp'] - range_gui['aoacmag_mean']
        deltas = np.hstack([acq_deltas, gui_deltas])
        plt.figure(figsize=(6, 4))
        plt.hist(deltas, bins=np.arange(-2.3, 2.3, .025), log=True)
        plt.xlim(-2.3, 2.3)
        plt.ylabel('N Stars')
        plt.xlabel('Catalog Mag - Observed Mag')
        plt.title('Acq and Guide Catalog vs Observed Mag\n Year {}'.format(year))
        plt.tight_layout()
        plt.savefig('mag_diff_hist_{}.png'.format(year))
