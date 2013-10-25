"""
Plot the guide star limiting magnitude vs time.  This defines limiting magnitude as
the brightness for which the ensemble no_trak fraction is 50%.  No_trak fraction is
defined as the fraction of guide stars that are not tracked for over 5% of the
observation.

TO DO: this should be renamed plot_guide_mag_limit (along with the function) for
consistency with plot_acq_mag_limit.
"""

import sys
sys.path.insert(0, '../dark_current_model')

from itertools import izip

import numpy as np
import matplotlib.pyplot as plt
from Ska.Matplotlib import plot_cxctime

import dark_models


def plot_limiting_mag():
    """
    Plot the future limiting magnitude for three cases of max T_ccd (-14, -11, -7).
    """
    warm_threshold = 100.

    plt.close(1)
    plt.figure(1, figsize=(6, 4))
    color = 'g'

    for t_ccd_max, offset, color in izip((-14, -11, -7), (0, 0.0025, 0.005), ('r', 'g', 'b')):
        dates, t_ccds = dark_models.ccd_temperature_model(t_ccd_max, start='2014:001',
                                                          stop='2018:001', n=20)
        warm_fracs = []
        for time, t_ccd in izip(dates.secs, t_ccds):
            warm_frac = dark_models.get_warm_fracs(warm_threshold, date=time, T_ccd=t_ccd)
            warm_fracs.append(warm_frac + offset)
        n100 = np.array(warm_fracs)
        mag_limit = np.log10(0.5 / n100) / 1.18 + 10.09
        plot_cxctime(dates.secs, mag_limit, color + '-',
                     label='T_ccd < {} C'.format(t_ccd_max), linewidth=1.5)

    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.03
    plt.xlim(x0 - dx, x1 + dx)
    plt.grid()
    plt.legend(loc='best', fontsize=12)
    plt.title('ACA Limiting magnitude vs. time')
    plt.ylabel('ACA Mag')
    plt.tight_layout()
    plt.savefig('limiting_mag.png')
