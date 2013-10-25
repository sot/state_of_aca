"""
Plot the acquisition star limiting magnitude vs time.  Limiting magnitude is defined as
the brightness for which the probability of success ID during acquisition is 50%. 
"""

import sys
sys.path.insert(0, '../dark_current_model')
from itertools import izip

import matplotlib.pyplot as plt
from Ska.Matplotlib import plot_cxctime
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


def plot_acq_mag_limit():
    """
    Plot the future acq star limiting magnitude for three cases of
    max T_ccd (-14, -11, -7).
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

        mag_limits = []
        for warm_frac in warm_fracs:
            mag_limit = optimize.brentq(acq_success_prob, 6.0, 12.0, args=(warm_frac, 0.5))
            mag_limits.append(mag_limit)

        plot_cxctime(dates.secs, mag_limits, color + '-',
                     label='T_ccd < {} C'.format(t_ccd_max), linewidth=1.5)

    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.03
    plt.xlim(x0 - dx, x1 + dx)
    plt.grid()
    plt.legend(loc='best', fontsize=12)
    plt.title('Acquisition magnitude limit vs. time')
    plt.ylabel('ACA Mag')
    plt.tight_layout()
    plt.savefig('acq_mag_limit.png')
