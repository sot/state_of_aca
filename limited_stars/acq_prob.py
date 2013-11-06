import os
import sys
from itertools import izip
from scipy import optimize
import numpy as np

sys.path.insert(0, '../dark_current_model')
import dark_models

import matplotlib.pyplot as plt


def acq_fail_prob(mag, warm_frac, prob_offset=0):
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


# I was a bit confused by FigureOfMerit in starcheck, so I did this via combinations
def prob_n_stars(p):
    """probability of exactly n stars"""
    acq_prob = []
    for n_star in range(0, len(p) + 1):
        from itertools import combinations
        n_prob = 0
        for comb in combinations(range(len(p)), n_star):
            # find every combination of n_star success
            comb_prob = 1.0
            for slot in range(0, len(p)):
                if slot in comb:
                    # if the slot in the success set, use the probability
                    comb_prob = comb_prob * p[slot]
                    #print " {} ".format(slot),
                else:
                    # if the slot wasn't in the success set, use the fail prob
                    comb_prob = comb_prob * (1 - p[slot])
                    #print "!{} ".format(slot),
            # the chance of n stars is the OR probability of the individual
            # combinations
            n_prob = n_prob + comb_prob
            #print comb_prob
        acq_prob.append(n_prob)
    #for n in range(len(p) + 1):
    #    print "Prob of {} is {}".format(n, acq_prob[n])
    return acq_prob



def obsid_probabilities():
    per_obsid_two = { -7: [],
                       -11: [],
                       -14: [] }
    per_obsid_three = { -7: [],
                        -11: [],
                        -14: [] }
    per_obsid_four = { -7: [],
                        -11: [],
                        -14: [] }

    target_year = 2018
    acq_stars = np.load('/proj/sot/ska/data/acq_stat_reports/acq_stats_with_temp_and_warmpix.npy')
    warm_threshold = 100

    #acq_mag_limit = {}
    for t_ccd_max, offset, color in izip((-14, -11, -7), (0, 0.0025, 0.005), ('r', 'g', 'b')):
        # load these from np save files if they exists
        if (os.path.exists('acq_two_{}.npy'.format(t_ccd_max))
            & os.path.exists('acq_three_{}.npy'.format(t_ccd_max))
            & os.path.exists('acq_four_{}.npy'.format(t_ccd_max))):
                per_obsid_two[t_ccd_max] = np.load('acq_two_{}.npy'.format(t_ccd_max))
                per_obsid_three[t_ccd_max] = np.load('acq_three_{}.npy'.format(t_ccd_max))
                per_obsid_four[t_ccd_max] = np.load('acq_four_{}.npy'.format(t_ccd_max))
                continue

        dates, t_ccds = dark_models.ccd_temperature_model(t_ccd_max,
                                                          start='{}:001'.format(target_year),
                                                          stop='{}:001'.format(target_year), n=1)

        warm_frac = (dark_models.get_warm_fracs(warm_threshold, date=dates.secs[0], T_ccd=t_ccds[0])
                     + offset)

        for time, t_ccd in izip(dates.secs, t_ccds):
            warm_frac = dark_models.get_warm_fracs(warm_threshold, date=time, T_ccd=t_ccd)
            #mag_limit = optimize.brentq(acq_fail_prob, 6.0, 12.0, args=(warm_frac + offset, 0.5))
            #acq_mag_limit[t_ccd_max] = mag_limit

        obsids = np.unique(acq_stars['obsid'])
        for obsid in np.unique(obsids):
            # ignore cal for the time being
            if obsid > 40000:
                continue
            acq = acq_stars[(acq_stars['obsid'] == obsid)]
            # ignore multi-obi
            if len(np.unique(acq['obi'])) > 1:
                continue
            if len(acq) < 8:
                continue
            acq_probs = [1 - acq_fail_prob(star['mag'], warm_frac) for star in acq]
            prob_n = prob_n_stars(acq_probs)
            cum = np.cumsum(prob_n)
            per_obsid_two[t_ccd_max].append(cum[2])
            per_obsid_three[t_ccd_max].append(cum[3])
            per_obsid_four[t_ccd_max].append(cum[4])
        np.save('acq_two_{}'.format(t_ccd_max), per_obsid_two[t_ccd_max])
        np.save('acq_three_{}'.format(t_ccd_max), per_obsid_three[t_ccd_max])
        np.save('acq_four_{}'.format(t_ccd_max), per_obsid_three[t_ccd_max])

    return per_obsid_two, per_obsid_three, per_obsid_four

def plot_n_or_fewer(n):
    if n not in [2, 3, 4]:
        raise ValueError("N must be 2, 3 or 4")
    per_obsid_two, per_obsid_three, per_obsid_four = obsid_probabilities()
    if n == 2:
        per_obsid = per_obsid_two
    if n == 3:
        per_obsid = per_obsid_three
    if n == 4:
        per_obsid = per_obsid_four
    plt.figure(figsize=(6, 4))
    plt.hist(per_obsid_three[-7], histtype='bar',
             color='red', log=True, bins=10**np.linspace(0,1,10))
    plt.hist(per_obsid_three[-11], histtype='bar',
             color='green', log=True, bins=10**np.linspace(0,1,10))
    plt.hist(per_obsid_three[-14], histtype='bar',
             color='blue', log=True, bins=10**np.linspace(0,1,10))
    plt.xscale('log')
    plt.xlim(0.88728187439428363, 10.854267376291101)
    plt.grid()
    r = plt.Rectangle((0, 0), 1, 1, fc="r")
    g = plt.Rectangle((0, 0), 1, 1, fc="g")
    b = plt.Rectangle((0, 0), 1, 1, fc="b")
    plt.legend([r, g, b], ['-7', '-11', '-14'], loc="upper right")
    plt.xlabel('Prob. of {} or fewer ACQ stars (10^n)'.format(n))
    plt.ylabel('N observations')
    plt.tight_layout()
    plt.show()
