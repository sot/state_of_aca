import numpy as np
from itertools import izip
import sys
sys.path.insert(0, '../dark_current_model')
from scipy import optimize
import dark_models

target_year = 2018
guide_stars = np.load('/proj/sot/ska/data/gui_stat_reports/guide_stats_with_temp_and_warmpix.npy')
acq_stars = np.load('/proj/sot/ska/data/acq_stat_reports/acq_stats_with_temp_and_warmpix.npy')

warm_threshold = 100


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
    """probability of at least n stars"""
    acq_prob = [0 for s in range(0,9)]
    # we'll just call the at least 0 stars = 1
    acq_prob[0] = 1
    for n_star in range(1, 9):
        from itertools import combinations
        n_prob = 0
        for comb in combinations(range(0, 8), n_star):
            comb_prob = 1.0
            for slot in comb:
                comb_prob = comb_prob * p[slot]
            n_prob = n_prob + comb_prob - (n_prob * comb_prob)
        acq_prob[n_star] = n_prob
    return acq_prob

#    raise ValueError
#
#    for i in range(0, 2**len(p)):
#        print i
#        total_prob = 1.0
#        n_acq = 0
#        prob = 0
#        for slot in range(0, len(p)):
#            if ((i & 2**slot) > 1):
#                prob = p[slot]
#                n_acq = n_acq + 1
#                print "slot {} hit {:.4f}".format(slot, prob),
#            else:
#                print "slot {} fail {:.4f}".format(slot, prob),
#                prob = 1 - p[slot]
#            total_prob = total_prob * prob
#        print "\n"
#        acq_prob[n_acq] = acq_prob[n_acq] + total_prob
#    raise ValueError
#    return acq_prob

def cum_prob(acq_prob):
    """probability of n or fewer stars"""
    return [1 - acq_prob[i + 1] for i in range(1, 7)]


#def cum_prob(acq_prob):
#    cum_prob = [0 for s in acq_prob]
#    exp = 0
#    for i in range(1, len(acq_prob)+1):
#        exp = exp + i * acq_prob[i]
#        for j in range(i, len(acq_prob) + 1):
#            cum_prob[i] = cum_prob[i] + acq_prob[j]
#    cum_prob_l = [log(1.0 - s)/log(10.0) for s in cum_prob]
#    raise ValueError
#    return cum_prob_l


per_obsid_four = { -7: [],
                    -11: [],
                    -14: [] }
per_obsid_three = { -7: [],
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
        mag_limit = optimize.brentq(acq_fail_prob, 6.0, 12.0, args=(warm_frac + offset, 0.5))
        acq_mag_limit[t_ccd_max] = mag_limit

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
        cum = cum_prob(prob_n)
        per_obsid_four[t_ccd_max].append(cum[4])
        per_obsid_three[t_ccd_max].append(cum[3])

plt.hist(per_obsid_three[-7], histtype='bar', color='red', log=True, bins=10**np.linspace(0,1,10))
plt.hist(per_obsid_three[-11], histtype='bar', color='green', log=True, bins=10**np.linspace(0,1,10))
plt.hist(per_obsid_three[-14], histtype='bar', color='blue', log=True, bins=10**np.linspace(0,1,10))
plt.xscale('log')
plt.xlim(0.88728187439428363, 10.854267376291101)
plt.grid()
r = Rectangle((0, 0), 1, 1, fc="r")
g = Rectangle((0, 0), 1, 1, fc="g")
b = Rectangle((0, 0), 1, 1, fc="b")
legend([r, g, b], ['-7', '-11', '-14'], loc="upper right")
plt.xlabel('Prob. of 3 or fewer ACQ stars (10^n)')
plt.ylabel('N observations')
plt.show()
