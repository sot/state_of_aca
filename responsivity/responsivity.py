"""
Plot the delta magnitude histories of the 4 stars which have been tracked more than N
(300) times.

>>> run -i responsivity
>>> repeats = get_repeats()
>>> plot_repeat_star_mags(repeats)
"""

import numpy as np
import matplotlib.pyplot as plt

from Ska.DBI import DBI
from Chandra.Time import DateTime


cache = {}


def plot_repeat_star_mags(repeats):
    plt.close(1)
    plt.figure(1, figsize=(6, 4))

    aca_db = DBI(server='sybase', dbi='sybase', user='aca_read')
    for agasc_id in repeats['id']:
        print agasc_id
        if agasc_id in cache:
            obsdata = cache[agasc_id]
        else:
            obsdata = aca_db.fetchall("select * from trak_stats_data where id = {} "
                                      "order by kalman_tstart".format(agasc_id))
            cache[agasc_id] = obsdata
        years = DateTime(obsdata['kalman_tstart']).frac_year
        scatter = np.random.uniform(-0.5, 0.5, size=len(years))
        dmags = obsdata['aoacmag_mean'] - np.median(obsdata['aoacmag_mean'])
        plt.plot(years + scatter, dmags, '.', label='ID {}'.format(agasc_id))
    aca_db.conn.close()

    plt.xlabel('Year')
    plt.ylabel('Delta Mag')
    plt.grid()
    plt.ylim(-0.1, 0.1)
    plt.legend(loc='upper left', fontsize=10)
    plt.title('ACA Responsivity')
    plt.tight_layout()
    plt.savefig('responsivity.png')


def get_repeats(n_repeats=300):  # default gives us 4 stars
    aca_db = DBI(server='sybase', dbi='sybase', user='aca_read')
    repeats = aca_db.fetchall("select id, count(id) as num_obs from trak_stats_data "
                              "group by id having (count(id) > {})".format(n_repeats))
    repeats = repeats[repeats['id'] >= 20]
    aca_db.conn.close()

    return repeats
