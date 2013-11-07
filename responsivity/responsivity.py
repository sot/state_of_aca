from Ska.DBI import DBI
from Chandra.Time import DateTime
import matplotlib.pyplot as plt

# plot the delta magnitude histories of the 4 stars which have been tracked more than N (300) times

def plot_repeat_star_mags():
    aca_db = DBI(server='sybase', dbi='sybase', user='aca_read')
    n_repeats = 300 # gives us 4 stars
    repeats = aca_db.fetchall("select id, count(id) as num_obs from trak_stats_data group by id having (count(id) > {})".format(n_repeats))
    plt.figure(figsize=(6, 3))
    for agasc_id in repeats['id']:
        if agasc_id < 20:
            continue
        print agasc_id
        obsdata = aca_db.fetchall("select * from trak_stats_data where id = {} order by kalman_tstart".format(agasc_id))
        plt.plot(DateTime(obsdata['kalman_tstart']).frac_year, obsdata['mag_exp'] - obsdata['aoacmag_mean'], '.')
    plt.xlabel('Year')
    plt.ylabel('Catalog Mag - Observed Mag')
    plt.grid()
    plt.ylim(-0.17, 0.17)
    plt.tight_layout()
    plt.savefig('responsivity.png')

