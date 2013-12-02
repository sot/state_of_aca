#!/usr/bin/env python

import os
import numpy as np

import Ska.DBI
from Chandra.Time import DateTime
import Ska.Numpy
from Ska.engarchive import fetch_sci

import warm_pix


def get_ccd_temps(acqs):
    """
    Get AACCCDPT temps for each acquisition
    """
    tempccd = fetch_sci.Msid('AACCCDPT', '2000:001:00:00:00.000', DateTime().date, stat='5min')
    return tempccd.vals[np.searchsorted(tempccd.times,
                                        acqs['tstart'])]

def get_warm_pix(limit, temps, frac_years):
    """
    Get n100 warm pixel values from table lookup.
    """
    # I did not figure out vectorization over the table
    return [warm_pix.est_warm_pix(limit, temp, frac_year) for temp, frac_year in zip(temps, frac_years)]



sqlaca = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read', database='aca', numpy=True)
min_acq_time = DateTime('2000:001:00:00:00.000')

cols_except_date = ('obsid', 'obi', 'tstart', 'tstop', 'slot', 'idx', 'cat_pos',
                    'type', 'agasc_id', 'obc_id', 'yang', 'zang', 'mag', 'color',
                    'halfw', 'mag_obs', 'yang_obs', 'zang_obs', 'y_offset', 'z_offset',
                    'd_mag', 'd_yang', 'd_zang', 'revision')

all_acq = sqlaca.fetchall('select {} from acq_stats_data where tstart >= {} order by tstart'.format(
                ",".join(cols_except_date),
                min_acq_time.secs))

# write out the warm pixel data to a np save file every time acq stats are updated
temps = get_ccd_temps(all_acq)
acq_with_temp = Ska.Numpy.add_column(all_acq, 'tempccd', np.array(temps))
frac_years = DateTime(all_acq['tstart']).frac_year
acq_with_frac_year = Ska.Numpy.add_column(acq_with_temp, 'frac_year', frac_years)
acq_with_warm_pix = acq_with_frac_year
for limit in [50, 75, 100, 125, 150, 200, 1000]:
    lim_warm_pix = get_warm_pix(limit, temps, frac_years)
    acq_with_warm_pix = Ska.Numpy.add_column(acq_with_warm_pix, 'n{}'.format(str(limit)), lim_warm_pix)
np.save('all_acq.npy', acq_with_warm_pix)
