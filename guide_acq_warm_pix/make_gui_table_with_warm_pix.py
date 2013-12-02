#!/usr/bin/env python

# Acquisition Statistics Report generation

import os
import sys
import numpy as np
import logging

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

import jinja2
import scipy.stats

import Ska.DBI
import Ska.Numpy
from Chandra.Time import DateTime
import Ska.Matplotlib
import Ska.report_ranges
from Ska.engarchive import fetch_sci
from star_error import high_low_rate
import warm_pix

task = 'gui_stat_reports'
TASK_SHARE = os.path.join(os.environ['SKA'],'share', task)
TASK_DATA = os.path.join(os.environ['SKA'], 'data', task)
#TASK_SHARE = "."

jinja_env = jinja2.Environment(
	loader=jinja2.FileSystemLoader(os.path.join(TASK_SHARE, 'templates')))

logger = logging.getLogger(task)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(message)s')


def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--webdir",
                      default="/proj/sot/ska/www/ASPECT/gui_stat_reports",
                      help="Output web directory")
    parser.add_option("--datadir",
                      default="/proj/sot/ska/data/gui_stat_reports",
                      help="Output data directory")
    parser.add_option("--url",
		      default="/mta/ASPECT/gui_stat_reports/")
    parser.add_option("--verbose",
                      type='int',
                      default=1,
                      help="Verbosity (0=quiet, 1=normal, 2=debug)")
    parser.add_option("--bad_thresh",
		      type='float',
		      default=0.05)
    parser.add_option("--obc_bad_thresh",
		      type='float',
		      default=0.05)
    parser.add_option("--days_back",
		      default=30,
		      type='int'),
    opt, args = parser.parse_args()
    return opt, args



def get_ccd_temps(trak):
    """                                                                                          
    Get AACCCDPT temps for each tracked object
    """
    tempccd = fetch_sci.Msid('AACCCDPT', '2000:001:00:00:00.000', DateTime().date, stat='5min')
    mid_time = np.mean([trak['kalman_tstart'], trak['kalman_tstop']], axis=0)
    return tempccd.vals[np.searchsorted(tempccd.times,
                                        mid_time)]

def get_warm_pix(limit, temps, frac_years):
    """                                                                                          
    Get n100 warm pixel values from table lookup.                                                
    """
    # I did not figure out over the table                                                        
    return [warm_pix.est_warm_pix(limit, temp, frac_year) for temp, frac_year in zip(temps, frac_years)\
]



def main(opt):
    """
    Update star statistics plots.  Mission averages are computed with all stars
    from 2003:001 to the end of the interval.
    """
    
    sqlaca = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read', database='aca', numpy=True)
    min_time = DateTime('2003:001:00:00:00.000')

    data_table = 'trak_stats_data'

    cols_except_date = ('obsid',
                        'obi',
                        'slot',
                        'idx',
                        'cat_pos',
                        'id',
                        'type',
                        'color',
                        'cyan_exp',
                        'czan_exp',
                        'mag_exp',
                        'kalman_datestart',
                        'kalman_datestop',
                        'kalman_tstart',
                        'kalman_tstop',
                        'aoacyan_min',
                        'aoacyan_mean',
                        'aoacyan_max',
                        'aoacyan_rms',
                        'aoacyan_median',
                        'aoacyan_5th',
                        'aoacyan_95th',
                        'aoaczan_min',
                        'aoaczan_mean',
                        'aoaczan_max',
                        'aoaczan_rms',
                        'aoaczan_median',
                        'aoaczan_5th',
                        'aoaczan_95th',
                        'aoacmag_min',
                        'aoacmag_mean',
                        'aoacmag_max',
                        'aoacmag_rms',
                        'aoacmag_median',
                        'aoacmag_5th',
                        'aoacmag_95th',
                        'n_samples',
                        'not_tracking_samples',
                        'bad_status_samples',
                        'obc_bad_status_samples',
                        'common_col_samples',
                        'sat_pix_samples',
                        'def_pix_samples',
                        'quad_bound_samples',
                        'ion_rad_samples',
                        'mult_star_samples',
                        'sample_interval_secs',
                        'revision')

    all_gui = sqlaca.fetchall('select {} from trak_stats_data where kalman_tstart >= {} order by kalman_tstart'.format(
                    ",".join(cols_except_date),
                    min_time.secs))
    temps = get_ccd_temps(all_gui)
    gui_with_temp = Ska.Numpy.add_column(all_gui, 'tempccd', np.array(temps))
    frac_years = DateTime(all_gui['kalman_tstart']).frac_year
    gui_with_frac_year = Ska.Numpy.add_column(gui_with_temp, 'frac_year', frac_years)
    gui_with_warm_pix = gui_with_frac_year
    limits = [50, 75, 100, 125, 150, 200, 1000]
    for limit in limits:
        warm_pix = get_warm_pix(limit, temps, frac_years)
        gui_with_warm_pix = Ska.Numpy.add_column(gui_with_warm_pix, 'n{}'.format(str(limit)), warm_pix)
    np.save('all_gui.npy', gui_with_warm_pix)



if __name__ == '__main__':
    opt, args = get_options()
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARN)
    if opt.verbose == 2:
	    ch.setLevel(logging.DEBUG)
    if opt.verbose == 0:
	    ch.setLevel(logging.ERROR)
    logger.addHandler(ch) 
    main(opt)
