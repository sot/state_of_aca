import numpy as np
import pickle

# Read in the warm pix lookup tables that live in the obc_bad_stat_warm_pix area
# and return warm pix estimates from the est_warm_pix method

warm_pix_dir = './pix_lookup_data/'
warm_pix_info = pickle.load(open('{}/warm_pix_00050.pkl'.format(warm_pix_dir)))
# rows are per cal date in the lookup table
lookup_table = {}
for limit_thresh in [50, 75, 100, 125, 150, 200, 1000]:
    lookup_table[limit_thresh] = np.load(open('{}/lookup_table_{:05}.np'.format(
                warm_pix_dir, limit_thresh)))

# this bookkeeping got convoluted...
bins = np.array([float(i) for i in lookup_table[100].dtype.names])

def est_warm_pix(limit, temp, frac_year):
    """
    find the closest temperature "bin" and then interpolate between the
    dark currents based on time to get an estimated fraction for this temperature
    and time

    :param temp: temperature in C
    :param frac_year: frac year like DateTime().frac_year
    :returns: temp in C
    """
    temp_bin_idx = np.argmin(np.abs(temp - bins))
    bin_name = lookup_table[limit].dtype.names[temp_bin_idx]
    temperature_time_slice = lookup_table[limit][bin_name]
    return np.interp(frac_year,
                     warm_pix_info['frac_year'].astype(float),
                     temperature_time_slice)
