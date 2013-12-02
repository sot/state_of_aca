
import numpy as np
import matplotlib.pyplot as plt
from  Ska.Matplotlib import cxctime2plotdate
dac_vals = np.array([])
dac_times = np.array([])
from Chandra.Time import DateTime
for year in range(2002, 2014):
    print year
    y_dac_vals = np.load('{}_vals.npy'.format(year))
    y_dac_mask = np.load('{}_mask.npy'.format(year))
    dac_vals = np.append(dac_vals, y_dac_vals[~y_dac_mask])
    del y_dac_vals
    y_dac_times = cxctime2plotdate(np.load('{}_times.npy'.format(year)))
    dac_times = np.append(dac_times, y_dac_times[~y_dac_mask])
    del y_dac_times
    del y_dac_mask

plt.figure(figsize=(6,4))
plt.plot_date(dac_times, dac_vals, 'b.', markersize=1)
plt.ylim(230, 530)
plt.xlabel('Year')
plt.ylabel('DAC level')
plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=30)
plt.tight_layout()
plt.grid()
plt.savefig('dac_level_over_time.png')
