#!/usr/bin/env python
"""
NOTE: this script is mostly deprecated in favor of using the fit_sota_model
IPython notebook (which does a simultaneous fit of unbinned data) in the
aca_stats project.
"""

from __future__ import division

import os
import sys
import numpy as np
import logging
from glob import glob
import json
import jinja2

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
        matplotlib.use('Agg')
import matplotlib.pyplot as plt


import mx.DateTime
from Chandra.Time import DateTime
from Ska.report_ranges import in_range, timerange, get_next

import time
nowdate = time.ctime()
print "---------- acq stat reports summary update at %s ----------" % (nowdate)

task = 'acq_stat_reports'
TASK_SHARE = os.path.join(os.environ['SKA'],'share', task)
TASK_DATA =  os.path.join(os.environ['SKA'],'data', task)
#TASK_SHARE = '.'


from star_error import high_low_rate

datadir = os.path.join(os.environ['SKA'], 'data', task)
if not os.path.exists(datadir):
    os.makedirs(datadir)
#datadir = os.path.join('/proj/sot/ska', 'data' , task )
#plotdir = os.path.join(os.environ['SKA'], 'www', 'ASPECT', task, 'summary')
plotdir = "./acq_summary_plots"
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

time_pad = .05

#data = { 'month': glob(os.path.join(datadir, '????', 'M??', 'rep.json')),
#         'quarter': glob(os.path.join(datadir, '????', 'Q?', 'rep.json')),
#         'semi': glob(os.path.join(datadir, '????', 'S?', 'rep.json')),
#         'year': glob(os.path.join(datadir, '????', 'YEAR', 'rep.json')),}

figmap = { 'fail_rate' : 1}

all_acq = np.load('all_acq.npy')
#all_acq = np.load(os.path.join(TASK_DATA, 'all_acq.npy'))
stars = all_acq


mag_ranges = {'08.00': {'bright': 8,
                      'faint': 8.5},
              '08.50': {'bright': 8.5,
                      'faint': 9},
              '09.00': {'bright': 9,
                      'faint': 9.25},
              '09.25': {'bright': 9.25,
                      'faint': 9.5},
              '09.50': {'bright': 9.5,
                      'faint': 9.75},
              '09.75': {'bright': 9.75,
                      'faint': 10},
              '10.00': {'bright': 10,
                      'faint': 10.3},
              '10.30': {'bright': 10.3,
                      'faint': 10.6},
              '10.60': {'bright': 10.6,
                      'faint': 11}}

mag_list = sorted(mag_ranges)


data = {}
for range_type in ['month', 'quarter', 'semi', 'year']:
    t0 = timerange(in_range(range_type, DateTime('1998:001').mxDateTime))
    now = mx.DateTime.now()
    data[range_type] = {}
    for mag in mag_ranges:
        t = t0.copy()
        data[range_type][mag] = []
        while t['stop'] < now:
            new_t = get_next(t)
            t = new_t
            range_acqs = all_acq[(all_acq['tstart'] >= DateTime(new_t['start']).secs)
                                 & (all_acq['tstop'] < DateTime(new_t['stop']).secs)
                                 & (all_acq['mag'] < mag_ranges[mag]['faint'])
                                 & (all_acq['mag'] >= mag_ranges[mag]['bright'])]
            good = range_acqs[range_acqs['obc_id'] == 'ID']
            bad = range_acqs[range_acqs['obc_id'] == 'NOID']
            n50_mean = np.mean(range_acqs['n50'])
            n75_mean = np.mean(range_acqs['n75'])
            n100_mean = np.mean(range_acqs['n100'])
            n125_mean = np.mean(range_acqs['n125'])
            n150_mean = np.mean(range_acqs['n150'])
            n200_mean = np.mean(range_acqs['n200'])
            n1000_mean = np.mean(range_acqs['n1000'])
            if len(range_acqs):
                err_high, err_low = high_low_rate(len(bad),
                                                  len(range_acqs))
                data[range_type][mag].append([((DateTime(new_t['start']).secs
                                               + DateTime(new_t['stop']).secs)/2),
                                              (DateTime((DateTime(new_t['start']).secs
                                               + DateTime(new_t['stop']).secs)/2).frac_year),
                                              new_t['year'], new_t['subid'],
                                              (len(bad) / len(range_acqs)),
                                              err_high,
                                              err_low,
                                              (n50_mean / (1024*1024)),
                                              (n75_mean / (1024*1024)),
                                              (n100_mean / (1024*1024)),
                                              (n125_mean / (1024*1024)),
                                              (n150_mean / (1024*1024)),
                                              (n200_mean / (1024*1024)),
                                              (n1000_mean / (1024*1024))])

        # just overwrite with a recarray
        data[range_type][mag] = np.rec.fromrecords(data[range_type][mag],names=[
                        'tstart', 'frac_year', 'year', 'range_id', 'bad_frac',
                        'err_high', 'err_low', 'n50', 'n75', 'n100', 'n125', 'n150', 'n200', 'n1000'])





import sherpa.ui as ui



ftype = 'fail_rate'
fit_info = {}
fit_start = 2007
for range_type in ['quarter', 'semi']:
    fit_info[range_type] = {}
    for mag in mag_ranges:
        fit_info[range_type][mag] = {}
        ok = data[range_type][mag]['frac_year'] > fit_start
        times = data[range_type][mag][ok]['frac_year']
        bad_frac = data[range_type][mag][ok]['bad_frac']
        err_high = data[range_type][mag][ok]['err_high']
        err_high[err_high == 0] = 0.001
        err_low = data[range_type][mag][ok]['err_low']
        err_low[err_low == 0] = 0.001
        #def my_chi2(data, model, staterror=None, syserror=None, weight=None):
        #    fvec = (data - model)
        #    fvec[fvec > 0] = fvec[fvec > 0] / (staterror[fvec > 0] * err_high[fvec > 0])
        #    fvec[fvec < 0] = fvec[fvec < 0] / (staterror[fvec < 0] * err_low[fvec < 0])
        #    stat = (fvec**2).sum()
        #    return (stat, fvec)
        #def my_err(data):
        #    return np.zeros_like(data)
        for limit in [50, 75, 100, 125, 150, 200, 1000]:
            warm_frac = data[range_type][mag][ok]['n{}'.format(limit)]
            print "range_type {}".format(range_type)
            print "mag {}".format(mag)
            print "limit is {}".format(limit)
            extent = np.max(warm_frac) - np.min(warm_frac)
            wp_min = np.min(warm_frac)
            warm_frac = warm_frac - wp_min
            def scaled_warm_frac(pars, x):
                scaled = pars[1] + warm_frac * pars[0]
                return scaled
            data_id = 1
            ui.set_method('simplex')
            ui.set_stat('chi2datavar')
            #ui.set_stat('leastsq')
            #ui.load_user_stat("chi2custom", my_chi2, my_err)
            #ui.set_stat(chi2custom)
            ui.load_user_model(scaled_warm_frac, 'model')
            ui.add_user_pars('model', ['scale', 'offset'])
            ui.set_model(data_id, 'model')
            ui.load_arrays(data_id,
                           np.array(times),
                           np.array(bad_frac))
            fmod = ui.get_model_component('model')
            fmod.scale.min = 1e-9
            max_err = np.max([data[range_type][mag][ok]['err_high'],
                              data[range_type][mag][ok]['err_low']], axis=0)
            ui.set_staterror(data_id, max_err)
            ui.fit(data_id)
            f = ui.get_fit_results()
            scale = f.rstat ** .5
            ui.set_staterror(data_id, max_err * scale)
            ui.fit()
            f = ui.get_fit_results()
            if f.rstat > 3:
                raise ValueError
            ui.confidence()
            conf = ui.get_confidence_results()
            fit_info[range_type][mag][limit] = dict(fit=str(f),
                                                    conf=str(conf),
                                                    fmod=fmod,
                                                    fit_orig=f,
                                                    conf_orig=conf)

            fig = plt.figure(figsize=(5,3))
            axa = plt.subplot(1, 1, 1)
            axa.plot(data[range_type][mag]['frac_year'],
                     100 * data[range_type][mag]['bad_frac'], 'b.')

            warm_frac_all = data[range_type][mag]['n{}'.format(limit)] - wp_min
            # plot the warm fraction in fail scale to help the autoscale
            axa.plot(data[range_type][mag]['frac_year'],
                     100 * ((warm_frac_all * fmod.scale.val)
                            + fmod.offset.val) , 'b', alpha=.5)

            plt.ylabel("Acq Fail Rate (%)", color='blue')


            axb = axa.twinx()
            axb.plot(data[range_type][mag]['frac_year'],
                     100 * warm_frac_all, 'r')
            plt.ylabel("N{} Warm Pixel CCD Frac (%)".format(limit), color='red')
            plt.tight_layout()
            plt.title("%s fail rate (%s)" % (range_type, mag))
            axa_ylim = axa.get_ylim()
            axb.set_autoscale_on(False)
            axb.set_ylim((axa_ylim[0] - (100 * fmod.offset.val)) / fmod.scale.val,
                         (axa_ylim[1] - (100 * fmod.offset.val)) / fmod.scale.val)

            plt.savefig(os.path.join(plotdir, "summary_%s_%s_%s_%s_wp.png" % (
                                    range_type, ftype, mag, limit)))


data = []
for mag in mag_list:
    range_stars = stars[(stars['tstart'] >= DateTime('{}:000:00:00:00.000'.format(fit_start)).secs)
                        & (stars['mag'] < mag_ranges[mag]['faint'])
                        & (stars['mag'] >= mag_ranges[mag]['bright'])]
    data_row = [mag_ranges[mag]['bright'],
                mag_ranges[mag]['faint'],
                float(np.mean(range_stars['mag'])),
                fit_info['semi'][mag][100]['fmod'].scale.val,
                fit_info['semi'][mag][100]['conf_orig'].parmins[0] or 0,
                fit_info['semi'][mag][100]['conf_orig'].parmaxes[0] or 0,
                fit_info['semi'][mag][100]['fmod'].offset.val,
                fit_info['semi'][mag][100]['conf_orig'].parmins[1] or 0,
                fit_info['semi'][mag][100]['conf_orig'].parmaxes[1] or 0]
    data.append(data_row)
cols = ['mag0', 'mag1', 'mag_mean',
        'scale', 'scale_parmin', 'scale_parmax',
        'offset', 'offset_parmin', 'offset_parmax']
data_rec = np.rec.fromrecords(data, names=cols)
np.save('acq_warm_fit', data_rec)
f = open('acq_warm_fit.dat', 'w')
f.write(",".join(cols))
f.write("\n")
for d in data_rec:
    f.write(",".join("{}".format(n) for n in d.tolist()))
    f.write("\n")
f.close()


ftype = 'fail_rate'
mags = [(mag_ranges[mag]['bright'] + mag_ranges[mag]['faint'])/2 for mag in mag_list]
scales = [fit_info['semi'][mag][100]['fmod'].scale.val for mag in mag_list]
scale_parmax = [fit_info['semi'][mag][100]['conf_orig'].parmaxes[0] for mag in mag_list]
scale_parmin = [fit_info['semi'][mag][100]['conf_orig'].parmins[0] for mag in mag_list]
offsets = [fit_info['semi'][mag][100]['fmod'].offset.val for mag in mag_list]
offset_parmax = [fit_info['semi'][mag][100]['conf_orig'].parmaxes[1] for mag in mag_list]
offset_parmin = [fit_info['semi'][mag][100]['conf_orig'].parmins[1] for mag in mag_list]
fig=plt.figure(figsize=(5, 3))
scale_parmin[0] = 0
plt.errorbar(mags, scales, [scale_parmax, np.array(scale_parmin)*-1.0], linestyle='')
plt.xlabel('Catalog magnitude')
plt.ylabel('Warm Pix Scale')
plt.tight_layout()
plt.grid()
plt.savefig('acq_fail_warm_scale.png')
fig=plt.figure(figsize=(5, 3))
plt.errorbar(mags, offsets, [offset_parmax, np.array(offset_parmin)*-1.0], linestyle='')
plt.xlabel('Catalog magnitude')
plt.ylabel('Warm Pix Offset')
plt.tight_layout()
plt.grid()
plt.savefig('acq_fail_warm_offset.png')
print "---------- acq stat reports summary update complete ----------"
