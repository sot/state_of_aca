#!/usr/bin/env python
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
import sherpa.ui as ui

import mx.DateTime
from Chandra.Time import DateTime
from Ska.report_ranges import in_range, timerange, get_next

task = 'gui_stat_reports'
TASK_SHARE = os.path.join(os.environ['SKA'],'share', task)
TASK_DATA = os.path.join(os.environ['SKA'], 'data', task)

from star_error import high_low_rate

datadir = os.path.join(os.environ['SKA'], 'data', task)
if not os.path.exists(datadir):
    os.makedirs(datadir)
#datadir = os.path.join('/proj/sot/ska', 'data' , task )
#plotdir = os.path.join(os.environ['SKA'], 'www', 'ASPECT', task, 'summary')
plotdir = "./gui_summary_plots"
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

time_pad = .05

data = { 'month': glob(os.path.join(datadir, '????', 'M??', 'rep.json')),
         'quarter': glob(os.path.join(datadir, '????', 'Q?', 'rep.json')),
         'semi': glob(os.path.join(datadir, '????', 'S?', 'rep.json')),
         'year': glob(os.path.join(datadir, '????', 'YEAR', 'rep.json')),}

figmap = { 'bad_trak' : 1,
           'obc_bad' : 2,
           'no_trak' : 3 }

all_trak = np.load(os.path.join(TASK_DATA, 'all_gui.npy'))
stars = all_trak[all_trak['type'] != 'FID']
del all_trak

mag_ranges = {'08.00': {'bright': 8,
                      'faint': 8.5},
              '08.50': {'bright': 8.5,
                      'faint': 9},
              '09.00': {'bright': 9,
                      'faint': 9.5},
              '09.50': {'bright': 9.5,
                      'faint': 10},
              '10.00': {'bright': 10,
                      'faint': 10.5},
              '10.50': {'bright': 10.5,
                      'faint': 11}}

mag_list = sorted(mag_ranges)

bad_thresh = 0.05
obc_bad_thresh = 0.05

warm_limits = [50, 75, 100, 125, 150, 200]

names = ['tstart', 'frac_year', 'year', 'range_id',
         'mag_mean',
         'mean_bad_trak',
         'bad_trak', 'bad_trak_err_high', 'bad_trak_err_low', 
         'obc_bad', 'obc_bad_err_high', 'obc_bad_err_low', 
         'no_trak', 'no_trak_err_high', 'no_trak_err_low']
for limit in warm_limits:
        names.append('n{}'.format(limit))
        

data = {}
for range_type in ['semi', 'quarter', 'month']:
    t0 = timerange(in_range(range_type, DateTime('1998:001').mxDateTime))
    now = mx.DateTime.now()
    data[range_type] = {}
    for mag in mag_ranges:
        t = t0.copy() 
        data[range_type][mag] = []
        while t['stop'] < now:
            new_t = get_next(t)
            t = new_t
            range_gui = stars[(stars['kalman_tstart'] >= DateTime(new_t['start']).secs)
                              & (stars['kalman_tstop'] < DateTime(new_t['stop']).secs)
                              & (stars['mag_exp'] < mag_ranges[mag]['faint'])
                              & (stars['mag_exp'] >= mag_ranges[mag]['bright'])]
            if len(range_gui) > 1:
                bad_trak_frac = np.mean(range_gui['not_tracking_samples'] / range_gui['n_samples'])
                bad_trak = range_gui[range_gui['not_tracking_samples'] / range_gui['n_samples'] > bad_thresh]
                bad_trak_err_high, bad_trak_err_low = high_low_rate(len(bad_trak), len(range_gui))
                obc_bad = range_gui[range_gui['obc_bad_status_samples'] / range_gui['n_samples'] > obc_bad_thresh]
                obc_bad_err_high, obc_bad_err_low = high_low_rate(len(obc_bad), len(range_gui))
                no_trak = range_gui[range_gui['not_tracking_samples'] == range_gui['n_samples']]
                no_trak_err_high, no_trak_err_low = high_low_rate(len(no_trak), len(range_gui))
                entry = [((DateTime(new_t['start']).secs 
                           + DateTime(new_t['stop']).secs)/2),
                         (DateTime((DateTime(new_t['start']).secs
                                    + DateTime(new_t['stop']).secs)/2).frac_year),
                         new_t['year'], new_t['subid'],
                         np.mean(range_gui['mag_exp']),
                         bad_trak_frac,
                         (len(bad_trak) / len(range_gui)),
                         bad_trak_err_high,
                         bad_trak_err_low,
                         (len(obc_bad) / len(range_gui)),
                         obc_bad_err_high,
                         obc_bad_err_low,
                         (len(no_trak) / len(range_gui)),
                         no_trak_err_high,
                         no_trak_err_low]
                for limit in warm_limits:
                    warm_mean = np.mean(range_gui['n{}'.format(limit)])
                    entry.append(warm_mean / (1024*1024))
                data[range_type][mag].append(entry)
        # just overwrite with a recarray
        data[range_type][mag] = np.rec.fromrecords(data[range_type][mag],names=names)



norm_start = DateTime('2007:001:00:00:00.000').secs
#for ftype in ['bad_trak', 'obc_bad', 'no_trak']:
fit_info = {}
fit_start = 2007
for range_type in ['semi', 'quarter', 'month']:
    fit_info[range_type] = {}
    for mag in mag_ranges:
        if mag not in data[range_type]:
            continue
        fit_info[range_type][mag] = {}

        ok = data[range_type][mag]['frac_year'] > fit_start
        times = data[range_type][mag][ok]['frac_year']
        for ftype in ['bad_trak']:
            fit_info[range_type][mag][ftype] = {}
            bad_frac = data[range_type][mag][ok][ftype]
            err_high = data[range_type][mag][ok]["{}_err_high".format(ftype)]
            err_low = data[range_type][mag][ok]["{}_err_low".format(ftype)]
            err_high[err_high == 0] = .0001
            err_low[err_low == 0] = .0001
            for limit in warm_limits:
                print "range type {}".format(range_type)
                print "mag {}".format(mag)
                print "limit is {}".format(limit)
                print "ftype {}".format(limit)
                warm_frac = data[range_type][mag][ok]['n{}'.format(limit)]
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
                fmod.offset.val = 0
                ui.freeze(fmod.offset)
                max_err = np.max([err_high, err_low], axis=0)
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
                fit_info[range_type][mag][ftype][limit] = dict(fit=str(f),
                                                               conf=str(conf),
                                                               fmod=fmod,
                                                               fit_orig=f,
                                                               conf_orig=conf,
                                                               mag_mean=np.mean(data[range_type][mag][ok]['mag_mean']))
                fig = plt.figure(figsize=(5,3))
                axa = plt.subplot(1, 1, 1)
                axa.plot(data[range_type][mag]['frac_year'],
                         100 * data[range_type][mag][ftype], 'b.')
                # plot the warm fraction in fail scale to help the autoscale
                warm_frac_all = data[range_type][mag]['n{}'.format(limit)] - wp_min
                axa.plot(data[range_type][mag]['frac_year'],
                         100 * ((warm_frac_all * fmod.scale.val)
                                + fmod.offset.val), 'b', alpha=.5)

                plt.ylabel("{} Rate (%)".format(ftype) , color='blue')
                axb = axa.twinx()
                axb.plot(data[range_type][mag]['frac_year'],
                         100 * warm_frac_all, 'r')
                plt.ylabel("N{} Warm Pixel CCD Frac (%)".format(limit), color='red')
                plt.tight_layout()    
                plt.title("%s %s (%s)" % (range_type, ftype, mag))
                axa_ylim = axa.get_ylim()
                axb.set_autoscale_on(False)
                axb.set_ylim((axa_ylim[0] - (100 * fmod.offset.val)) / fmod.scale.val,
                             (axa_ylim[1] - (100 * fmod.offset.val)) / fmod.scale.val)
                plt.savefig(os.path.join(plotdir, "summary_%s_%s_%s_%s_wp.png" 
                                         % (range_type, ftype, mag, limit)))
            #if ftype == 'mean_bad_trak':
            #   continue
            #
            #fig1 = plt.figure(figsize=(5,3))
            #ax1 = fig1.gca()
            #fig2 = plt.figure(figsize=(5,3))
            #ax2 = fig2.gca()
            #
            #ax1.plot(data[range_type][mag]['frac_year'],
            #         data[range_type][mag][ftype],
            #         color = 'black',
            #         linestyle='',
            #         marker='.',
            #         markersize=5)
            #ax1.grid()
            #ax2.errorbar(data[range_type][mag]['frac_year'],
            #             data[range_type][mag][ftype],
            #             yerr = np.array([data[range_type][mag]['{}_err_low'.format(ftype)],
            #                              data[range_type][mag]['{}_err_high'.format(ftype)]]),
            #             color = 'black',
            #             linestyle='',
            #             marker='.',
            #             markersize=5)
            #ax2.grid()
            #fit_file = open(os.path.join(TASK_DATA, "%s_fitfile.json" % ftype), 'r')
        #fit#_file = open(os.path.join(datadir, "%s_fitfile.json" % ftype), 'r')
            #fit_text = fit_file.read()
            #fit = json.loads(fit_text)
            #trend_s_mxd = DateTime(fit['datestart']).mxDateTime
            #trend_start_frac = trend_s_mxd.year + (trend_s_mxd.day_of_year * 1.0 / 365)
            #m = fit['m']
            #b = fit['b']
            #now_mxd = DateTime().mxDateTime
            #now_frac = now_mxd.year + (now_mxd.day_of_year * 1.0 / 365)
            #for ax in [ax1, ax2]:
            #    ax.plot( [trend_start_frac,
            #              now_frac + 1],
            #             [ b,
            #               m * ((now_frac + 1) - trend_start_frac) + b],
            #             'r-')
            #
            #ax2_ylim = ax2.get_ylim()
# pad a bit #below 0 relative to ylim range
            #ax2.set_ylim(ax2_ylim[0] - 0.025*(ax2_ylim[1] - ax2_ylim[0]))
            #ax1.set_ylim(ax2.get_ylim())
            #
            #for ax in [ax1, ax2]:
            #    dxlim = now_frac - 2000
            #    ax.set_xlim(2000,
            #                now_frac + time_pad * dxlim)
        #   # ax = fig.get_axes()[0]
            #    labels = ax.get_xticklabels() + ax.get_yticklabels()
            #    for label in labels:
            #        label.set_size('small')
            #        ax.set_ylabel('Rate', fontsize=12)
            #        ax.set_title("%s %s" % (range_type,ftype), fontsize=12)
            #
            #fig1.subplots_adjust(left=.15)
            #fig2.subplots_adjust(left=.15)
            #fig1.savefig(os.path.join(plotdir, "summary_%s_%s_%s.png" % (range_type, ftype, mag)))
            #fig2.savefig(os.path.join(plotdir, "summary_%s_%s_%s_eb.png" % (range_type, ftype, mag)))
            #plt.close('all')

#for d in data.keys():
#
#
#    data[d].sort()
#    rates =  dict([ (ftype, dict(time=np.array([]),
#                                rate=np.array([]),
#                                err_h=np.array([]),
#                                err_l=np.array([])))
#                    for ftype in ['bad_trak', 'no_trak', 'obc_bad']])
#
#   
#    for p in data[d]:
#        rep_file = open(p, 'r')
#        rep_text = rep_file.read()
#        rep = json.loads(rep_text)
#        for ftype in rates.keys():
#            mxd = DateTime( (DateTime(rep['datestart']).secs
#                             +  DateTime(rep['datestop']).secs) / 2).mxDateTime
#            frac_year = mxd.day_of_year * 1.0 / 365
#            rates[ftype]['time'] = np.append(rates[ftype]['time'],
#                                                mxd.year + frac_year) 
#            for fblock in rep['fail_types']:
#                if fblock['type'] == ftype:
#                    rates[ftype]['rate'] = np.append(rates[ftype]['rate'],
#                                                     fblock['rate'])
#                    rates[ftype]['err_h'] = np.append(rates[ftype]['err_h'],
#                                                      fblock['rate_err_high'])
#                    rates[ftype]['err_l'] = np.append(rates[ftype]['err_l'],
#                                                      fblock['rate_err_low'])
#
#    
#
#    for ftype in ['no_trak', 'bad_trak', 'obc_bad',]:
#
#        
#        fig1 = plt.figure(1,figsize=(5,3))
#        ax1 = fig1.gca()
#        fig2 = plt.figure(2,figsize=(5,3))
#        ax2 = fig2.gca()
#
#        ax1.plot(rates[ftype]['time'],
#                 rates[ftype]['rate'],
#                 color = 'black',
#                 linestyle='',
#                 marker='.',
#                 markersize=5)
#        ax1.grid()
#        ax2.errorbar(rates[ftype]['time'],
#                     rates[ftype]['rate'],
#                     yerr = np.array([rates[ftype]['err_l'],
#                                     rates[ftype]['err_h']]),
#                     color = 'black',
#                     linestyle='',
#                     marker='.',
#                     markersize=5)
#        ax2.grid(
#        fit_file = open(os.path.join(datadir, "%s_fitfile.json" % ftype), 'r')
#        fit_text = fit_file.read()
#        fit = json.loads(fit_text)
#        trend_start_frac = DateTime(fit['datestart']).frac_year
#        m = fit['m']
#        b = fit['b']
#        now_frac = DateTime().frac_year
#        for ax in [ax1, ax2]:
#            ax.plot( [trend_start_frac,
#                      now_frac + 1],
#                     [ b,
#                       m * ((now_frac + 1) - trend_start_frac) + b],
#                     'r-')
#        ax2_ylim = ax2.get_ylim()
#        # pad a bit below 0 relative to ylim range
#        ax2.set_ylim(ax2_ylim[0] - 0.025*(ax2_ylim[1] - ax2_ylim[0]))
#        ax1.set_ylim(ax2.get_ylim())
#        
#        for ax in [ax1, ax2]:
#            dxlim = now_frac - 2000
#            ax.set_xlim(2000,
#                        now_frac + time_pad * dxlim)
#            #    ax = fig.get_axes()[0]
#            labels = ax.get_xticklabels() + ax.get_yticklabels()
#            for label in labels:
#                label.set_size('small')
#            ax.set_ylabel('Rate', fontsize=12)
#            ax.set_title("%s %s" % (d,ftype), fontsize=12)
#
#        fig1.subplots_adjust(left=.15)
#        fig2.subplots_adjust(left=.15)
#        fig1.savefig(os.path.join(plotdir, "summary_%s_%s.png" % (d, ftype)))
#        fig2.savefig(os.path.join(plotdir, "summary_%s_%s_eb.png" % (d, ftype)))
#        plt.close(fig1)
#        plt.close(fig2)                                 
#
#
#
#
#
jinja_env = jinja2.Environment(
	loader=jinja2.FileSystemLoader(os.path.join(TASK_SHARE, 'templates')))
#
outfile = os.path.join(plotdir, 'guide_summary.html')
template = jinja_env.get_template('summary.html')
page = template.render(fit_info=fit_info)
f = open(outfile, 'w')
f.write(page)
f.close()
#                

plt.close('all')
#(mag_ranges[mag]['bright'] + mag_ranges[mag]['faint'])/2,
data = []
for mag in mag_list:
    range_stars = stars[(stars['kalman_tstart'] >= DateTime('{}:000:00:00:00.000'.format(fit_start)).secs)
                        & (stars['mag_exp'] < mag_ranges[mag]['faint'])
                        & (stars['mag_exp'] >= mag_ranges[mag]['bright'])]

    data_row = [mag_ranges[mag]['bright'],
                mag_ranges[mag]['faint'],
                float(np.mean(range_stars['mag_exp'])),
                fit_info['semi'][mag]['bad_trak'][100]['fmod'].scale.val,
                fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmins[0],
                fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmaxes[0]]
                #fit_info['semi'][mag]['bad_trak'][100]['fmod'].offset.val,
                #fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmins[1],
                #fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmaxes[1]]
    data.append(data_row)
cols = ['mag0', 'mag1', 'mag_mean', 'scale', 'scale_parmin', 'scale_parmax']
data_rec = np.rec.fromrecords(data, names=cols)
f = open('gui_warm_fit.dat', 'w')
f.write(",".join(cols))
f.write("\n")
for d in data_rec:
    f.write(",".join("{}".format(n) for n in d.tolist()))
    f.write("\n")
f.close()

#                                           'fit_offset', 'fit_offset_min', 'fit_offset_max'])
#np.save('gui_warm_fit', data_rec)
#mags = [(mag_ranges[mag]['bright'] + mag_ranges[mag]['faint'])/2 for mag in mag_list]
#scales = [fit_info['semi'][mag]['bad_trak'][100]['fmod'].scale.val for mag in mag_list]
#scale_parmax = [fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmaxes[0] for mag in mag_list]
#scale_parmin = [fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmins[0] for mag in mag_list]
#offsets = [fit_info['semi'][mag]['bad_trak'][100]['fmod'].offset.val for mag in mag_list]
#offset_parmax = [fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmaxes[1] for mag in mag_list]
#offset_parmin = [fit_info['semi'][mag]['bad_trak'][100]['conf_orig'].parmins[1] for mag in mag_list]
#fig=plt.figure(figsize=(5, 3))
#plt.errorbar(mags, scales, [scale_parmax, np.array(scale_parmin)*-1.0], linestyle='')
#plt.xlabel('Catalog magnitude')
#plt.ylabel('Warm Pix Scale')
#plt.tight_layout()
#plt.grid()
#plt.savefig('bad_trak_warm_scale.png')
#fig=plt.figure(figsize=(5, 3))
#plt.errorbar(mags, offsets, [offset_parmax, np.array(offset_parmin)*-1.0], linestyle='')
#plt.xlabel('Catalog magnitude')
#plt.ylabel('Warm Pix Offset')
#plt.grid()
#plt.tight_layout()
#plt.savefig('bad_trak_warm_offset.png')

