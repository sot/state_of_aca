import re
import os
from glob import glob
import pickle

import matplotlib.pyplot as plt
import numpy as np
from sherpa import ui

import dark_models

sbp = None  # for pychecker
g1 = None

method = 'levmar'
ui.set_stat('cash')
ui.set_method('simplex')
ui.load_user_model(dark_models.smooth_broken_pow, 'sbp')
ui.add_user_pars('sbp', ('gamma1', 'gamma2', 'x_b', 'x_r', 'ampl1'))


def fit_gauss_sbp():
    g1 = ui.gauss1d.g1
    ui.set_model(sbp + g1)
    ui.set_method('simplex')

    g1.fwhm = 5.0
    g1.pos = 7.0
    g1.ampl = 30000.
    ui.freeze(sbp.gamma1)
    ui.freeze(sbp.gamma2)
    ui.freeze(sbp.x_b)
    ui.freeze(sbp.x_r)
    ui.freeze(g1.fwhm)
    ui.freeze(g1.pos)
    ui.thaw(g1.ampl)
    ui.fit()

    ui.thaw(g1.fwhm)
    ui.thaw(g1.pos)
    ui.fit()

    ui.thaw(sbp)
    ui.freeze(sbp.x_r)
    ui.fit()


def fit_sbp():
    ui.set_model(sbp)

    ui.thaw(sbp)
    ui.freeze(sbp.x_r)
    ui.freeze(sbp.gamma1)
    ui.fit()


def fit_evol(dateglob='20?????', rootdir='darkhist_peaknorm', outroot='', xmin=25.0, xmax=4000,
             conf=True, gauss=False):
    results = {}
    fileglob = os.path.join(rootdir, '{}.dat'.format(dateglob))

    for i, filename in enumerate(glob(fileglob)):
        filedate = re.search(r'(\d{7})', filename).group(1)
        print "\n\n*************** {} *****************".format(filename)
        plt.figure(1)
        ui.load_data(1, filename, 2)
        data = ui.get_data()
        ui.ignore(None, xmin)
        ui.ignore(xmax, None)

        dark_models.xall = data.x
        # dark_models.imin = np.where(xall > xmin)[0][0]
        # dark_models.imax = np.where(xall > xmax)[0][0]

        sbp.gamma1 = 0.05
        sbp.gamma2 = 3.15
        sbp.gamma2.min = 2.
        sbp.gamma2.max = 4.
        sbp.x_b = 130.
        sbp.x_b.min = 100.
        sbp.x_b.max = 160.
        sbp.x_r = 50.
        ok = (data.x > 40) & (data.x < 60)
        sbp.ampl1 = np.mean(data.y[ok])

        if gauss:
            fit_gauss_sbp()
        else:
            fit_sbp()

        pars = (sbp.gamma1.val, sbp.gamma2.val, sbp.x_b.val, sbp.x_r.val, sbp.ampl1.val)
        fit_y = dark_models.smooth_broken_pow(pars, data.x)

        if conf:
            ui.set_conf_opt('numcores', 1)
            ui.conf()
            res = ui.get_conf_results()
            result = dict((x, getattr(res, x))
                          for x in ('parnames', 'parmins', 'parvals', 'parmaxes'))
            result['x'] = data.x
            result['y'] = data.y
            result['y_fit'] = fit_y
            results[filedate] = result

        if outroot is not None:
            ui.notice(0, xmax)
            ui.set_xlog()
            ui.set_ylog()
            ui.plot_fit()
            plt.xlim(1, 1e4)
            plt.ylim(0.5, 1e5)
            plt.grid(True)
            plt.xlabel('Dark current (e-/sec)')
            outfile = os.path.join(rootdir, '{}{}.png'.format(outroot, filedate))
            print 'Writing', outfile
            plt.savefig(outfile)

            if conf:
                outfile = os.path.join(rootdir, '{}{}.pkl'.format(outroot, filedate))
                print 'Writing', outfile
                pickle.dump(result, open(outfile, 'w'), protocol=-1)

    if outroot is not None:
        outfile = os.path.join(rootdir, '{}fits.pkl'.format(outroot))
        print 'Writing', outfile
        pickle.dump(results, open(outfile, 'w'), protocol=-1)

    return results
