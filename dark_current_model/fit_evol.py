import time
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pickle
import dark_models
from sherpa import ui

sbp = None  # for pychecker

method = 'levmar'
ui.set_stat('cash')
ui.load_user_model(dark_models.smooth_broken_pow, 'sbp')
ui.add_user_pars('sbp', ('gamma1', 'gamma2', 'x_b', 'x_r', 'ampl1'))


def fit_evol(fileglob='darkhist/20?????.dat', outfile=None, xmin=35, xmax=6000):
    results = dict()
    for i, filename in enumerate(glob(fileglob)):
        print "\n\n*************** {} *****************".format(filename)
        plt.figure(1)
        ui.load_data(1, filename, 2)
        data = ui.get_data()
        xall = data.x
        plt.loglog(data.x, data.y)
        ui.ignore(None, xmin)
        ui.ignore(xmax, None)

        dark_models.xall = xall
        dark_models.imin = np.where(xall > xmin)[0][0]
        dark_models.imax = np.where(xall > xmax)[0][0]

        ui.set_model(sbp)
        ui.set_method('simplex')

        sbp.gamma1 = 0.16
        sbp.gamma2 = 3.15
        sbp.x_b = 125.
        sbp.x_r = 155.
        sbp.ampl1 = 12000.
        ui.freeze(sbp.gamma1)
        ui.freeze(sbp.gamma2)
        ui.freeze(sbp.x_b)
        ui.freeze(sbp.x_r)

        ui.fit()
        # ui.plot_fit()
        # ui.set_xlog()
        # ui.set_ylog()
        ui.set_conf_opt('numcores', 1)
        ui.conf()
        res = ui.get_conf_results()
        results[filename] = dict((x, getattr(res, x))
                                 for x in ('parnames', 'parmins', 'parvals', 'parmaxes'))
    if outfile:
        pickle.dump(results, open(outfile, 'w'))
    return results
