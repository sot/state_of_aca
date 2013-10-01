import numpy as np
from glob import glob
import pickle
import dark_models

method = 'levmar'
set_stat('cash')
# set_stat('chi2gehrels')
load_user_model(dark_models.smooth_broken_pow, 'sbp')
add_user_pars('sbp', ('gamma1', 'gamma2', 'x_b', 'x_r', 'ampl1'))

xmin = 20
xmax = 6000

fileglob = 'darkhist/200????.dat'
results = dict()
for filename in glob(fileglob):
    load_data(1, filename, 2)
    data = get_data()
    xall = data.x
    ignore(None, xmin)
    ignore(xmax, None)

    dark_models.xall = xall
    dark_models.imin = np.where(xall > xmin)[0][0]
    dark_models.imax = np.where(xall > xmax)[0][0]

    set_model(sbp)

    sbp.gamma1 = 0.16
    sbp.gamma2 = 3.15
    sbp.x_b = 125.
    sbp.x_r = 155.
    sbp.ampl1 = 12000.
    # freeze(sbp.gamma1)
    # freeze(sbp.gamma2)
    # freeze(sbp.x_b)
    freeze(sbp.x_r)

    fit()
    plot_fit()
    log_scale(XY_AXIS)
    projection()
    res = get_projection_results()
    results[filename] = dict((x, getattr(res, x)) for x in ('parnames', 'parmins', 'parvals', 'parmaxes'))
                             
# pickle.dump(results, open('fits.pickle', 'w'))

def smooth_broken_line(pars, x):
    """Smooth broken line function"""
    b1 = 3.5
    x_r = 2.5
    x_b = 2.5
    m1 = 0.0
    m2 = -2.0
    b2 = (m1-m2) * x_r + b1

    y = m1 * x + b1
    ok = x > x_b
    y[ok] = m2 * x[ok] + b2
    return y

