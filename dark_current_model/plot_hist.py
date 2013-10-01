import numpy as np
import pyfits
from glob import glob
import os
import re
from darkbins import x0, x1, dx, bins

clf()
# for darkdir in glob('aca_dark_cal/2??????'):
for darkdir in glob('aca_dark_cal/2008???'):
    print darkdir
    hdus = pyfits.open(os.path.join(darkdir,'imd.fits'))
    scale = 1.0
    if darkdir < 'aca_dark_cal/2006330':
        scale *= 0.61
    if darkdir < 'aca_dark_cal/2003120':
        scale *= 0.61
    dark = hdus[0].data.flatten() * scale
    dark += np.random.uniform(-0.5, 0.5, 1024**2)
    y, xb = histogram(dark, bins)
    x = (xb[1:] + xb[:-1]) / 2
    loglog(x, y)
##     out = open('darkhist/' + re.sub(r'aca_dark_cal/', '', darkdir) + '.dat', 'w')
##     for i in range(len(x)):
##         print >>out, x[i], y[i], 1+y[i]*0.1
##     out.close()
