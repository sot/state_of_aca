from matplotlib.pyplot import *
from matplotlib import cm
from astropy.io import fits

if not 'dark' in globals():
    dark = fits.open('2013191_dark.fits')[0].data
    star = fits.open('obs890_adat41.fits')[1].data

# 237 == 1
# 279 == 100
# dt = 42 degC => 10*2 => 21 degC per decade

star_img = star[400]['img_corr'] / 1.7  # e-/sec from a 1.7 sec readout
dark_img = dark[400:440, 400:440]  # e-/sec from calibration
dark_t0 = -15.31
star_mag = 10.3

rc("axes", labelsize=10, titlesize=10)
rc("xtick", labelsize=10)
rc("ytick", labelsize=10)
rc("font", size=10)

def draw_star(t=-16.0, mag=10.3):
    clf()
    img = dark_img.copy() * 10 ** ((t - dark_t0) / 21.0)
    nn = dark_img.shape[0] / 2
    img[nn -3:nn + 3, nn -3:nn + 3] += star_img * 10 ** ((star_mag - mag) * 0.4)

    imshow(img, interpolation='nearest', cmap=cm.gray)
    colorbar()


def draw_star_dark_range():
    figure(figsize=(3.5, 3.5))
    for t in [-15, -10, -5, 0]:
        for mag in [10.6, 10.0]:
            draw_star(t, mag)
            title('T = {} degC, Mag_ACA = {}'.format(t, mag))
            ax = gca()
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            tight_layout()
            savefig('star_img_mag{}_T{}.png'.format(mag, t))

def draw_star_dark2000():
    figure(figsize=(3.5, 3.5))
    t = -20.3
    mag = 10.6

    draw_star(t, mag)
    title('T = -10 degC, Mag_ACA = {}'.format(t, mag))
    ax = gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    tight_layout()
    savefig('star_img_2000_mag{}_T{}.png'.format(mag, t))
