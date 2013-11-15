import matplotlib.pyplot as plt

from Chandra.Time import DateTime
from Ska.engarchive import fetch_sci as fetch
from Ska.Matplotlib import plot_cxctime
from kadi import events

events.safe_suns.interval_pad = (10000, 200000)
events.normal_suns.interval_pad = (10000, 200000)
bad_intervals = events.safe_suns | events.normal_suns


def plot_aca_ccd_mission():
    plt.close(1)
    fig = plt.figure(1, figsize=(6, 4))
    dat = fetch.Msid('aacccdpt', '2000:001', stat='daily')
    dat.remove_intervals(bad_intervals)
    plot_cxctime(dat.times, dat.mins, '-b', fig=fig)
    plot_cxctime(dat.times, dat.maxes, '-b', fig=fig)
    plt.grid()
    plt.legend(loc='upper right', fontsize=12)
    plt.ylabel('Temperature (C)')
    plt.title('ACA CCD temperature')
    plt.tight_layout()
    plt.savefig('aca_ccd_mission.png')


def plot_aca_ccd_cases():
    plt.close(1)
    fig = plt.figure(1, figsize=(6, 4))
    dat = fetch.Msid('aacccdpt', DateTime() - 4 * 365, stat='daily')
    dat.remove_intervals(bad_intervals)
    plot_cxctime(dat.times, dat.vals, '-b', fig=fig)
    plot_cxctime(DateTime(['2012:001', '2018:001']).secs, [-19, -7], '-r', lw=2.5, label='-7 C')
    plot_cxctime(DateTime(['2016:001', '2018:001']).secs, [-11, -11], '-g', lw=2.5, label='-11 C')
    plot_cxctime(DateTime(['2014:185', '2018:001']).secs, [-14, -14], '-b', lw=2.5, label='-14 C')
    plt.legend(loc='upper left', fontsize=12)
    plt.grid()
    plt.ylabel('Temperature (C)')
    plt.title('ACA CCD temperature cases')
    plt.tight_layout()
    plt.savefig('aca_ccd_cases.png')
