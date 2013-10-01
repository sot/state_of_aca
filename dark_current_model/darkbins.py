"""Define bins for dark current distribution calculations"""
import numpy as np

x0 = 2.
x1 = 20000.
dx = 0.05
bins = 10**np.arange(np.log10(x0), np.log10(x1), dx)
