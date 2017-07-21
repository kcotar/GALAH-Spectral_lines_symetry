import os, imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.table import Table

imp.load_source('helper', '../Carbon-Spectra/helper_functions.py')
from helper import get_spectra_dr52, get_spectra_dr51, move_to_dir

GALAH_BLUE  = [4718, 4903]
GALAH_GREEN = [5649, 5873]
GALAH_RED   = [6481, 6739]
GALAH_IR    = [7590, 7890]


def get_band(wvl):
    if wvl >= GALAH_BLUE[0] and wvl <= GALAH_BLUE[1]:
        return 0
    elif wvl >= GALAH_GREEN[0] and wvl <= GALAH_GREEN[1]:
        return 1
    elif wvl >= GALAH_RED[0] and wvl <= GALAH_RED[1]:
        return 2
    elif wvl >= GALAH_IR[0] and wvl <= GALAH_IR[1]:
        return 3
    else:
        return np.nan

galah_data_dir = '/home/klemen/GALAH_data/'
spectra_data_dir = '/media/storage/HERMES_REDUCED/dr52/'
galah_param = Table.read(galah_data_dir+'sobject_iraf_52_reduced.fits')
abs_linelist = Table.read(galah_data_dir+'GALAH_Cannon_linelist.csv')  # Element, line_centre, lne_start, line_end

# get_spectra_dr52(spectra_data_dir, bands=[1, 2, 3, 4], root=spectra_data_dir, extension=4, individual=False)

print get_band(4251)