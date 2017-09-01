import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.table import Table, join

imp.load_source('helper', '../Carbon-Spectra/helper_functions.py')
from helper import get_spectra_dr52, get_spectra_dr51, move_to_dir

# ------------------------------------------------------------------
# -------------------- Data reading and initial handling -----------
# ------------------------------------------------------------------
galah_data_dir = '/home/klemen/GALAH_data/'
spectra_dir_2 = '/media/storage/HERMES_REDUCED/dr5.2/'
galah_param = Table.read(galah_data_dir+'sobject_iraf_52_reduced.fits')
galah_abund = Table.read(galah_data_dir+'sobject_iraf_cannon_1.2.fits')

final_symmetry_fits = 'results_symmetry.fits'
final_ew_fits = 'results_ew.fits'
symmetry_results_all = Table.read(final_symmetry_fits)
ew_results_all = Table.read(final_ew_fits)

move_to_dir('Plots')

# join galah and subset of all data
galah_data = join(galah_param, galah_abund, keys='sobject_id', join_type='inner')
idx_sym_cols = np.in1d(symmetry_results_all['sobject_id'], galah_data['sobject_id'])
symmetry_results = symmetry_results_all[idx_sym_cols]
ew_results = ew_results_all[idx_sym_cols]

# prints to asses connection between datasets
# print galah_data[:20]['sobject_id']
# print symmetry_results[:20]['sobject_id']
# print galah_data[-20:]['sobject_id']
# print symmetry_results[-20:]['sobject_id']

abund_lines = symmetry_results.colnames[1:]
abund_cols = [col for col in galah_abund.colnames if '_abund' in col and not '_e_' in col]
for a_col in abund_cols:
    print a_col
    element = (a_col.split('_')[0]).capitalize()
    # get all lines for this element
    element_lines = [line for line in abund_lines if element in line]
    for abs_line in element_lines:
        print abs_line
        plt.title(a_col)
        c_data = galah_data[a_col]
        plt.scatter(ew_results[abs_line], symmetry_results[abs_line], lw=0, s=0.2, c=c_data, cmap='jet',
                    vmin=np.percentile(c_data, 1), vmax=np.percentile(c_data, 99))
        plt.xlim((np.nanpercentile(ew_results_all[abs_line], 1), np.nanpercentile(ew_results_all[abs_line], 99)))
        plt.ylim((np.nanpercentile(symmetry_results_all[abs_line], 1), np.nanpercentile(symmetry_results_all[abs_line], 99)))
        plt.xlabel('Line EW value')
        plt.ylabel('Symmetry coefficient value')
        plt.colorbar()
        plt.savefig(abs_line + '_ew_'+a_col+'.png', dpi=300)
        plt.close()
