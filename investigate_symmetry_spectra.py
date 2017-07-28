import imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.table import Table

imp.load_source('helper', '../Carbon-Spectra/helper_functions.py')
from helper import get_spectra_dr52, get_spectra_dr51, move_to_dir

# ------------------------------------------------------------------
# -------------------- Data reading and initial handling -----------
# ------------------------------------------------------------------
galah_data_dir = '/home/klemen/GALAH_data/'
spectra_dir_2 = '/media/storage/HERMES_REDUCED/dr5.2/'
galah_param = Table.read(galah_data_dir+'sobject_iraf_52_reduced.fits')
final_fits = 'symmetry_results.fits'
symmetry_results = Table.read(final_fits)

telurics_csv = 'flaglabels201707.csv'
tellutic_overlap = pd.read_csv(telurics_csv, sep=',', na_values='NULL')

# create subset of data
idx_possible_telluric = tellutic_overlap['Ca6508.8496'] > 0.05
tellutic_overlap_sub = tellutic_overlap[idx_possible_telluric]
symmetry_results = symmetry_results[np.in1d(galah_param['sobject_id'], tellutic_overlap_sub['sobject_id'])]

# plt.scatter(tellutic_overlap_sub['Ca6508.8496'], symmetry_results['Ca_6508.8496'], lw=0, s=1, c='black', alpha=0.2)
# plt.show()


# ------------------------------------------------------------------
# -------------------- Line settings -------------------------------
# ------------------------------------------------------------------
col_name = 'Ca_6508.8496'
symmetry_coefficients_bad = -0.1
symmetry_coefficients_ok = 0.0
symmetry_coefficients_bad2 = 0.1
symmetry_range = 0.02
i_band = 3

wvl_center = float(col_name.split('_')[1])
wvl_range = 1


# ------------------------------------------------------------------
# -------------------- Plots ---------------------------------------
# ------------------------------------------------------------------
idx_bad = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_bad - symmetry_range),
                                  symmetry_results[col_name] < (symmetry_coefficients_bad + symmetry_range)))[0]
idx_ok = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_ok - symmetry_range),
                                 symmetry_results[col_name] < (symmetry_coefficients_ok + symmetry_range)))[0]
idx_bad2 = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_bad2 - symmetry_range),
                                   symmetry_results[col_name] < (symmetry_coefficients_bad2 + symmetry_range)))[0]

print 'Found bad:', len(idx_bad)
print 'Found bad2:', len(idx_bad2)
print 'Found ok:', len(idx_ok)

move_to_dir('Plots')
n_plots = 300

plot_bad_cols = idx_bad[np.int64(np.random.random(n_plots)*len(idx_bad))]
plot_bad_cols2 = idx_bad2[np.int64(np.random.random(n_plots)*len(idx_bad2))]
plot_ok_cols = idx_ok[np.int64(np.random.random(n_plots)*len(idx_ok))]

for s_id in galah_param[plot_bad_cols]['sobject_id']:
    s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
    if len(s2) == 0:
        continue
    plt.plot(w2[0], s2[0]+0.1, color='red', alpha=0.1, lw=0.75)

for s_id in galah_param[plot_ok_cols]['sobject_id']:
    s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
    if len(s2) == 0:
        continue
    plt.plot(w2[0], s2[0]-0.2, color='black', alpha=0.1, lw=0.75)

for s_id in galah_param[plot_bad_cols]['sobject_id']:
    s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
    if len(s2) == 0:
        continue
    plt.plot(w2[0], s2[0]-0.5, color='blue', alpha=0.1, lw=0.75)

plt.ylim((0.1, 1.2))
plt.xlim((wvl_center-wvl_range, wvl_center+wvl_range))
plt.title('Red spectra: symmetry ~ {:0.2f}     black spectra: symmetry ~ {:0.2f}'.format(symmetry_coefficients_bad,
                                                                                         symmetry_coefficients_ok))
plt.ylabel('Flux')
plt.xlabel('Wavelength')
# plt.show()
plt.savefig(col_name+'_spectra.png', dpi=200)
