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

final_symmetry_fits = 'results_symmetry.fits'
final_ew_fits = 'results_ew.fits'
symmetry_results = Table.read(final_symmetry_fits)
ew_results = Table.read(final_ew_fits)

telurics_csv = 'flaglabels201707.csv'
tellutic_overlap = pd.read_csv(telurics_csv, sep=',', na_values='NULL')

move_to_dir('Plots')

# print spectra of all lines at peak of the symmetry coefficient
# symmetry_range = 0.02
# wvl_range = 1
# n_plots = 200
# for col_name in symmetry_results.colnames[1:]:
#     print col_name
#     hist_val, hist_bin = np.histogram(symmetry_results[col_name], range=(-1, 1), bins=200)
#     max_symmetry = hist_bin[np.argmax(hist_val)]
#     idx_use = np.where(np.logical_and(symmetry_results[col_name] > (max_symmetry - symmetry_range),
#                                       symmetry_results[col_name] < (max_symmetry + symmetry_range)))[0]
#     plot_use_cols = idx_use[np.int64(np.random.random(n_plots)*len(idx_use))]
#     wvl_center = float(col_name.split('_')[1])
#     for s_id in galah_param[plot_use_cols]['sobject_id']:
#         i_band = int(wvl_center/1000.)-3
#         s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
#         if len(s2) == 0:
#             continue
#         plt.plot(w2[0], s2[0], color='black', alpha=0.1, lw=0.75)
#     plt.ylim((0.1, 1.2))
#     plt.xlim((wvl_center - wvl_range, wvl_center + wvl_range))
#     plt.axvline(x=wvl_center, lw=1, color='black', ls='dashed')
#     plt.title('Spectra at peak symmetry')
#     plt.ylabel('Flux')
#     plt.xlabel('Wavelength')
#     plt.savefig(col_name + '_spectra_central.png', dpi=300)
#     plt.close()
#
# raise SystemExit

# ------------------------------------------------------------------
# -------------------- Line settings -------------------------------
# ------------------------------------------------------------------
col_name = 'Ti_6716.6660'
symmetry_coefficients_bad = 0.0
symmetry_coefficients_ok = -1.
symmetry_coefficients_bad2 = 0.
symmetry_range = 0.1

wvl_center = float(col_name.split('_')[1])
i_band = int(wvl_center/1000)-3
wvl_range = 1

# subset of data
idx_sobjects = np.logical_and(np.logical_and(galah_param['teff_guess'] < 4800., galah_param['feh_guess'] < -1.5),
                              ew_results[col_name] > 0.01)
# idx_sobjects = np.logical_and(galah_param['teff_guess'] > 6500., ew_results[col_name] > 0.04)
# idx_sobjects = ew_results[col_name] > 0.05
# idx_sobjects = galah_param['teff_guess'] < 6000.
symmetry_results = symmetry_results[idx_sobjects]
galah_param = galah_param[idx_sobjects]


# ------------------------------------------------------------------
# -------------------- Plots ---------------------------------------
# ------------------------------------------------------------------
idx_bad = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_bad - symmetry_range),
                                  symmetry_results[col_name] < (symmetry_coefficients_bad + symmetry_range)))[0]
idx_ok = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_ok - symmetry_range),
                                 symmetry_results[col_name] < (symmetry_coefficients_ok + symmetry_range)))[0]
# idx_bad2 = np.where(np.logical_and(symmetry_results[col_name] > (symmetry_coefficients_bad2 - symmetry_range),
#                                    symmetry_results[col_name] < (symmetry_coefficients_bad2 + symmetry_range)))[0]

print 'Found bad:', len(idx_bad)
# print 'Found bad2:', len(idx_bad2)
print 'Found ok:', len(idx_ok)

n_plots = 300

plot_bad_cols = idx_bad[np.int64(np.random.random(n_plots)*len(idx_bad))]
# plot_bad_cols2 = idx_bad2[np.int64(np.random.random(n_plots)*len(idx_bad2))]
plot_ok_cols = idx_ok[np.int64(np.random.random(n_plots)*len(idx_ok))]

# for s_id in galah_param[plot_bad_cols_2]['sobject_id']:
#     s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
#     if len(s2) == 0:
#         continue
#     plt.plot(w2[0], s2[0]+0.1, color='blue', alpha=0.1, lw=0.75)

s_ids_ok = galah_param[plot_ok_cols]['sobject_id']
for s_id in galah_param[plot_ok_cols]['sobject_id']:
    s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
    if len(s2) == 0:
        continue
    plt.plot(w2[0], s2[0], color='black', alpha=0.1, lw=0.75)
print ','.join([str(s) for s in s_ids_ok[:30]])

# for s_id in galah_param[plot_bad_cols]['sobject_id']:
#     s2, w2 = get_spectra_dr52(str(s_id), bands=[i_band], root=spectra_dir_2, individual=False, extension=4)
#     if len(s2) == 0:
#         continue
#     plt.plot(w2[0], s2[0]-0.3, color='red', alpha=0.1, lw=0.75)

plt.axvline(x=wvl_center, lw=1, color='black', ls='dashed')
plt.ylim((0.1, 1.2))
plt.xlim((wvl_center-wvl_range, wvl_center+wvl_range))
plt.title('Red spectra: symmetry ~ {:0.2f}     black spectra: symmetry ~ {:0.2f}'.format(symmetry_coefficients_bad,
                                                                                         symmetry_coefficients_ok))
plt.ylabel('Flux')
plt.xlabel('Wavelength')
# plt.show()
# plt.savefig(col_name+'_spectra.png', dpi=300)
plt.savefig(col_name+'_feh-low_teff-low.png', dpi=300)
# plt.savefig(col_name+'_spectra_offset.png', dpi=300)

