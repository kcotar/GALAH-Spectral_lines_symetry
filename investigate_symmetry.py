import os, imp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.table import Table

imp.load_source('helper', '../Carbon-Spectra/helper_functions.py')
from helper import get_spectra_dr52, get_spectra_dr51, move_to_dir

# ------------------------------------------------------------------
# -------------------- Constants and settings ----------------------
# ------------------------------------------------------------------

GALAH_BLUE  = [4710, 4910]
GALAH_GREEN = [5640, 5880]
GALAH_RED   = [6480, 6740]
GALAH_IR    = [7590, 7890]

# ------------------------------------------------------------------
# -------------------- Functions -----------------------------------
# ------------------------------------------------------------------


def get_band(wvl):
    # based on wvl determine in which color arm does it lie
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


def get_line_halfwidth(w_c, w_b, w_e):
    # get minimum width of the spectral absorbtion line
    return min(w_c - w_b, w_e - w_c)


def _is_boundary(wvl_obs, wvl_bin, wvl_target):
    return np.logical_and(wvl_target < (wvl_obs + wvl_bin / 2.),
                          wvl_target > (wvl_obs - wvl_bin / 2.))


def _determine_bin_percentage(wvl_obs, wvl_bin, wvl_target, lower=False, upper=False):
    if lower:
        return (wvl_target - (wvl_obs - wvl_bin/2.)) / wvl_bin
    if upper:
        return ((wvl_obs + wvl_bin/2.) - wvl_target) / wvl_bin


def _ew_partial(flux):
    return 1. - flux  # not the best option for noisy data, shallow lines and continuum problems


def _integrate(spectra, wvl, wvl_bin, wvl_begin, wvl_end):
    ew_integral = 0
    stage = 0
    for i_p in range(len(spectra)):
        wvl_obs = wvl[i_p]
        if stage == 0:
            if _is_boundary(wvl_obs, wvl_bin, wvl_begin):
                boundary_integ = (_ew_partial(spectra[i_p]) * wvl_bin) * _determine_bin_percentage(wvl_obs, wvl_bin, wvl_begin, upper=True)
                ew_integral += boundary_integ
                stage = 1
            else:
                continue
        elif stage == 1:
            if _is_boundary(wvl_obs, wvl_bin, wvl_end):
                boundary_integ = (_ew_partial(spectra[i_p]) * wvl_bin) * _determine_bin_percentage(wvl_obs, wvl_bin, wvl_end, lower=True)
                ew_integral += boundary_integ
                break
            else:
                ew_integral += _ew_partial(spectra[i_p]) * wvl_bin
    # TODO: how to handle negative EW integrals
    return ew_integral


def integrate_line_halves(spectra, wvl, line_center, line_halfwidth):
    wvl_bin = wvl[1] - wvl[0]
    if (line_center - line_halfwidth) < wvl[0] or (line_center + line_halfwidth) > wvl[-1]:
        return np.nan, np.nan
    idx_use = np.logical_and(wvl >= (line_center - line_halfwidth - wvl_bin),
                             wvl <= (line_center + line_halfwidth + wvl_bin))
    spectra_sub = spectra[idx_use]
    wvl_sub = wvl[idx_use]
    return (_integrate(spectra_sub, wvl_sub, wvl_bin, line_center - line_halfwidth, line_center),
            _integrate(spectra_sub, wvl_sub, wvl_bin, line_center, line_center + line_halfwidth))

# ------------------------------------------------------------------
# -------------------- Data reading and initial handling -----------
# ------------------------------------------------------------------

galah_data_dir = '/home/klemen/GALAH_data/'
spectra_data_dir = '/media/storage/HERMES_REDUCED/dr5.2/'
galah_param = Table.read(galah_data_dir+'sobject_iraf_52_reduced.fits')
abs_linelist = Table.read(galah_data_dir+'GALAH_Cannon_linelist.csv')  # Element, line_centre, lne_start, line_end

# test subsets
# galah_param = galah_param[5000:10000]

# Create a Table that will hold results of the analysis
abs_line_names = [str(line['Element'])+'_{:.4f}'.format(line['line_centre']) for line in abs_linelist]
empty_nan_line = [np.nan for i_l in range(len(abs_linelist))]

ew_results = Table(np.zeros((len(galah_param), len(abs_linelist))), names=abs_line_names)
symmetry_results = Table(np.zeros((len(galah_param), len(abs_linelist))), names=abs_line_names)

final_symmetry_fits = 'results_symmetry.fits'
final_ew_fits = 'results_ew.fits'
if not os.path.isfile(final_symmetry_fits):
    for i_g in range(len(galah_param)):
        galah_object = galah_param[i_g]
        print 'Reading data for sobject_id: ', galah_object['sobject_id']
        spectra, wavelength = get_spectra_dr52(str(galah_object['sobject_id']), bands=[1, 2, 3, 4], root=spectra_data_dir)
        if len(spectra) == 0:
            symmetry_results[i_g] = empty_nan_line
            continue
        # print ' - working on spectral absorption lines'
        for i_l in range(len(abs_linelist)):
            absorption_line = abs_linelist[i_l]
            absorption_line_halfwidth = get_line_halfwidth(absorption_line['line_centre'],
                                                           absorption_line['line_start'], absorption_line['line_end'])
            i_band = get_band(absorption_line['line_centre'])
            ew_left, ew_right = integrate_line_halves(spectra[i_band], wavelength[i_band],
                                                      absorption_line['line_centre'], absorption_line_halfwidth)
            # compute observed parameters
            ew = (ew_left + ew_right)
            asymetry = (ew_left - ew_right) / (ew_left + ew_right)
            # save to final table(s)
            ew_results[abs_line_names[i_l]][i_g] = ew
            symmetry_results[abs_line_names[i_l]][i_g] = asymetry

    # add sobject_id column to results
    symmetry_results.add_column(galah_param['sobject_id'], index=0)
    # save results
    if os.path.isfile(final_symmetry_fits):
        os.remove(final_symmetry_fits)
        os.remove(final_ew_fits)
    symmetry_results.write(final_symmetry_fits)
    ew_results.write(final_ew_fits)
else:
    print 'Reading precomputed data'
    symmetry_results = Table.read(final_symmetry_fits)
    ew_results = Table.read(final_ew_fits)

move_to_dir('Plots')
# plot results as histograms
print 'Plotting graphs'
# txt = open('distribution_summary.csv', 'w')
for abs_line in abs_line_names:
    print abs_line
    plt.title(abs_line)
    median_hist = np.nanmedian(symmetry_results[abs_line])
    hist_vals, hist_bins, _ = plt.hist(symmetry_results[abs_line], bins=300, range=(median_hist - 0.15, median_hist + 0.15))
    # hist_vals, hist_bins, _ = plt.hist(symmetry_results[abs_line], bins=300, range=(-2, 2.))
    # plt.xlim((-2, 2.))
    # plt.ylim((0, 35000))
    plt.xlabel('Symmetry coefficient value')
    plt.ylabel('Symmetry distribution')
    plt.savefig(abs_line+'_4.png', dpi=200)
    plt.close()
    # txt.write(abs_line+' '+str(np.max(hist_vals))+' '+str(np.nanstd(symmetry_results[abs_line]))+' '+str(np.nanmedian(symmetry_results[abs_line]))+'\n')
# txt.close()
