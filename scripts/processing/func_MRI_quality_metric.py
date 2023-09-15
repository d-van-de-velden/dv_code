#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
#
import os
import os.path as op
import numpy as np
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.viz.plot_tsnr import plot_tsnr


RAS_AXIS_ORDER = {"x": 0, "y": 1, "z": 2}

def gsr(epi_data, mask, direction="y", ref_file=None, out_file=None):
    # sourcery skip: raise-specific-error
    """
    ORIGINAL from mriqc
    Compute the :abbr:`GSR (ghost to signal ratio)` [Giannelli2010]_.

    The procedure is as follows:

    #. Create a Nyquist ghost mask by circle-shifting the original mask by :math:`N/2`.

    #. Rotate by :math:`N/2`

    #. Remove the intersection with the original mask

    #. Generate a non-ghost background

    #. Calculate the :abbr:`GSR (ghost to signal ratio)`


    .. warning ::

    This should be used with EPI images for which the phase
    encoding direction is known.

    :param str epi_file: path to epi file
    :param str mask_file: path to brain mask
    :param str direction: the direction of phase encoding (x, y, all)
    :return: the computed gsr

    """
    direction = direction.lower()
    if direction[-1] not in ["x", "y", "all"]:
        raise Exception(
            f"Unknown direction {direction}, should be one of x, -x, y, -y, all"
        )

    if direction == "all":
        result = []
        for newdir in ["x", "y"]:
            ofile = None
            if out_file is not None:
                fname, ext = op.splitext(ofile)
                if ext == ".gz":
                    fname, ext2 = op.splitext(fname)
                    ext = ext2 + ext
                ofile = "{0}_{1}{2}".format(fname, newdir, ext)
            result += [gsr(epi_data, mask, newdir, ref_file=ref_file, out_file=ofile)]
        return result

    # Roll data of mask through the appropriate axis
    axis = RAS_AXIS_ORDER[direction]
    n2_mask = np.roll(mask, mask.shape[axis] // 2, axis=axis)

    # Step 3: remove from n2_mask pixels inside the brain
    n2_mask = n2_mask * (1 - mask)

    # Step 4: non-ghost background region is labeled as 2
    n2_mask = n2_mask + 2 * (1 - n2_mask - mask)

    # Step 5: signal is the entire foreground image
    ghost = np.mean(epi_data[n2_mask == 1]) - np.mean(epi_data[n2_mask == 2])
    signal = np.median(epi_data[n2_mask == 0])

    return float(ghost / signal)


def motion_evaluation(fname_par=None, fd_threshold=0.3):
#
#
#
#
#
    fdir_qc = f'{os.path.dirname(fname_par)}/QC/'
    check_make_dir(fdir_qc)
    fname_qc = f'QC_motion_{os.path.basename(fname_par)[:-4]}.png'
    fname_fd = f'QC_FD1D_{os.path.basename(fname_par)[:-4]}.txt'

    motion_params = np.genfromtxt(fname_par).T

    rotations = np.transpose(np.abs(np.diff(motion_params[0:3, :])))
    translations = np.transpose(np.abs(np.diff(motion_params[3:6, :])))
    fd = np.sum(translations, axis=1) + \
        (50 * np.pi / 180) * np.sum(rotations, axis=1)
    fd = np.insert(fd, 0, 0)

    out_file = fdir_qc + fname_fd
    np.savetxt((out_file), fd)


    figure, axis = plt.subplots(3,1)
    figure.set_size_inches(10, 10)
    figure.set_dpi(100)
    figure.set_facecolor('grey')
    n = 7
    colors = plt.cm.jet(np.linspace(0,1,n))

    axis[0].plot(motion_params[0], color=colors[0], label='x', linewidth=2.5)
    axis[0].plot(motion_params[1], color=colors[1], label='y', linewidth=2.5)
    axis[0].plot(motion_params[2], color=colors[2], label='z', linewidth=2.5)
    axis[0].set_title("MCFLIRT estimated rotations (radians)")
    axis[0].set_xlabel('Scan volume [TR:2]')
    axis[0].set_facecolor("lightgrey")
    axis[0].legend(loc='upper right', ncol=1)

    axis[1].plot(motion_params[3], color=colors[3], label='x', linewidth=2.5)
    axis[1].plot(motion_params[4], color=colors[4], label='y', linewidth=2.5)
    axis[1].plot(motion_params[5], color=colors[5], label='z', linewidth=2.5)
    axis[1].set_title('MCFLIRT estimated translations (mm)')
    axis[1].set_xlabel('Scan volume [TR:2]')
    axis[1].set_facecolor("lightgrey")
    axis[1].legend(loc='upper right', ncol=1)                

    axis[2].plot(fd, color=colors[6], label='FD', linewidth=4)
    axis[2].axhline(fd_threshold, color='red', linestyle='--', lw=2, alpha=0.8)
    x = np.arange(fd.shape[0])
    trans = mtransforms.blended_transform_factory(axis[2].transData, axis[2].transAxes)
    axis[2].fill_between(x, 0, 1, where=fd >= fd_threshold,
                            facecolor='red', alpha=0.5, transform=trans)
    axis[2].set_title("Framewise displacement ")
    axis[2].set_xlabel('Scan volume [TR:2]')
    axis[2].set_facecolor("lightgrey")
    axis[2].legend(loc='upper right', ncol=1)


    figure.set_tight_layout('rect')
    fname_plot = (fdir_qc + fname_qc)
    plt.savefig(fname_plot)

    return fname_plot


def calc_tSNR(fname, params):
    
    tmp_fname_func = os.path.splitext(os.path.basename(fname))[0]
    fname_func = os.path.splitext(os.path.basename(tmp_fname_func))[0]

    idents  = fname_func.split('_')
    subjID  = idents[0]
    session = idents[1]
    run     = idents[2]

    fdir_func_r  = (params.get('fdir_proc_pre') 
                    + '/' + subjID
                    + '/' + session
                    + '/func'
                    )

    fdir_derivatives_anat = (params.get('fdir_proc_pre')
                            + '/' + subjID 
                            + '/' + session + '/anat'
                            )

    fname_func_mean = f'{fdir_func_r}/{fname_func}_mean.nii.gz'
    fname_func_std = f'{fdir_func_r}/{fname_func}_std.nii.gz'
    fname_func_tSNR = f'{fdir_func_r}/{fname_func}_tSNR.nii.gz'
    fname_func_r_tSNR = f'{fdir_func_r}/{fname_func}_tSNR_r.nii.gz'
    print(f'# Calcualting mean functional image \nFrom: {fname}\n   To: {fname_func_mean}')
    if os.system(f'FSL fslmaths {fname} -Tmean {fname_func_mean}') == 0: print('done successfully..')
    print(f'# Calcualting std functional image \nFrom: {fname}\n   To: {fname_func_std}')
    if os.system(f'FSL fslmaths {fname} -Tstd {fname_func_std}') == 0: print('done successfully..')
    
    print(f'# Calcualting tSNR for functional image \nFrom: mean / std \n   To: {fname_func_tSNR}')
    if os.system(f'FSL fslmaths {fname_func_mean} -div {fname_func_std} {fname_func_tSNR}') == 0: print('done successfully..')

    
    fname_T1w = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_nu.nii.gz'
    
    if os.path.exists(fname_T1w) == False:
        print(f'!! No anatomical derivatives present for [ {subjID} | {session} ] ')
    else:
        print('')
    fname_T1w_brain = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_nu_brain.nii'

    print(f'# Extract brain from T1w image \nFrom: {fname_T1w}\n   To: {fname_T1w_brain}')
    if os.system(f'FSL bet {fname_T1w} {fname_T1w_brain}') == 0: print('done successfully..')
    print(f'# Coregister functional MR to anatomical T1w \nFrom: {fname_func_tSNR}\n   To: {fname_func_r_tSNR}')
    if os.system(f'FSL flirt -in {fname_func_tSNR} -ref {fname_T1w} -out {fname_func_r_tSNR}') == 0: print('done successfully..')
    
    
    fname_parc = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_aparc+aseg.nii.gz'
    parc = nib.load(fname_parc)
    dat_aparc = np.array(parc.dataobj)
    func_tSNR = nib.load(fname_func_r_tSNR)
    dat_func_tSNR = np.array(func_tSNR.dataobj)
    dat_func_lh_tSNR = dat_func_tSNR
    dat_func_lh_tSNR[dat_aparc < 1000 ] = 0    
    dat_func_lh_tSNR[dat_aparc > 1037 ] = 0
    
    dat_func_tSNR = np.array(func_tSNR.dataobj)
    dat_func_rh_tSNR = dat_func_tSNR
    dat_func_rh_tSNR[dat_aparc < 2000 ] = 0
    dat_func_rh_tSNR[dat_aparc > 2037 ] = 0

    dat_func_tSNR2 = dat_func_lh_tSNR + dat_func_rh_tSNR
    dat_func_tSNR[dat_func_tSNR == 0] = 'nan'
    

    tSNR_median = round( np.nanmedian( dat_func_tSNR ), 2)
    print(f'The median tSNR in the grey matter is: {tSNR_median}')


    return fname_func_r_tSNR


def get_tSNR(fname_tsnr, fname_parc, params, doPlot=False):
#
#
#
#
    tmp_fname_tsnr = os.path.splitext(os.path.basename(fname_tsnr))[0]
    tmp_fname_tsnr = os.path.splitext(os.path.basename(tmp_fname_tsnr))[0]

    idents  = tmp_fname_tsnr.split('_')
    subjID  = idents[0]
    session = idents[1]
    
    parc = nib.load(fname_parc)
    dat_aparc = np.array(parc.dataobj)
    
    dat = nib.load(fname_tsnr).get_fdata()
    dat_func_lh_tSNR = dat
    dat_func_lh_tSNR[dat_aparc < 1000 ] = 0    
    dat_func_lh_tSNR[dat_aparc > 1037 ] = 0
    
    
    dat = nib.load(fname_tsnr).get_fdata()
    dat_func_rh_tSNR = dat
    dat_func_rh_tSNR[dat_aparc < 2000 ] = 0
    dat_func_rh_tSNR[dat_aparc > 2037 ] = 0

    dat_func_tSNR2 = dat_func_lh_tSNR + dat_func_rh_tSNR
    
    dat_func_tSNR2[dat_func_tSNR2 == 0] = 'nan'
    dat_func_rh_tSNR[dat_func_rh_tSNR == 0] = 'nan'
    dat_func_lh_tSNR[dat_func_lh_tSNR == 0] = 'nan'   


    dat = nib.load(fname_tsnr).get_fdata()
    dat[dat_aparc == 0] = 'nan'
    
    fdir_qc = f'{os.path.dirname(fname_tsnr)}/QC/'
    check_make_dir(fdir_qc)
    fname_qc = f'QC_tSNR_{tmp_fname_tsnr}.png'
    fname_plot = (fdir_qc + fname_qc)
    
    
    data = [dat, dat_func_tSNR2, dat_func_rh_tSNR, dat_func_lh_tSNR]
    ROI_label = ['wholebrain', 'grey matter','grey matter lh','grey matter rh' ]
    
    if doPlot == True:
        plot_tsnr(fname_tsnr, data, ROI_label, dat_idx=1)

    tSNR_median = []
    for jdat in data:
        tSNR_median.append(round( np.nanmedian( jdat ), 2))

    return fname_plot, tSNR_median, ROI_label


