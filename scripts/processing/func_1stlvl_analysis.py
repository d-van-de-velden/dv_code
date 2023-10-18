#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
###############################################################

import json
import os
from dv_code.scripts.processing.func_MRI_quality_metric import calc_tSNR_aligned, read_table_tSNR

import nibabel as nib
import numpy as np
import pandas as pd
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.misc.analysis_feedback import loadingBar
from dv_code.scripts.misc.select_filetype import get_only_nifti
from nilearn import plotting
from nilearn.glm.first_level import (FirstLevelModel,
                                    make_first_level_design_matrix)
from nilearn.glm import threshold_stats_img
import matplotlib.pyplot as plt
from nilearn.image import binarize_img
from nilearn.maskers import NiftiMasker
import seaborn as sns


def func_apply_glm_treshholded(participants=None, params=None, smoothing_fwhm=8, tSNR_tresh=30, contrast=None):
    
    """_summary_
    Iterate across subjects, sessions, runs. Calculate tSNR and only use run for GLM 
    if it passes given tSNR threshhold. Calculate all contrast maps and average across subjects.
    Save all contrast maps as nifti for each subject, session. 
    Average all contrast maps (niftis) of subjects, per session.
    """
    
    if participants is None:
        participants = []
    check_make_dir(params.get('fdir_proc_pre'))
        
    print(
        f'####\nPreprocessing of functional MR data is initiated for: {len(participants)} subjects\n####'
    )

    all_sessions = list()
    runs_count = 0
    for iSubj in range(len(participants)):
        print(f'# Identifying sessions from "derivatives/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_proc_pre') + '/' + participants[iSubj]

        if os.path.exists(tmp_fdir_subj):
            sessions = os.listdir(tmp_fdir_subj)
            for session in sessions:
                print(f'  # Performing calculation on session:   {session}')
                all_sessions.append(sessions)
                fdir_derivatives_func = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/func'
                                        )
                fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/anat'
                                        )
                if os.path.exists(fdir_derivatives_func) and os.path.exists(fdir_derivatives_anat) :
                    runs = os.listdir(fdir_derivatives_func)
                    runs = get_only_nifti(runs, 1)
                    
                    proc_lvl = ['_st_mcf_topUP_r_w.nii.gz']

                    runs[:] = [proc_stage for proc_stage in runs if any(item in proc_stage for item in proc_lvl)]
                    runs.sort()
                    if len(runs) > 0:
                        for run in runs:
                            
                            # Grab information from filename
                            indices = run.split('_')
                            subjID = indices[0]
                            ses = indices[1]
                            runID = indices[2]
                            print(f'    # Performing calculation on run:   {runID}')
                            run_reduced = run.replace(proc_lvl[0], "")
                            
                            # Get functional data as nifti file
                            tmp_fname_final_funcMR = fdir_derivatives_func + '/' + run
                            if os.path.exists(tmp_fname_final_funcMR):
                                print(f'     Found functional data :   {tmp_fname_final_funcMR}.')
                                                                
                            # Perform image quality metric analysis
                            fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                + '/' + subjID 
                                + '/' + session + '/anat'
                                )
                            fname_parc = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_aparc+aseg_w.nii.gz'
    
                            # read the table containing the tSNR information of given 
                            # functional MR file
                            tSNR_dict = read_table_tSNR(tmp_fname_final_funcMR, params)
                            
                            if len(tSNR_dict) == 0:
                                # if no such file exists, calculate tSNR
                                tmp_tSNR_median = calc_tSNR_aligned(tmp_fname_final_funcMR, fname_parc, params)
                            else:
                                tmp_tSNR_median = float(tSNR_dict[0]['tSNR'])
                                
                            if tmp_tSNR_median >= tSNR_tresh:
                                
                                # Count the runs that surpassed the tSNR threshold and pass 
                                # that information to the final session average filename 
                                runs_count = runs_count + 1
                            
                                # Get motion parameters for functional data as par file
                                tmp_fname_motion_funcMR = fdir_derivatives_func + '/' + run_reduced + '_st_mcf.par'
                                if os.path.exists(tmp_fname_motion_funcMR):
                                    print(f'     Found motion parameter data :   {tmp_fname_motion_funcMR}.')
                                
                                
                                # Get json file to functional data
                                tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]
                                tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'
                                tmp_fname_raw_funcMR_json = (tmp_fdir_func
                                                            + '/' + run_reduced + '.json')
                                if os.path.exists(tmp_fname_raw_funcMR_json):
                                    print(f'     Found .json file for functional data :   {tmp_fname_raw_funcMR_json}.')
                                with open(f'{tmp_fname_raw_funcMR_json}') as f:
                                    func_json = json.load(f)
                                RepetitionTime = func_json.get('RepetitionTime')
                                print(f'        Found TR of : {RepetitionTime}s for functional data.')
                                
                                # Get event file to functional data
                                eventspath = (tmp_fdir_func 
                                            + '/' + run_reduced
                                            + '.tsv')
                                if os.path.exists(eventspath):
                                    print(f'     Found evt data for functional data : {eventspath}.')
                                

                                print('     # Load data for processing functional data for stimuli...')
                                # Load data
                                func    = nib.load(tmp_fname_final_funcMR)
                                tr      = RepetitionTime
                                n_scans = func.shape[3]
                                
                                events  = pd.read_csv(eventspath,sep='\t')
                                
                                baseline = {'onset': [events.onset[47] + events.duration[47]],
                                            'duration': [( n_scans * tr ) - events.onset[47] + events.duration[47]],
                                            'trial_type': ['baseline']}
                                df_baseline = pd.DataFrame(baseline)
                                
                                events = pd.concat([events, df_baseline], ignore_index=True, sort=False)

    
                                motion_par = ["tx", "ty", "tz", "rx", "ry", "rz"]
                                motion = np.array( pd.read_csv(tmp_fname_motion_funcMR,
                                                    delim_whitespace=True,
                                                    names=motion_par) 
                                                )

                                func_fit_glm(func, tr, events, motion, motion_par, 
                                            subjID, smoothing_fwhm, params, contrast,
                                            tSNR_tresh, ses, runID)

                        
                        
                        # Average all contrasts and runs
                        fdir_der_firstlvl_final = ( params.get('fdir_proc') 
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/'
                            )
                        check_make_dir(fdir_der_firstlvl_final)
                        print("Average for each contrast over all runs")
                        for i, contrast_id in enumerate(contrast):
                            loadingBar(i, len(contrast), contrast_id)
                            str_contrast = contrast_id.replace(" ", "_")

                            image_array = list()
                            for j, run in enumerate(runs):
                                # Identify identifier
                                indices = run.split('_')
                                subjID = indices[0]
                                ses = indices[1]
                                runID = indices[2]
                                
                                fdir_der_firstlvl = ( params.get('fdir_proc') 
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/' + runID + '/'
                            )
                                fname_zmap = fdir_der_firstlvl + f"{subjID}_{ses}_{runID}_{str_contrast}_tSNR{tSNR_tresh}_zscore.nii"
                                if os.path.exists(fname_zmap):
                                    tmp_zmap   = nib.load(fname_zmap)
                                    image_array.append( tmp_zmap.get_fdata() )
                                                            
                            
                            if len(image_array) > 0:
                                image_array = np.array(image_array)
                                avg_zmap = np.median(image_array, axis=0)
                                fname_avg_zmap = fdir_der_firstlvl_final + f"{subjID}_{ses}_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                                
                                avg_zmap_img = nib.Nifti1Image(avg_zmap, func.affine)
                                nib.save(img=avg_zmap_img, filename=fname_avg_zmap)
                                
                                sd_zmap = np.std(image_array, axis=0)
                                fname_sd_zmap = fdir_der_firstlvl_final + f"{subjID}_{ses}_sd_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                                
                                sd_zmap_img = nib.Nifti1Image(sd_zmap, func.affine)
                                nib.save(img=sd_zmap_img, filename=fname_sd_zmap)



    # Average all contrasts and runs of all sessions across subjects 
    all_sessions = list( np.unique(all_sessions))
    for session in all_sessions:
        print(session)
        fdir_group_firstlvl_final = ( params.get('fdir_proc') 
                                + '/1st_level/Group/' 
                                + session + '/'
                            )
        check_make_dir(fdir_group_firstlvl_final)
        for i, contrast_id in enumerate(contrast):
            loadingBar(i, len(contrast), contrast_id)
            str_contrast = contrast_id.replace(" ", "_")
            
            ALL_subj_image_array = list()
            subj_count = 0
            for j, subjID in enumerate(participants):
                tmp_fdir_firstlvl = ( params.get('fdir_proc') 
                                        + '/1st_level/'
                                        + subjID + '/' + ses + '/'
                                    )
                if os.path.exists(tmp_fdir_firstlvl):
                    fname_zmap2 = tmp_fdir_firstlvl + f"{subjID}_{ses}_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                    if os.path.exists(fname_zmap2):
                        tmp_zmap2   = nib.load(fname_zmap2)
                        ALL_subj_image_array.append( tmp_zmap2.get_fdata() )
                        subj_count = subj_count +1
                    
            
            # Z scores              
            ALL_subj_image_array = np.array(ALL_subj_image_array)
            all_avg_zmap = np.mean(ALL_subj_image_array, axis=0)
            fname_avg_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_Nruns_{runs_count}_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                            
            all_avg_zmap_img = nib.Nifti1Image(all_avg_zmap, func.affine)
            nib.save(img=all_avg_zmap_img, filename=fname_avg_zmap)
            
            #all_sd_zmap = np.std(ALL_subj_image_array, axis=0)
            #fname_sd_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_Nruns_{runs_count}_sd_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                            
            #all_sd_zmap_img = nib.Nifti1Image(all_sd_zmap, func.affine)
            #nib.save(img=all_sd_zmap_img, filename=fname_sd_zmap)
            
            vol_data_init   = nib.load(fname_avg_zmap)
            all_avg_zmap_FDR, threshold = threshold_stats_img(stat_img=vol_data_init, alpha=0.05,
                                        height_control='fdr',
                                        threshold=1.96,
                                        cluster_threshold=0,
                                        two_sided=False)
            print(f"Contrast maps are thresholded by:\n   {'fdr'} p<{0.05:.3f} threshold: {threshold:.3f}")
            fname_avg_FDR_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_Nruns_{runs_count}_avg_z_FDR_{str_contrast}_tSNR{tSNR_tresh}.nii"
            nib.save(img=all_avg_zmap_FDR, filename=fname_avg_FDR_zmap)
    return

def plot_beta_map(fmri_map=None, anat=None, contrastID=None):
    
    """_summary_
    Plot the fmri map as figure, no saving pictures
    """
    plotting.plot_stat_map(
            fmri_map,
            bg_img=anat,
            threshold=3,
            display_mode="ortho",
            black_bg=True,
            title=contrastID,
        )
    plotting.show()
    print('.')
    return

def func_fit_glm(func=None, tr=2, events=None, 
                motion=None, motion_par=None, 
                subjID='subj-DUMMY', smoothing_fwhm=8, 
                params=None, contrast='',
                tSNR_tresh=0,
                sesID='', runID=''):

    """_summary_
    Using all input variables to perform GL modelling and save all contrast map as Niftis
    """
    n_scans = func.shape[3]
    events['onset'] = events.loc[:, 'onset'] - tr / 2
    
    # make design matrix
    frame_times = tr * np.arange(n_scans)
    dm = make_first_level_design_matrix(
        frame_times=frame_times,
        hrf_model='spm',
        events=events,
        add_regs=motion,
        add_reg_names=motion_par,
        min_onset=0)

    # init model
    fmri_glm = FirstLevelModel(
        t_r=tr,
        subject_label=subjID,
        minimize_memory=True,
        high_pass=0.008,
        smoothing_fwhm=smoothing_fwhm)

    # fit model
    print("     # Fitting a GLM")
    fmri_maps = fmri_glm.fit(func,
                            design_matrices=dm)

    fdir_der_firstlvl = (params.get('fdir_proc')
                        + '/1st_level/'
                        + subjID + '/' + sesID + '/' + runID + '/'
                            )
    check_make_dir(fdir_der_firstlvl)

    print("     # Make nifti for each contrast")
    # Make nifti for each contrast
    for i, contrast_id in enumerate(contrast):
        loadingBar(i, len(contrast), contrast_id)
        z_map = fmri_maps.compute_contrast(contrast_id, output_type="z_score")
        str_contrast = contrast_id.replace(" ", "_")
        
        fname_zmap = fdir_der_firstlvl + f"{subjID}_{sesID}_{runID}_{str_contrast}_tSNR{tSNR_tresh}_zscore.nii"
        nib.save(z_map, fname_zmap)

    print("     # Done\n")

    return

def func_apply_glm_fixed_effect(participants=None, params=None, smoothing_fwhm=8, tSNR_tresh=30, contrast=None):
    
    """_summary_
    Iterate across subjects, sessions, runs. Calculate tSNR and only use run for GLM 
    if it passes given tSNR threshhold. Calculate all contrast maps and average across subjects.
    Save all contrast maps as nifti for each subject, session. 
    Average all contrast maps (niftis) of subjects, per session.
    """
    
    if participants is None:
        participants = []
    check_make_dir(params.get('fdir_proc_pre'))
        
    print(
        f'####\nPreprocessing of functional MR data is initiated for: {len(participants)} subjects\n####'
    )

    all_sessions = list()
    runs_count = 0
    for iSubj in range(len(participants)):
        print(f'# Identifying sessions from "derivatives/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_proc_pre') + '/' + participants[iSubj]

        if os.path.exists(tmp_fdir_subj):
            sessions = os.listdir(tmp_fdir_subj)
            sessions.sort()
            
            sessions = ['ses-1']
            
            for session in sessions:
                print(f'  # Performing calculation on session:   {session}')
                all_sessions.append(sessions)
                fdir_derivatives_func = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/func'
                                        )
                fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/anat'
                                        )
                if os.path.exists(fdir_derivatives_func) and os.path.exists(fdir_derivatives_anat) :
                    runs = os.listdir(fdir_derivatives_func)
                    runs = get_only_nifti(runs, 1)
                    
                    proc_lvl = ['_st_mcf_topUP_r_w.nii.gz']

                    runs[:] = [proc_stage for proc_stage in runs if any(item in proc_stage for item in proc_lvl)]
                    runs.sort()
                    fmri_runs       = []
                    design_matrices = []
                    if len(runs) > 0:
                        
                        for run in runs:
                            
                            # Grab information from filename
                            indices = run.split('_')
                            subjID = indices[0]
                            ses = indices[1]
                            runID = indices[2]
                            print(f'    # Performing calculation on run:   {runID}')
                            run_reduced = run.replace(proc_lvl[0], "")
                            
                            # Get functional data as nifti file
                            tmp_fname_final_funcMR = fdir_derivatives_func + '/' + run
                            if os.path.exists(tmp_fname_final_funcMR):
                                print(f'     Found functional data :   {tmp_fname_final_funcMR}.')
                                                                
                            # Perform image quality metric analysis
                            fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                + '/' + subjID 
                                + '/' + session + '/anat'
                                )
                            fname_parc = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_aparc+aseg_w.nii.gz'
    
                            # read the table containing the tSNR information of given 
                            # functional MR file
                            tSNR_dict = read_table_tSNR(tmp_fname_final_funcMR, params)
                            
                            if len(tSNR_dict) == 0:
                                # if no such file exists, calculate tSNR
                                tmp_tSNR_median = calc_tSNR_aligned(tmp_fname_final_funcMR, fname_parc, params)
                            else:
                                tmp_tSNR_median = float(tSNR_dict[0]['tSNR'])
                            
                            
                            if tmp_tSNR_median >= tSNR_tresh:
                                
                                # Count the runs that surpassed the tSNR threshold and pass 
                                # that information to the final session average filename 
                                runs_count = runs_count + 1
                            
                                # Get motion parameters for functional data as par file
                                tmp_fname_motion_funcMR = fdir_derivatives_func + '/' + run_reduced + '_st_mcf.par'
                                if os.path.exists(tmp_fname_motion_funcMR):
                                    print(f'     Found motion parameter data :   {tmp_fname_motion_funcMR}.')
                                
                                
                                # Get json file to functional data
                                tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]
                                tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'
                                tmp_fname_raw_funcMR_json = (tmp_fdir_func
                                                            + '/' + run_reduced + '.json')
                                if os.path.exists(tmp_fname_raw_funcMR_json):
                                    print(f'     Found .json file for functional data :   {tmp_fname_raw_funcMR_json}.')
                                with open(f'{tmp_fname_raw_funcMR_json}') as f:
                                    func_json = json.load(f)
                                RepetitionTime = func_json.get('RepetitionTime')
                                print(f'        Found TR of : {RepetitionTime}s for functional data.')
                                
                                # Get event file to functional data
                                eventspath = (tmp_fdir_func 
                                            + '/' + run_reduced
                                            + '.tsv')
                                if os.path.exists(eventspath):
                                    print(f'     Found evt data for functional data : {eventspath}.')
                                

                                print('     # Load data for processing functional data for stimuli...')
                                # Load data
                                func    = nib.load(tmp_fname_final_funcMR)
                                tr      = RepetitionTime
                                n_scans = func.shape[3]
                                
                                events  = pd.read_csv(eventspath,sep='\t')
                                
                                baseline = {'onset': [events.onset[47] + events.duration[47]],
                                            'duration': [( n_scans * tr ) - events.onset[47] + events.duration[47]],
                                            'trial_type': ['baseline']}
                                df_baseline = pd.DataFrame(baseline)
                                
                                events = pd.concat([events, df_baseline], ignore_index=True, sort=False)

    
                                motion_par = ["tx", "ty", "tz", "rx", "ry", "rz"]
                                motion = np.array( pd.read_csv(tmp_fname_motion_funcMR,
                                                    delim_whitespace=True,
                                                    names=motion_par) 
                                                )

                                
                                
                        
                                n_scans = func.shape[3]
                                events['onset'] = events.loc[:, 'onset'] - tr / 2
                                
                                # make design matrix
                                frame_times = tr * np.arange(n_scans)
                                dm = make_first_level_design_matrix(
                                    frame_times=frame_times,
                                    hrf_model='spm',
                                    events=events,
                                    add_regs=motion,
                                    add_reg_names=motion_par,
                                    min_onset=0)
                                
                                fmri_runs.append(func)
                                design_matrices.append(dm)
                                
                        if len(fmri_runs) > 0:
                            # init model
                            fmri_glm = FirstLevelModel(
                                t_r=tr,
                                minimize_memory=False,
                                signal_scaling=False,
                                high_pass=0.008,
                                smoothing_fwhm=smoothing_fwhm)

                            # fit model
                            print("     # Fitting a GLM")
                            fmri_maps = fmri_glm.fit(fmri_runs,
                                                    design_matrices=design_matrices)

                            fdir_der_firstlvl = (params.get('fdir_proc')
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/' + runID + '/'
                                                )
                            check_make_dir(fdir_der_firstlvl)
                            

                            print("     # Make nifti for each contrast")
                            # Make nifti for each contrast
                            for i, contrast_id in enumerate(contrast):
                                loadingBar(i, len(contrast), contrast_id)
                                z_map = fmri_maps.compute_contrast(contrast_id, output_type="z_score")
                                str_contrast = contrast_id.replace(" ", "_")
                                

                                fdir_der_firstlvl_final = (params.get('fdir_proc')
                                                            + '/1st_level/'
                                                            + subjID + '/' + ses + '/'
                                                            )
                                fname_avg_zmap = fdir_der_firstlvl_final + f"{subjID}_{ses}_FE_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                                nib.save(z_map, fname_avg_zmap)
                                
                            print("     # Done\n")
                    


    # Average all contrasts and runs of all sessions across subjects 
    all_sessions = list( np.unique(all_sessions))
    for session in all_sessions:
        print(session)
        fdir_group_firstlvl_final = ( params.get('fdir_proc') 
                                + '/1st_level/Group/' 
                                + session + '/'
                            )
        check_make_dir(fdir_group_firstlvl_final)
        for i, contrast_id in enumerate(contrast):
            loadingBar(i, len(contrast), contrast_id)
            str_contrast = contrast_id.replace(" ", "_")
            
            ALL_subj_image_array = list()
            subj_count = 0
            for j, subjID in enumerate(participants):
                tmp_fdir_firstlvl = ( params.get('fdir_proc') 
                                        + '/1st_level/'
                                        + subjID + '/' + ses + '/'
                                    )
                if os.path.exists(tmp_fdir_firstlvl):
                    fname_zmap2 = tmp_fdir_firstlvl + f"{subjID}_{ses}_FE_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                    if os.path.exists(fname_zmap2):
                        tmp_zmap2   = nib.load(fname_zmap2)
                        ALL_subj_image_array.append( tmp_zmap2.get_fdata() )
                        subj_count = subj_count +1
                    
            
            # Z scores              
            ALL_subj_image_array = np.array(ALL_subj_image_array)
            all_avg_zmap = np.mean(ALL_subj_image_array, axis=0)
            fname_avg_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_FE_avg_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                            
            all_avg_zmap_img = nib.Nifti1Image(all_avg_zmap, func.affine)
            nib.save(img=all_avg_zmap_img, filename=fname_avg_zmap)
            
            #all_sd_zmap = np.std(ALL_subj_image_array, axis=0)
            #fname_sd_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_Nruns_{runs_count}_sd_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                            
            #all_sd_zmap_img = nib.Nifti1Image(all_sd_zmap, func.affine)
            #nib.save(img=all_sd_zmap_img, filename=fname_sd_zmap)
            
            vol_data_init   = nib.load(fname_avg_zmap)
            all_avg_zmap_FDR, threshold = threshold_stats_img(stat_img=vol_data_init, alpha=0.05,
                                        height_control='fdr',
                                        threshold=1.96,
                                        cluster_threshold=0,
                                        two_sided=False)
            print(f"Contrast maps are thresholded by:\n   {'fdr'} p<{0.05:.3f} threshold: {threshold:.3f}")
            fname_avg_FDR_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_FE_avg_z_FDR_{str_contrast}_tSNR{tSNR_tresh}.nii"
            nib.save(img=all_avg_zmap_FDR, filename=fname_avg_FDR_zmap)
            
    return

def func_apply_glm_groupaverage_BOLD(participants=None, params=None, smoothing_fwhm=8, tSNR_tresh=30, contrast=None):
    
    """_summary_
    Iterate across subjects, sessions, runs. Calculate tSNR and only use run for GLM 
    if it passes given tSNR threshhold. Calculate all contrast maps and average across subjects.
    Save all contrast maps as nifti for each subject, session. 
    Average all contrast maps (niftis) of subjects, per session.
    """
    
    if participants is None:
        participants = []
    check_make_dir(params.get('fdir_proc_pre'))
    
    participants = list( np.unique(participants))
    print(
        f'| TITLE | #### Constructing average BOLD response per cond. per subj. of functional MR data ####'
    )
    sessions = ['ses-1', 'ses-2', 'ses-3', 'ses-4', 'ses-5', 'ses-6', 'ses-7', 'ses-8', 'ses-9', 'ses-10', 'ses-11', 'ses-12' ]
    # Make the masks based on group-averaged contrast results
    print(f'| INFO | ~~~ Creating the masks based on group-averaged contrast results ~~~')
    contrast_masks = {}
    for session in sessions:
        for i, contrast_id in enumerate(contrast):
            loadingBar(i, len(contrast), contrast_id)
            str_contrast = contrast_id.replace(" ", "_")

            fdir_group_firstlvl_final = ( params.get('fdir_proc') 
                                        + '/1st_level/Group/' 
                                        + session + '/'
                                        )
            if os.path.isdir(fdir_group_firstlvl_final):
                fname_avg_FDR_zmap = (fdir_group_firstlvl_final 
                                    + f"GroupN38_{session}_Nruns_104_FE_avg_z_FDR_{str_contrast}_tSNR{tSNR_tresh}.nii")
                Groupavg_FDR_zmap  = nib.load(fname_avg_FDR_zmap)
                
                Mask_Groupavg_FDR_zmap = binarize_img(Groupavg_FDR_zmap, threshold=0.1)

                contrast_masks[session, i,1] = Mask_Groupavg_FDR_zmap
            else:
                print('| INFO | Skip session...' )
    print(f'| INFO | -> DONE ~~~ Creating the masks based on group-averaged contrast results ~~~')


    all_sessions = list()
    runs_count = 0
    for iSubj in range(len(participants)):
        print(f'| INFO | # Identifying sessions from "derivatives/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_proc_pre') + '/' + participants[iSubj]

        if os.path.exists(tmp_fdir_subj):
            sessions = os.listdir(tmp_fdir_subj)
            sessions.sort()
            for session in sessions:
                print(f'  # Performing calculation on session:   {session}')
                all_sessions.append(sessions)
                fdir_derivatives_func = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/func'
                                        )
                fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                        + '/' + participants[iSubj] 
                                        + '/' + session + '/anat'
                                        )
                if os.path.exists(fdir_derivatives_func) and os.path.exists(fdir_derivatives_anat) :
                    runs = os.listdir(fdir_derivatives_func)
                    runs = get_only_nifti(runs, 1)
                    
                    proc_lvl = ['_st_mcf_topUP_r_w.nii.gz']

                    runs[:] = [proc_stage for proc_stage in runs if any(item in proc_stage for item in proc_lvl)]
                    runs.sort()
                    fmri_runs       = []
                    design_matrices = []
                    if len(runs) > 0:
                        
                        all_events = pd.DataFrame()
                        all_runID = []
                        for run in runs:
                            
                            # Grab information from filename
                            indices = run.split('_')
                            subjID = indices[0]
                            ses = indices[1]
                            runID = indices[2]
                            print(f'    # Performing calculation on run:   {runID}')
                            run_reduced = run.replace(proc_lvl[0], "")
                            
                            # Get functional data as nifti file
                            tmp_fname_final_funcMR = fdir_derivatives_func + '/' + run
                            if os.path.exists(tmp_fname_final_funcMR):
                                print(f'     Found functional data :   {tmp_fname_final_funcMR}.')
                                                                
                            # Perform image quality metric analysis
                            fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                + '/' + subjID 
                                + '/' + session + '/anat'
                                )
                            fname_parc = f'{fdir_derivatives_anat}/{subjID}_{session}_T1w_aparc+aseg_w.nii.gz'
    
                            # read the table containing the tSNR information of given 
                            # functional MR file
                            tSNR_dict = read_table_tSNR(tmp_fname_final_funcMR, params)
                            
                            if len(tSNR_dict) == 0:
                                # if no such file exists, calculate tSNR
                                tmp_tSNR_median = calc_tSNR_aligned(tmp_fname_final_funcMR, fname_parc, params)
                            else:
                                tmp_tSNR_median = float(tSNR_dict[0]['tSNR'])
                            
                            
                            if tmp_tSNR_median >= tSNR_tresh:
                                
                                # Count the runs that surpassed the tSNR threshold and pass 
                                # that information to the final session average filename 
                                runs_count = runs_count + 1
                            
                                # Get motion parameters for functional data as par file
                                tmp_fname_motion_funcMR = fdir_derivatives_func + '/' + run_reduced + '_st_mcf.par'
                                if os.path.exists(tmp_fname_motion_funcMR):
                                    print(f'     Found motion parameter data :   {tmp_fname_motion_funcMR}.')
                                
                                
                                # Get json file to functional data
                                tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]
                                tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'
                                tmp_fname_raw_funcMR_json = (tmp_fdir_func
                                                            + '/' + run_reduced + '.json')
                                if os.path.exists(tmp_fname_raw_funcMR_json):
                                    print(f'     Found .json file for functional data :   {tmp_fname_raw_funcMR_json}.')
                                with open(f'{tmp_fname_raw_funcMR_json}') as f:
                                    func_json = json.load(f)
                                RepetitionTime = func_json.get('RepetitionTime')
                                print(f'        Found TR of : {RepetitionTime}s for functional data.')
                                
                                # Get event file to functional data
                                eventspath = (tmp_fdir_func 
                                            + '/' + run_reduced
                                            + '.tsv')
                                if os.path.exists(eventspath):
                                    print(f'     Found evt data for functional data : {eventspath}.')
                                

                                print('     # Load data for processing functional data for stimuli...')
                                # Load data
                                func    = nib.load(tmp_fname_final_funcMR)
                                tr      = RepetitionTime
                                n_scans = func.shape[3]
                                
                                events  = pd.read_csv(eventspath,sep='\t')
                                
                                baseline = {'onset': [events.onset[47] + events.duration[47]],
                                            'duration': [( n_scans * tr ) - events.onset[47] + events.duration[47]],
                                            'trial_type': ['baseline']}
                                df_baseline = pd.DataFrame(baseline)
                                
                                events = pd.concat([events, df_baseline], ignore_index=True, sort=False)
                                

    
                                motion_par = ["tx", "ty", "tz", "rx", "ry", "rz"]
                                motion = np.array( pd.read_csv(tmp_fname_motion_funcMR,
                                                    delim_whitespace=True,
                                                    names=motion_par) 
                                                )

                                n_scans = func.shape[3]
                                
                                # make design matrix
                                frame_times = tr * np.arange(n_scans)
                                dm = make_first_level_design_matrix(
                                    frame_times=frame_times,
                                    hrf_model='spm',
                                    events=events,
                                    add_regs=motion,
                                    add_reg_names=motion_par,
                                    min_onset=0)
                                
                                events['run'] = runID
                                all_events = pd.concat([all_events,events], ignore_index=True, sort=False)
                                fmri_runs.append(func)
                                design_matrices.append(dm)
                                
                        if len(fmri_runs) > 0:
                            # init model
                            fmri_glm = FirstLevelModel(
                                t_r=tr,
                                minimize_memory=False,
                                signal_scaling=(0, 1),    # 0,1,(0, 1)
                                smoothing_fwhm=smoothing_fwhm)

                            # fit model
                            print("     # Fitting a GLM")
                            fmri_maps = fmri_glm.fit(fmri_runs,
                                                    design_matrices=design_matrices)

                            fdir_der_firstlvl = (params.get('fdir_proc')
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/' + runID + '/'
                                                )
                            check_make_dir(fdir_der_firstlvl)
                            

                            print("     # Make subject-average BOLD response for each trial type")
                            
                            pre_stim  = -1
                            post_stim = +10
                            # Make nifti for each contrast
                            ident_runs = list (np.unique(all_events['run']))
                            for i, contrast_id in enumerate(contrast):
                                #loadingBar(i, len(contrast), contrast_id)
                                #z_map = fmri_maps.compute_contrast(contrast_id, output_type="z_score")
                                str_contrast = contrast_id.replace(" ", "_")

                                tmp_trial_type = str_contrast.removesuffix('_-_null_event')
                                
                                got_session = False
                                for i_mask in contrast_masks:
                                    if i_mask[0] == session:
                                        got_session = True
                                    
                                N_evts_trial_type = sum ( all_events['trial_type'].str.contains(tmp_trial_type))
                                if (N_evts_trial_type > 0) & ( got_session == True ):
                                    
                                    tmp_contrast_masks = contrast_masks[session, i,1]
                                                                        
                                    ROI_masker = NiftiMasker(tmp_contrast_masks,
                                                            standardize='psc', t_r=tr,
                                                            detrend=True)
                                    ROI_nifti = ROI_masker.mask_img
                                    
                                    if sum ( sum( sum( ROI_nifti.get_fdata() ) )) > 0:
                                        tmp_BOLD_tc  = list()
                                        tmp_BOLD_tc2 = list()
                                        for irun in np.arange(len(fmri_runs)):
                                            tmp_observed_timeseries = ROI_masker.fit_transform(fmri_runs[irun])
                                            tmp_observed_timeseries = np.mean( tmp_observed_timeseries, 1)
                                            
                                            tmp_predicted_timeseries = ROI_masker.fit_transform(fmri_maps.predicted[irun])
                                            tmp_predicted_timeseries = np.mean( tmp_predicted_timeseries, 1)
                                                                                        
                                            tmp_evt_cond = (all_events['run'] == ident_runs[irun])
                                            tmp_con_cond = (all_events['trial_type'] == tmp_trial_type)
                                            
                                            tmp_events   = all_events[tmp_evt_cond].onset.values
                                            tmp_con_cond = tmp_con_cond[tmp_evt_cond]
                                            tmp_events   = ( tmp_events / tr )
                                            
                                            tmp_contrast_idx = np.floor( tmp_events[tmp_con_cond] )
                                            
                                            
                                            for j, evt_idx in enumerate(tmp_contrast_idx):
                                                evt_idx = int(evt_idx)
                                                
                                                contrast_idx_start = evt_idx + pre_stim
                                                contrast_idx_end   = evt_idx + post_stim
                                                
                                                if contrast_idx_end > len(tmp_predicted_timeseries):
                                                    contrast_idx_end = len(tmp_predicted_timeseries)
                                                
                                                
                                                tmp_BOLD_tc_loop  = tmp_predicted_timeseries[ contrast_idx_start:contrast_idx_end]

                                                tmp_BOLD_tc_loop2 = tmp_observed_timeseries[ contrast_idx_start:contrast_idx_end]
                                                
                                                tmp_BOLD_tc.append(tmp_BOLD_tc_loop)
                                                tmp_BOLD_tc2.append(tmp_BOLD_tc_loop2)
                                        
                                        tmp_BOLD_tc = tmp_BOLD_tc[1::]   
                                        tmp_BOLD_tc2 = tmp_BOLD_tc2[1::]
                                        tmp_BOLD_tc  = np.array(tmp_BOLD_tc)
                                        tmp_BOLD_tc2 = np.array(tmp_BOLD_tc2)
                                        
                                        
                                        avg_BOLD_tc = np.median(tmp_BOLD_tc,0)
                                        avg_BOLD_tc2 = np.median(tmp_BOLD_tc2,0)
                                                                                
                                                                                
                                        fdir_der_firstlvl = (params.get('fdir_proc')
                                                    + '/1st_level/'
                                                    + subjID + '/' + ses + '/')
                                                    
                                        fname_avg_BOLD_tc = ( fdir_der_firstlvl
                                                    + f"{subjID}_{ses}_avg_BOLD_{tmp_trial_type}_tSNR{tSNR_tresh}" )
                                        np.save(fname_avg_BOLD_tc, avg_BOLD_tc,
                                                allow_pickle=True, fix_imports=True)
                            print("     # Done\n")

    # Average all contrasts and runs of all sessions across subjects 
    all_sessions = sessions
    for session in all_sessions:
        print(session)
        fdir_group_firstlvl_final = ( params.get('fdir_proc') 
                                + '/1st_level/Group/' 
                                + session + '/'
                            )
        
        
        BOLD_fMRI_dataframe = pd.DataFrame(columns=['subject', 'timepoint','signal', 'condition', 'session'])
        
        check_make_dir(fdir_group_firstlvl_final)
        for i, contrast_id in enumerate(contrast):
            loadingBar(i, len(contrast), contrast_id)
            str_contrast = contrast_id.replace(" ", "_")

            tmp_trial_type = str_contrast.removesuffix('_-_null_event')
                                
            ALL_tmp_avgBOLD_subj = list()
            subj_count = 0
            for j, subjID in enumerate(participants):
                tmp_fdir_firstlvl = ( params.get('fdir_proc') 
                                        + '/1st_level/'
                                        + subjID + '/' + ses + '/'
                                    )
                if os.path.exists(tmp_fdir_firstlvl):
                    fname_avgBOLD_subj = tmp_fdir_firstlvl = ( tmp_fdir_firstlvl
                                                + f"{subjID}_{ses}_avg_BOLD_{tmp_trial_type}_tSNR{tSNR_tresh}.npy" )
                    print(fname_avgBOLD_subj)
                    if os.path.exists(fname_avgBOLD_subj):
                        tmp_avgBOLD_subj   = np.load(fname_avgBOLD_subj)

                        ALL_tmp_avgBOLD_subj.append( tmp_avgBOLD_subj )
                        subj_count = subj_count +1
                        
            if subj_count > 0:
                # Group Average BOLD response         
                ALL_avgBOLD_subj = np.array(ALL_tmp_avgBOLD_subj)
                
                all_avg_BOLD = np.mean(ALL_avgBOLD_subj, axis=0)
                #all_std_BOLD = np.std(ALL_avgBOLD_subj, axis=0)
                
                timepoints = np.arange(-1, all_avg_BOLD.shape[0]-1)
                
                idx_counter = 0
                tmp_dataframe = pd.DataFrame(columns=['subject', 'timepoint','signal', 'condition', 'session'])
                for kk, tmp_data_BOLD in enumerate(ALL_avgBOLD_subj):
                    for k, i_tp in enumerate(timepoints):
                        tmp_dataframe.loc[idx_counter] = {'subject':str(kk), 'timepoint':i_tp, 
                                        'signal': tmp_data_BOLD[k], 'condition': tmp_trial_type,
                                        'session':session}
                        
                        idx_counter = idx_counter + 1 
                

                BOLD_fMRI_dataframe = pd.concat([BOLD_fMRI_dataframe, tmp_dataframe], ignore_index=True)

                
                fname_all_avg_BOLD = ( fdir_group_firstlvl_final +
                    f"GroupN{subj_count}_{session}_Nruns_{runs_count}_FE_avg_BOLD_{tmp_trial_type}_tSNR{tSNR_tresh}")
                                    
                np.save(fname_all_avg_BOLD, BOLD_fMRI_dataframe)
                
            # TODO: Add relation plots to visualize session depended changes (https://seaborn.pydata.org/examples/timeseries_facets.html)
            plt.figure
            sns.lineplot(x="timepoint", y="signal",
                hue="condition",data=BOLD_fMRI_dataframe)
            plt.xlim((-2, 10))
            plt.title(f'BOLD response (N:{subj_count})')
            plt.ylabel('BOLD signal')
            plt.xlabel('Time (scans)')
            fname_fig_avg_BOLD = ( fdir_group_firstlvl_final +
                    f"GroupN{subj_count}_{session}_Nruns_{runs_count}_FE_avg_BOLD_{tmp_trial_type}_tSNR{tSNR_tresh}")
                
            plt.savefig(fname_fig_avg_BOLD)
            
    return
