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
from nilearn import glm, image, masking, plotting
from nilearn.glm.first_level import (FirstLevelModel,
                                    make_first_level_design_matrix)
from nilearn.reporting import make_glm_report


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
                            runs_count = runs_count + 1
                            
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
    
                            
                            tSNR_dict = read_table_tSNR(tmp_fname_final_funcMR, params)
                            
                            if len(tSNR_dict) == 0:
                                tmp_tSNR_median = calc_tSNR_aligned(tmp_fname_final_funcMR, fname_parc, params)
                            else:
                                tmp_tSNR_median = float(tSNR_dict[0]['tSNR'])
                                
                            if tmp_tSNR_median >= tSNR_tresh:
                                
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
                            image_array2 = list()
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
                                avg_zmap = np.mean(image_array, axis=0)
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
            ALL_subj_image_array2 = list()
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
            
            all_sd_zmap = np.std(ALL_subj_image_array, axis=0)
            fname_sd_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_Nruns_{runs_count}_sd_z_{str_contrast}_tSNR{tSNR_tresh}.nii"
                            
            all_sd_zmap_img = nib.Nifti1Image(all_sd_zmap, func.affine)
            nib.save(img=all_sd_zmap_img, filename=fname_sd_zmap)

    return

def func_apply_glm(participants=None, params=None, smoothing_fwhm=8, contrast=None):
    

    if participants is None:
        participants = []
    check_make_dir(params.get('fdir_proc_pre'))
        
    print(
        f'####\nPreprocessing of functional MR data is initiated for: {len(participants)} subjects\n####'
    )

    all_sessions = list()
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
                            
                            # Get anatomical data as nifti file
                            tmp_fname_MR = fdir_derivatives_anat + f'/{subjID}_{ses}_T1w_nu_w.nii.gz'
                            if os.path.exists(tmp_fname_MR):
                                print(f'     Found anatomical data :   {tmp_fname_MR}.')
                                    
                            # Get motion parameters for functional data as par file
                            tmp_fname_motion_funcMR = fdir_derivatives_func + '/' + run_reduced + '_st_mcf.par'
                            
                            if os.path.exists(tmp_fname_motion_funcMR):
                                print(f'     Found motion parameter data :   {tmp_fname_motion_funcMR}.')
                                
                            # Get json file to functional data
                            tmp_fdir_subj = params.get(
                                'fdir_data') + participants[iSubj]
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



                            # Load data
                            print('     # Load data for processing functional data for stimuli...')
                            func = nib.load(tmp_fname_final_funcMR)
                            tr = RepetitionTime
                            n_scans = func.shape[3]

                            events = pd.read_csv(eventspath, sep='\t')

                            baseline = {'onset': [events.onset[47] + events.duration[47]],
                                        'duration': [(n_scans * tr) - events.onset[47] + events.duration[47]],
                                        'trial_type': ['baseline']}
                            df_baseline = pd.DataFrame(baseline)

                            events = pd.concat(
                                [events, df_baseline], ignore_index=True, sort=False)

                            motion_par = ["tx", "ty", "tz",
                                        "rx", "ry", "rz"]
                            motion = np.array(pd.read_csv(tmp_fname_motion_funcMR,
                                                        delim_whitespace=True,
                                                        names=motion_par)
                                            )

                            # make design matrix
                            frame_times = tr * np.arange(n_scans)
                            dm = make_first_level_design_matrix(
                                frame_times=frame_times,
                                hrf_model='spm + derivative',
                                events=events,
                                add_regs=motion,
                                add_reg_names=motion_par,
                                min_onset=0)

                            # init model
                            # TODO drift model wanted?
                            fmri_glm = FirstLevelModel(
                                t_r=tr,
                                subject_label=participants[iSubj],
                                minimize_memory=True,
                                high_pass=0.008,
                                smoothing_fwhm=smoothing_fwhm)

                            # fit model
                            print("     # Fitting a GLM")
                            fmri_maps = fmri_glm.fit(func,
                                                    design_matrices=dm)

                            fdir_der_firstlvl = (params.get('fdir_proc')
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/' + runID + '/'
                                                )
                            check_make_dir(fdir_der_firstlvl)

                            print("     # Make nifti for each contrast")
                            # Make nifti for each contrast
                            for i, contrast_id in enumerate(contrast):
                                loadingBar(i, len(contrast), contrast_id)
                                z_map = fmri_maps.compute_contrast(
                                    contrast_id, output_type="z_score")
                                str_contrast = contrast_id.replace(" - ", "_vs_")
                                
                                fname_zmap = fdir_der_firstlvl + f"{subjID}_{ses}_{runID}_{str_contrast}_zscore.nii"
                                nib.save(z_map, fname_zmap)

                                #p_map = fmri_maps.compute_contrast(
                                #        contrast_id, output_type="p_value")
                                #fname_pmap = fdir_der_firstlvl + f"{subjID}_{ses}_{runID}_{str_contrast}_pvalue.nii"
                                #nib.save(p_map, fname_pmap)

                            print("     # Done\n")

                        # Average all contrasts and runs
                        fdir_der_firstlvl_final = (params.get('fdir_proc')
                                                + '/1st_level/'
                                                + subjID + '/' + ses + '/'
                                                )
                        check_make_dir(fdir_der_firstlvl_final)
                        print("Average for each contrast over all runs")
                        for i, contrast_id in enumerate(contrast):
                            loadingBar(i, len(contrast), contrast_id)
                            str_contrast = contrast_id.replace(" - ", "_vs_")

                            image_array = list()
                            image_array2 = list()
                            for j, run in enumerate(runs):
                                # Identify identifier
                                indices = run.split('_')
                                subjID = indices[0]
                                ses = indices[1]
                                runID = indices[2]

                                fdir_der_firstlvl = (params.get('fdir_proc')
                                                    + '/1st_level/'
                                                    + subjID + '/' + ses + '/' + runID + '/'
                                                    )
                                fname_zmap = fdir_der_firstlvl + f"{subjID}_{ses}_{runID}_{str_contrast}_zscore.nii"
                                tmp_zmap = nib.load(fname_zmap)
                                image_array.append(tmp_zmap.get_fdata())

                                #fname_pmap = fdir_der_firstlvl + f"{subjID}_{ses}_{runID}_{str_contrast}_pvalue.nii"
                                #tmp_pmap = nib.load(fname_pmap)
                                #image_array2.append(tmp_pmap.get_fdata())

                            if len(image_array) > 0:
                                image_array = np.array(image_array)
                                avg_zmap = np.mean(image_array, axis=0)
                                fname_avg_zmap = fdir_der_firstlvl_final + f"{subjID}_{ses}_avg_z_{str_contrast}.nii"

                                avg_zmap_img = nib.Nifti1Image(avg_zmap, func.affine)
                                nib.save(img=avg_zmap_img,filename=fname_avg_zmap)

                            #if len(image_array2) > 0:
                            #    image_array2 = np.array(image_array2)
                            #    avg_pmap = np.mean(image_array2, axis=0)
                            #    fname_avg_pmap = fdir_der_firstlvl_final + f"{subjID}_{ses}_avg_p_{str_contrast}.nii"

                            #    avg_pmap_img = nib.Nifti1Image(avg_pmap, func.affine)
                            #    nib.save(img=avg_pmap_img,filename=fname_avg_pmap)


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
            str_contrast = contrast_id.replace(" - ", "_vs_")
            
            ALL_subj_image_array = list()
            ALL_subj_image_array2 = list()
            subj_count = 0
            for j, subjID in enumerate(participants):
                tmp_fdir_firstlvl = ( params.get('fdir_proc') 
                                        + '/1st_level/'
                                        + subjID + '/' + ses + '/'
                                    )
                if os.path.exists(tmp_fdir_firstlvl):
                    fname_zmap2 = tmp_fdir_firstlvl + f"{subjID}_{ses}_avg_z_{str_contrast}.nii"
                    if os.path.exists(fname_zmap2):
                        tmp_zmap2   = nib.load(fname_zmap2)
                        ALL_subj_image_array.append( tmp_zmap2.get_fdata() )
                        subj_count = subj_count +1
                    
                    #fname_pmap2 = tmp_fdir_firstlvl + f"{subjID}_{ses}_avg_p_{str_contrast}.nii"
                    #if os.path.exists(fname_pmap2):
                    #    tmp_pmap2   = nib.load(fname_pmap2)
                    #    ALL_subj_image_array2.append( tmp_pmap2.get_fdata() )
                    #    subj_count = subj_count +1
            
            # Z scores              
            ALL_subj_image_array = np.array(ALL_subj_image_array)
            all_avg_zmap = np.mean(ALL_subj_image_array, axis=0)
            fname_avg_zmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_avg_z_{str_contrast}.nii"
                            
            all_avg_zmap_img = nib.Nifti1Image(all_avg_zmap, func.affine)
            nib.save(img=all_avg_zmap_img, filename=fname_avg_zmap)
            
            ## P values   
            #ALL_subj_image_array2 = np.array(ALL_subj_image_array2)
            #all_avg_pmap = np.mean(ALL_subj_image_array2, axis=0)
            #fname_avg_pmap = fdir_group_firstlvl_final + f"GroupN{subj_count}_{ses}_avg_p_{str_contrast}.nii"
                            
            #all_avg_pmap_img = nib.Nifti1Image(all_avg_pmap, func.affine)
            #nib.save(img=all_avg_pmap_img, filename=fname_avg_pmap)


                            

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
    
    # make design matrix
    frame_times = tr * np.arange(n_scans)
    dm = make_first_level_design_matrix(
        frame_times=frame_times,
        hrf_model='spm + derivative',
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
        slice_time_ref=0.5,
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
                                
