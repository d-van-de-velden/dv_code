"""Main analysis script."""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
###############################################################
#  Import section
import os

from dv_code.scripts.dcm2BIDS_conversion import (convert_beh2BIDS,
                                                convert_dcm2BIDS,
                                                convert_evt2BIDS)
from dv_code.scripts.download_raw import (download_datashare_behavioral,
                                        download_datashare_dcm)
from dv_code.scripts.group.analyze_group import overview_tSNR
from dv_code.scripts.initialize_study import initialize_study
from dv_code.scripts.misc.read_analysis_params import read_analysis_params
from dv_code.scripts.preprocessing.func_preproc import do_func_preproc
from dv_code.scripts.preprocessing.struc_preproc import do_FS_recon
from dv_code.scripts.processing.func_1stlvl_analysis import func_apply_glm_groupaverage_BOLD, func_apply_glm_fixed_effect
from dv_code.scripts.work_participants import get_participants
from dv_code.scripts.processing.func_2ndlvl_analysis import make_overlap_contrastmap_surface



#  Establish study (-> RunIt=True), already done it? (->RunIt=False)
initialize_study(RunIt=False)

#  Load study framework parameter
params = read_analysis_params(f'{os.getcwd()}/analysis_params.json')

fdir_analysis = params.get('fdir_analysis')
fdir_datashare = params.get('fdir_url_datashare')
fdir_sourcedata = params.get('fdir_sourcedata')
participants = get_participants(params)


#  Grab data (DICOM and/or logfiles) from external cloud bases (e.g.: DataShare)
download_datashare_dcm(params)
download_datashare_behavioral(params)

# Grab all the subjects in the folder and create a list and add subject info
convert_dcm2BIDS(params)
convert_evt2BIDS(params)
convert_beh2BIDS(params)

# Run first pre-processing (structural)
participants = get_participants(params)
do_FS_recon(fdir_analysis, participants, use_HPC = 0, params=params)


# Run first pre-processing (functional)
participants  = get_participants(params)
do_func_preproc(participants, params, forceRun=False, use_HPC=False)


# Make group quality overview (functional)
participants = get_participants(params)
overview_tSNR(participants, params, tSNR_threshold=40, N_runsOK=2)


# Make 1st level analysis (functional)
contrast = [
        'unimodal_images - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'unimodal_audios - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'unimodal_audios - unimodal_images',                                # C1_vis; Unimodality specific activation; p(FDR)<0.001
        'unimodal_images - unimodal_audios',                                # C1_aud; Unimodality specific activation; p(FDR)<0.001
        '(congruent > unimodal_audios) & (congruent > unimodal_images)',    # C2
        'congruent - incongruent',                                          # Congruency contrast map; p(FDR) < 0.05
        'congruent - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'incongruent - null_event', 
]
participants = get_participants(params)
for itSNR in [40]:
        func_apply_glm_fixed_effect(participants, params, smoothing_fwhm=8, tSNR_tresh=itSNR, contrast=contrast)

# Make 1st level analysis (functional)
contrast = [
        'unimodal_images - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'unimodal_audios - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'congruent - null_event',                                     # Unimodality activation; q(FDR) < 0.0001
        'incongruent - null_event', 
]
participants = get_participants(params)
for itSNR in [40]: 
        func_apply_glm_groupaverage_BOLD(participants, params, smoothing_fwhm=8, tSNR_tresh=itSNR, contrast=contrast)



#%%

# 1. Plot unimodal contrasts
#    a) auditory activation --> [unimodal auditory > unimodal visual, q(FDR) < 0.001]

fdir = "/data/p_02825/cocoa/data/derivatives/1st_level/Group/ses-1/"
fname_in       = [fdir + "GroupN21_ses-1_FE_avg_z_FDR_unimodal_audios_-_null_event_tSNR40.nii",]
plotting_title = 'auditory activation'
fname_out      = '/data/p_02825/cocoa/data/derivatives/2nd_level/Group/ses-1/Average_FDRcorrected_UniAudio.png'
make_overlap_contrastmap_surface(fname_in, plotting_title,
                        alpha=0.05, mcc='fpr', cluster_ext=2, 
                        fname_out=fname_out)

#    b) visual activation -- > [unimodal visual > unimodal auditory, q(FDR) < 0.001]
fname_in       = [fdir + "GroupN21_ses-1_FE_avg_z_unimodal_images_-_null_event_tSNR40.nii"]
plotting_title = 'visual activation'
fname_out      = '/data/p_02825/cocoa/data/derivatives/2nd_level/Group/ses-1/Average_UniVisual.png'
make_overlap_contrastmap_surface(fname_in, plotting_title,
                        alpha=0.1, mcc='fpr', cluster_ext=2, 
                        fname_out=fname_out)

#    c) congruent activation -- > 
fname_in       = [fdir + "GroupN21_ses-1_FE_avg_z_FDR_congruent_-_null_event_tSNR40.nii"]
plotting_title = 'congruent activation'
fname_out      = '/data/p_02825/cocoa/data/derivatives/2nd_level/Group/ses-1/Avergae_FDRcorrected_BimodalCongruent.png'
make_overlap_contrastmap_surface(fname_in, plotting_title,
                        alpha=0.05, mcc='fpr', cluster_ext=2, 
                        fname_out=fname_out)

#    d) incongruent activation -- > 
fname_in       = [fdir + "GroupN21_ses-1_FE_avg_z_FDR_incongruent_-_null_event_tSNR40.nii"]
plotting_title = 'incongruent activation'
fname_out      = '/data/p_02825/cocoa/data/derivatives/2nd_level/Group/ses-1/Average_FDRcorrected_BimodalINcongruent.png'
make_overlap_contrastmap_surface(fname_in, plotting_title,
                        alpha=0.05, mcc='fpr', cluster_ext=2, 
                        fname_out=fname_out)

# 2. Make C2 contrast --> Overlap('congruent - unimodal_audios' AND 'congruent - unimodal_images')
fname_in       = [fdir + "GroupN21_ses-1_FE_avg_z_FDR_(congruent_>_unimodal_audios)_&_(congruent_>_unimodal_images)_tSNR40.nii"]
plotting_title = 'C2 contrast'
fname_out      = '/data/p_02825/cocoa/data/derivatives/2nd_level/Group/ses-1/Average_FDRcorrected_C2.png'
make_overlap_contrastmap_surface(fname_in, plotting_title,
                        alpha=0.05, mcc='fpr', cluster_ext=2, 
                        fname_out=fname_out)


