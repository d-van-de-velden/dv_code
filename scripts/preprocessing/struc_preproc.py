#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os
import stat
import numpy as np
import nibabel as nib
from subprocess import call

from dv_code.scripts.misc.feedback_messages import error_message
from dv_code.scripts.misc import check_make_dir


#fdir_analysis= fdir_analysis
#participants = participants
#use_HPC = 0

def do_FS_recon(fdir_analysis="", participants=None, use_HPC = 0, params=None):
    if participants is None:
        participants = []

    check_make_dir(params.get('fdir_fs'))

    if not fdir_analysis:
        input_str = "Variable 1 (fdir_analysis) is not provided. ABORT"
        error_message(input_str)
    
    env_path = os.getenv('PATH')
    host = '.'
    if 'cbs.mpg' in env_path:
        host = 'mpg'

    print(f'FreeSurfer recon_all is initiated for: {len(participants)} subjects')

    for iSubj in range(len(participants)):
        print(f'Identifying sessions from "data/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

        sessions = os.listdir(tmp_fdir_subj)
        for session in sessions:

            tmp_fname_T1w = (params.get('fdir_data') + participants[iSubj] + '/' + session 
                            + '/anat/' + participants[iSubj] + '_' + session + '_T1w.nii.gz')
            if os.path.exists(tmp_fname_T1w) == True:
                
                
                fdir_FS_session = params.get('fdir_fs') + '/' + session
                check_make_dir(fdir_FS_session)
                
                fname_FS_done = (fdir_FS_session + f'/{participants[iSubj]}/scripts/recon-all.done')
                if os.path.exists(fname_FS_done) == False:
                    print(tmp_fname_T1w)

                    fdir_bash_script = (params.get('fdir_bash') 
                                        + '/' + participants[iSubj]
                                        + '/' + session
                                        + '/FS_reconall/'
                                        )
                    check_make_dir(fdir_bash_script)
                    fname_bash_script = f'{fdir_bash_script}/run_fs_recon_all.sh'
                    fname_bash_script_ = f'{fdir_bash_script}/runner_fs_recon_all.sh'
                    check_make_dir(fdir_bash_script)

                    with open(fname_bash_script, 'w') as rsh:
                            rsh.write('''
                            #! /bin/sh
                            echo "I start FREESURFER now...."
                            FREESURFER bash ''' + fname_bash_script_
                            )

                    with open(fname_bash_script_, 'w') as rsh:
                            rsh.write('''\
                            #! /bin/bash
                            echo "I configure FREESURFER now...."
                            export SUBJECTS_DIR=''' + fdir_FS_session + ''' 
                            cd $SUBJECTS_DIR
                            recon-all -subjid ''' + participants[iSubj] + ''' -i ''' + tmp_fname_T1w + ''' -all
                            exit 0
                            ''')

                    os.chmod(fname_bash_script, stat.S_IRWXU)
                    os.chmod(fname_bash_script_, stat.S_IRWXU)
                    if use_HPC == 1:
                        print('Use HPC cluster for computation...')
                    else:
                        print('Use local machine for computation...')
                        fname_block_file = f'{fdir_bash_script}/blocked.txt'

                        if os.path.exists(fname_block_file) == True:
                            print('Job already running....')
                        else:
                            with open(fname_block_file, 'w') as rsh:
                                rsh.write("BLOCKED JOB\nfname= " + fname_bash_script)

                            if host == 'mpg':
                                os.system(f"xterm -e bash {fname_bash_script}")
                            else:
                                os.system(f"xterm -e bash {fname_bash_script_}")
                            os.remove(fname_block_file)

                else:
                    print(f'FreeSurfer recon-all already done for : {participants[iSubj]}')
                    
                #
                # TODO: Maybe add Afni Suma
                #'''@SUMA_Make_Spec_FS                                            \
                #-fs_setup                                                 \
                #-NIFTI                                                    \
                #-sid    ${subj}                                           \
                #-fspath ${dir_fs}/${subj}                                 \
                #|& tee  ${dir_echo}/o.01_suma_makespec_${subj}.txt
                #'''

                print('# INFO # Bringing anatomical derivative files to folder..')
                fdir_derivatives = (params.get('fdir_proc_pre')
                                    + '/' + participants[iSubj]
                                    + '/' + session + '/anat/'
                                    )
                check_make_dir(fdir_derivatives)
                # Source
                fname_T1w_nu0 = (fdir_FS_session + '/'
                                + participants[iSubj] + '/mri/nu.mgz')
                fname_T1w_parc0 = (fdir_FS_session + '/'
                                + participants[iSubj] + '/mri/aparc+aseg.mgz')

                # Target
                fname_T1w_nu1 = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_nu.mgz')
                fname_T1w_parc1 = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_aparc+aseg.mgz')

                os.system(f'cp {fname_T1w_nu0} {fname_T1w_nu1}')
                os.system(f'cp {fname_T1w_parc0} {fname_T1w_parc1}')

                # Converted
                fname_T1w_nu2 = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_nu.nii.gz')
                fname_T1w_parc2 = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_aparc+aseg.nii.gz')

                options = '-it mgz -ot nii'
                if os.path.exists(fname_T1w_nu2) == False:
                    os.system(f'FREESURFER mri_convert {options} {fname_T1w_nu1} {fname_T1w_nu2}')
                    os.remove(fname_T1w_nu1)
                if os.path.exists(fname_T1w_parc2) == False :
                    os.system(f'FREESURFER mri_convert {options} {fname_T1w_parc1} {fname_T1w_parc2}')
                    os.remove(fname_T1w_parc1)

                options = ''
                os.system(f'FSL fslreorient2std {options} {fname_T1w_nu2} {fname_T1w_nu2}')
                os.system(f'FSL fslreorient2std {options} {fname_T1w_parc2} {fname_T1w_parc2}')

                # Create masks
                fname_T1w_c1_m = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_nu_c1.nii.gz')
                fname_T1w_c2_m = (fdir_derivatives + participants[iSubj] + '_' + session + '_T1w_nu_c2.nii.gz')

                dat_T1_nu    = nib.load(fname_T1w_nu2)

                #   Grey matter   
                dat_T1_aparc = nib.load(fname_T1w_parc2)
                dat_T1_aparc = np.array(dat_T1_aparc.dataobj)
                dat_T1_c1_mask_rh = dat_T1_aparc
                dat_T1_c1_mask_rh[dat_T1_c1_mask_rh > 2036] = 0
                dat_T1_c1_mask_rh[dat_T1_c1_mask_rh < 2000] = 0

                dat_T1_aparc = nib.load(fname_T1w_parc2)
                dat_T1_aparc = np.array(dat_T1_aparc.dataobj)
                dat_T1_c1_mask_lh = dat_T1_aparc
                dat_T1_c1_mask_lh[dat_T1_c1_mask_lh > 1036] = 0
                dat_T1_c1_mask_lh[dat_T1_c1_mask_lh < 1000] = 0

                dat_T1_c1_mask = dat_T1_c1_mask_rh + dat_T1_c1_mask_lh
                dat_T1_c1_mask[dat_T1_c1_mask >= 1] = 1

                img = nib.Nifti1Image(dat_T1_c1_mask, dat_T1_nu.affine)
                img.get_data_dtype() == np.dtype(np.int16)
                nib.save(img, fname_T1w_c1_m)


                #   White matter   
                dat_T1_aparc = nib.load(fname_T1w_parc2)
                dat_T1_aparc = np.array(dat_T1_aparc.dataobj)
                dat_T1_c2_mask_rh  = dat_T1_aparc
                dat_T1_c2_mask_rh[dat_T1_c2_mask_rh > 4036] = 0
                dat_T1_c2_mask_rh[dat_T1_c2_mask_rh < 4000] = 0

                dat_T1_aparc = nib.load(fname_T1w_parc2)
                dat_T1_aparc = np.array(dat_T1_aparc.dataobj)
                dat_T1_c2_mask_lh = dat_T1_aparc
                dat_T1_c2_mask_lh[dat_T1_c2_mask_lh > 3036] = 0
                dat_T1_c2_mask_lh[dat_T1_c2_mask_lh < 3000] = 0

                dat_T1_c2_mask = dat_T1_c1_mask_rh + dat_T1_c1_mask_lh
                dat_T1_c2_mask[dat_T1_c2_mask >= 1] = 1

                img = nib.Nifti1Image(dat_T1_c2_mask, dat_T1_nu.affine)
                img.get_data_dtype() == np.dtype(np.int16)
                nib.save(img, fname_T1w_c2_m)


    return print('\n[ DONE ]\n\n')
