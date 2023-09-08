#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os
import re
import stat
import json

from PIL import Image

from dv_code.scripts.misc.select_filetype import get_only_nifti, get_only_json
from dv_code.scripts.misc.feedback_messages import error_message
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.processing.func_MRI_quality_metric import calc_tSNR, motion_evaluation, get_tSNR
from dv_code.scripts.viz.combine_plots import viz_concat_v


def do_func_preproc(fdir_analysis="", participants=None, 
                    params=None, forceRun=True, use_HPC=False):
    if participants is None:
        participants = []
    check_make_dir(params.get('fdir_proc_pre'))

    if not fdir_analysis:
        input_str = "Variable 1 (fdir_analysis) is not provided. ABORT"
        error_message(input_str)

    env_path = os.getenv('PATH')
    host = '.'
    if 'cbs.mpg' in env_path:
        host = 'mpg'
        
    print(
        f'Preprocessing of functional MR data is initiated for: {len(participants)} subjects'
    )

    for iSubj in range(len(participants)):
        print(f'Identifying sessions from "data/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

        sessions = os.listdir(tmp_fdir_subj)
        for session in sessions:
            print(session)
            tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'
            tmp_fdir_fmap = f'{tmp_fdir_subj}/{session}/fmap/'

            tmp_fs = os.listdir(tmp_fdir_fmap)
            tmp_fs = get_only_json(tmp_fs)
            with open(f'{tmp_fdir_fmap}/{tmp_fs[0]}') as f:
                fmap_json = json.load(f)
            TotalReadOutTime = fmap_json.get('TotalReadoutTime')

            runs = os.listdir(tmp_fdir_func)
            runs = get_only_nifti(runs, 1)
            print(runs)
            runs.sort()
            
            if len(runs) > 0:
                fdir_derivatives_func = (params.get('fdir_proc_pre')
                                    + '/' + participants[iSubj] 
                                    + '/' + session + '/func'
                                    )
                check_make_dir(fdir_derivatives_func)
                fdir_derivatives_anat = (params.get('fdir_proc_pre')
                                    + '/' + participants[iSubj] 
                                    + '/' + session + '/anat'
                                    )
                tmp_fname_T1w = f'{fdir_derivatives_anat}/{participants[iSubj]}_{session}_T1w_nu.nii.gz'

            if os.path.exists(tmp_fname_T1w):
            
                for run in runs:
                    tmp_fname_funcMR = f'{tmp_fdir_func}/{run}'

                    tmp_fname_final_funcMR = (fdir_derivatives_func
                                                + '/' + str(run)[:-7]
                                                + '_st_mcf_topUP.nii.gz')
                    if (os.path.exists(tmp_fname_final_funcMR) == False) or (forceRun == True):
                        print('Preprocessed file for functional data does not exist...')

                        if os.path.exists(tmp_fname_final_funcMR) == True and forceRun == True:
                            print('''!!! Preprocessed file for functional data DID exist.
                                    \nFORCING rerun ....''')

                        if os.path.exists(tmp_fname_funcMR):
                            print(tmp_fname_funcMR)

                            fdir_bash_script = (params.get('fdir_bash') 
                                                + '/' + participants[iSubj]
                                                + '/' + session
                                                + '/func_preprocessing/')
                            fdir_lib_bash = params.get('fdir_lib_bash')
                            check_make_dir(fdir_bash_script)
                            fname_bash_script = (fdir_bash_script
                                                + '/run_funcPreproc.sh' 
                                                )
                            fname_bash_script_ = (fdir_lib_bash
                                                + '/runner_func_preproc.sh'
                                                )
                            check_make_dir(fdir_bash_script)

                            with open(fname_bash_script, 'w') as rsh:
                                    rsh.write('''
                                    #! /bin/sh
                                    echo "I start AFNI + FSL now...."
                                    FREESURFER AFNI FSL bash ''' + fname_bash_script_
                                    + ''' -s ''' + participants[iSubj]
                                    + ''' -e ''' + session
                                    + ''' -f ''' + tmp_fname_funcMR
                                    + ''' -u ''' + tmp_fdir_func
                                    + ''' -l ''' + tmp_fdir_fmap
                                    + ''' -t ''' + str(TotalReadOutTime)
                                    + ''' -a ''' + tmp_fname_T1w
                                    + ''' -o ''' + fdir_derivatives_func 
                                    + ''' \n exit \n quit \n q 
                                    ''')

                            os.chmod(fname_bash_script, stat.S_IRWXU)
                            os.chmod(fname_bash_script_, stat.S_IRWXU)
                            if use_HPC == 1:
                                print('Use HPC cluster for computation...')
                            else:
                                print('Use local machine for computation...')
                                fname_block_file = f'{fdir_bash_script}/blocked.txt'

                                if os.path.exists(fname_block_file):
                                    print('Job already running....')
                                else:
                                    with open(fname_block_file, 'w') as rsh:
                                        rsh.write("BLOCKED JOB\nfname= " + fname_bash_script)

                                    # To debug

                                    if host == 'mpg':
                                        os.system(f'xterm -geometry 60x80 -e bash {fname_bash_script}')
                                    else:
                                        os.system(f'xterm -geometry 60x80 -e bash {fname_bash_script_}')
                                    os.remove(fname_block_file)

                                    prefix_del1 = os.path.basename(f'{tmp_fname_final_funcMR[:-7]}_r_')
                                    prefix_del2 = os.path.basename(f'{tmp_fname_final_funcMR[:-3]}_r')
                                    for f in os.listdir(fdir_derivatives_func):
                                        if re.match(prefix_del1, f):
                                            os.remove(f'{fdir_derivatives_func}/{f}')
                                        if re.match(prefix_del2, f):
                                            os.remove(f'{fdir_derivatives_func}/{f}')


                    else:
                        print('Preprocessed file for functional data does exist...')


                    # Perform image quality metric analysis
                    fname_par = (fdir_derivatives_func
                                + '/' + str(run)[:-7]
                                + '_st_mcf.par')

                    # Perform image quality metric analysis
                    fname_tsnr0 = calc_tSNR(tmp_fname_funcMR, params)
                    fname_tsnr1 = calc_tSNR(tmp_fname_final_funcMR, params)

                    fname_plot0 = motion_evaluation(fname_par, fd_threshold=0.3)
                    fname_plot1, tSNR_median1, ROI_label = get_tSNR(fname_tsnr0, params, doPlot=True)
                    fname_plot2, tSNR_median2, ROI_label = get_tSNR(fname_tsnr1, params, doPlot=True)

                    img0 = Image.open(fname_plot0)
                    img1 = Image.open(fname_plot1)
                    img2 = Image.open(fname_plot2)

                    fdir_ALL_qc = f'{os.path.dirname(fname_par)}/QC/'
                    fname_ALL_qc = f'QC_ALL_{os.path.basename(fname_par)[:-4]}.png'
                    img_a = viz_concat_v([img0, img1, img2])
                    img_a.save(( fdir_ALL_qc + fname_ALL_qc))


    return print('[ DONE ]')