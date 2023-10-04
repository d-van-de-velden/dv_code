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
import pathlib
import shutil
import zipfile
from dv_code.scripts.preprocessing.func_extract_evts import extract_events_csv
import numpy as np

from dv_code.scripts.work_participants import update_participants
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.misc.use_utils import unpack_compressed
from dv_code.scripts.work_participants import get_participants
from dv_code.scripts.misc.analysis_feedback import loadingBar

def convert_dcm2BIDS(params):
    
    
    fname_config =  params.get('fdir_analysis') + '/conf/configs/dcm2bids_config.json'

    # ----------   Check if nifits are present for all subjects -----
    print('Identifying sessions from "sourcedata"....')
    sessions = os.listdir(params.get('fdir_sourcedata'))

    for session in sessions:

        # Loop over files for the current session
        if session[:4] == 'ses-':
            tmp_dir = params.get('fdir_sourcedata') + '/' + session
            subjects_dir_entry = os.listdir(tmp_dir)
            subjects_dir_entry = np.sort( subjects_dir_entry )

            fname_dcm_dat = []
            for item in subjects_dir_entry:
                if (item[-4:] == '.zip' and item[-9:] != '_logs.zip') and 'Behavioral' not in item:
                    print('   Subject : ' + item[:4] + '  |  Session : ' + session)

                    fname_dcm_dat.append(item)
                    
                    source_file = tmp_dir + '/' + item
                    print(f'# Unpacking: {source_file}')
                    with zipfile.ZipFile(source_file,"r") as zip_ref:
                        fname_inside_dcm = tmp_dir + '/' + zip_ref.filelist[0].filename
                        zip_ref.extractall(tmp_dir)
                    print('done...\n')

                    tmp_dir_subj_dicom = fname_inside_dcm + 'DICOM'
                    tmp_dir_subj_nifti = fname_inside_dcm + '/niftidir'
                    if os.path.exists(tmp_dir_subj_nifti) == False:
                        check_make_dir(tmp_dir_subj_nifti)

                    cmd_convert2BIDS = ('dcm2bids' 
                                        + ' -d ' + tmp_dir_subj_dicom
                                        + ' -p ' + item[:4] 
                                        + ' -s ' + session 
                                        + ' -c ' + fname_config
                                        + ' -o ' + params.get('fdir_data')
                                        + ' --clobber'
                                        + ' --force_dcm2bids'
                                        )



                    fdir_BIDS_subj = params.get('fdir_data') + 'sub-' + item[:4] + '/' + session + '/'
                    fdir_BIDS_dir_entry = os.listdir(fdir_BIDS_subj)
                    do_convert = False
                    if os.path.exists(fdir_BIDS_subj) == False:
                        do_convert = True
                        check_make_dir(fdir_BIDS_subj)
                    elif len(fdir_BIDS_subj) == 0:
                        do_convert = True
                    elif 'func' not in fdir_BIDS_dir_entry or 'fmap' not in fdir_BIDS_dir_entry:
                        do_convert = True


                    if do_convert:
                        print('Now convert Nifties to BIDS conform directory structure...')
                        d = os.system(cmd_convert2BIDS)

                    print('\nCleaning directory from unzipped files...')
                    shutil.rmtree(fname_inside_dcm)
                    print('done..\n')

    update_participants(params)

    return print('\n [ DONE ]\n\n')


def convert_evt2BIDS(params):
    
    # ----------   Check if data is present for all subjects -----
    print('Identifying sessions from "sourcedata"....')
    sessions = os.listdir(params.get('fdir_sourcedata'))

    for session in sessions:

        # Loop over files for the current session
        if session[:4] == 'ses-':
            tmp_dir = params.get('fdir_sourcedata') + '/' + session
            subjects_dir_entry = os.listdir(tmp_dir)
            subjects_dir_entry = np.sort( subjects_dir_entry )

            fname_evt_dat = []
            for item in subjects_dir_entry:
                if item[-9:] == '_logs.zip' and 'Behavioral' not in item:
                    
                    indices = item.split('_')
                    subjID    = 'sub-' + indices[0]
                    ses     = indices[1].replace('-0','-')
                    study   = indices[2] 
                    print(f'   Subject : {subjID}  |  Session : {ses} |  Study : {study} ')

                    tmp_fdir_subj = params.get('fdir_data') + subjID
                    tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'
                    check_make_dir(tmp_fdir_func)
                    fname_evt_dat.append(item)
                    
                    source_file = tmp_dir + '/' + item
                    print(f'# Unpacking: {source_file}')
                    with zipfile.ZipFile(source_file,"r") as zip_ref:
                        zip_ref.extractall(tmp_dir)
                    print('done...\n')
                    
                    tmp_csv   = []
                    tmp_trash = []
                    for jtem in zip_ref.filelist:
                        
                        tmp_jtem = tmp_dir + '/' + jtem.filename
                        if '.csv' in tmp_jtem:
                            tmp_csv.append(tmp_dir + '/' + jtem.filename)
                        else:
                            tmp_trash.append(tmp_dir + '/' + jtem.filename)
                    
                    for jtem in tmp_trash:
                        os.remove(jtem)
                    
                    fname_evt_zip = source_file
                    extract_events_csv(fname_evt_zip, subjID, ses, params)
                    
                    for jtem in tmp_csv:
                        os.remove(jtem)
                    
                    
    
    
    return


def convert_beh2BIDS(params):
    
    fdir_data      = params.get('fdir_data')
    
    # ----------   Check if nifits are present for all subjects -----
    print('Identifying sessions from "sourcedata"....')
    sessions = os.listdir(params.get('fdir_sourcedata'))

    for session in sessions:

        # Loop over files for the current session
        if session[:4] == 'ses-':
            tmp_dir = params.get('fdir_sourcedata') + '/' + session
            subjects_dir_entry = os.listdir(tmp_dir)
            subjects_dir_entry = np.sort( subjects_dir_entry )

            for item in subjects_dir_entry:
                if (item[-4:] == '.zip' or item[-7:] == '.tar.gz') and 'Behavioral' in item:
                    print(item)
                    
                    tmp_fname_idents = os.path.splitext(os.path.basename(item))[0]
                    tmp_fname_idents = os.path.splitext(os.path.basename(tmp_fname_idents))[0]

                    idents  = tmp_fname_idents.split('_')
                    subjID  = idents[0]
                    tmp_str_BIDS_session = idents[1].replace("-0", "-")
                    studyID = idents[2]
                    studyID = studyID + '-Behavioral'
                    
                    fdir_BIDS_behave = f'{fdir_data}sub-{subjID}/{tmp_str_BIDS_session}/beh/'
                    shutil.rmtree(fdir_BIDS_behave, ignore_errors=True) 
                    check_make_dir(fdir_BIDS_behave)
                    
                    file_input = tmp_dir + '/' + item
                    unpack_compressed(finput=file_input, fdir_out=fdir_BIDS_behave)
                    
                    for jtem in os.listdir(fdir_BIDS_behave):
                        tmp_fname_idents0 = os.path.splitext(os.path.basename(jtem))[0]
                        tmp_fname_idents  = os.path.splitext(os.path.basename(tmp_fname_idents0))[0]
                        fend              = pathlib.Path(jtem).suffix

                        idents  = tmp_fname_idents.split('_')
                                                
                        if idents[0] == 'AE':
                            print(f'WRONG LOG FILE FOUND !!!! \n# Removing it: \n     File : {jtem}')
                            wrong_file = os.path.join(fdir_BIDS_behave, jtem)
                            os.remove(wrong_file)
                            
                        if idents[0] == 'BE':
                            subjID  = idents[1]
                            tmp_str_BIDS_session = idents[2].replace("-0", "-")
                            studyID = idents[3]
                            studyID = studyID + '-Behavioral'
                            
                        else:
                            subjID  = idents[0]
                            tmp_str_BIDS_session = idents[1].replace("-0", "-")
                            studyID = idents[2]
                            studyID = studyID + '-Behavioral'
                            
                        fname_rename = f'sub-{subjID}_{session}_{studyID}{fend}'
                        old_filename = os.path.join(fdir_BIDS_behave, jtem)
                        new_filename = os.path.join(fdir_BIDS_behave, fname_rename)
                        print('# Renameing:')
                        print(f'    From : {old_filename}')
                        print(f'    To   : {new_filename}')

                        os.rename(old_filename, new_filename)
                        


    update_participants(params)

    participants = get_participants(params)
    check_make_dir(params.get('fdir_proc_stat'))

    fdir_results  = params.get('fdir_proc_stat') + '/group/'
    check_make_dir(fdir_results)      

    list_indicator_session = []
    for iSubj in range(len(participants)):
        tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

        sessions = os.listdir(tmp_fdir_subj)

        for session in sessions:
            list_indicator_session.append(session)
            tmp_fdir_beh_src  = f'{tmp_fdir_subj}/{session}/beh/'
            tmp_fdir_beh_trgt = f'{fdir_results}/{session}/beh/'
            check_make_dir(tmp_fdir_beh_trgt)  

            if os.path.exists(tmp_fdir_beh_src) == True:
                fnames = os.listdir(tmp_fdir_beh_src)

                if fnames:
                    for fname in fnames:
                        # Only use .csv files
                        if fname[-4:] == '.csv':
                            loadingBar(iSubj, len(participants), task_part=participants[iSubj])
                            tmp_fname = fname
                            tmp_fname_beh = os.path.splitext(os.path.basename(tmp_fname))[0]
                            fname_beh = os.path.splitext(os.path.basename(tmp_fname_beh))[0]

                            idents  = fname_beh.split('_')
                            subjID  = idents[0]
                            session = idents[1]

                            # Copy the csv file
                            fname_csv_src  = ( tmp_fdir_beh_src  + tmp_fname)
                            fname_csv_trgt = ( tmp_fdir_beh_trgt + tmp_fname)
                            
                            shutil.copyfile(fname_csv_src, fname_csv_trgt)
                            

    return print('\n [ DONE ]\n\n')