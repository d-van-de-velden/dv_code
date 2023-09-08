#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#          Alexander Enge          (XX@YY.de)
#          Anne-Sophie Kieslinger  (XX@YY.de)
#
# License: BSD (3-clause)
# 
# 
#
# 
import os
import numpy as np
import json


def make_subjectlist(params):
    
    fdir_analysis = params.get('fdir_data')
    
    BIDS_subjects = os.listdir(params.get('fdir_data'))
    BIDS_subjects = np.sort( BIDS_subjects[1::] )
    # Browse through BIDS directory and gather information and paths from subjects 
    for i_subj in BIDS_subjects:
        tmp_dir_subj0      = fdir_analysis + '/' + i_subj
        tmp_session = os.listdir(tmp_dir_subj0)
        session = []
        for idx, i_item in enumerate(tmp_session):
            if 'ses-' in i_item:
                session.append(i_item)

        subject = dict()
        # Browse through the subject and create a subject info dictionary
        #subject = {"subject_ID": i_subj, "study_id" : "COCOA",
        #        "subject_Nses": len(session), "age": 33, "sex": 0}
        
        for i_sess in session:
            print('Subject :    ' + i_subj)
            print('Session :    ' + i_sess)
            tmp_dir_subj       = tmp_dir_subj0 + '/' + i_sess
            tmp_dir_subj_anat  = tmp_dir_subj + '/anat'
            tmp_dir_subj_func  = tmp_dir_subj + '/func'
            tmp_dir_subj_fmap  = tmp_dir_subj + '/fmap'

            tmp_ses_info = {"meas_date": "DUMMY", 
                            "fname_t1w": '[]', "fname_t2w": '[]',
                            "fname_func": '[]', "fname_fmap": '[]'}

            files_anat = os.listdir(tmp_dir_subj_anat)
            for i_file in files_anat:
                if "FLAIR" not in i_file:
                    if ".nii" in i_file:
                        tmp_ses_info["fname_t2w"] = tmp_dir_subj_anat + '/' + i_file
                if "T1w" not in i_file:
                    if ".nii" in i_file:
                        tmp_ses_info["fname_t1w"] = tmp_dir_subj_anat + '/' + i_file

            files_fmap = os.listdir(tmp_dir_subj_fmap)
            tmp_ = []
            for i_file in files_fmap:
                if ".nii" in i_file:
                    tmp_.append(tmp_dir_subj_fmap + '/' + i_file)
            tmp_ses_info["fname_fmap"] = tmp_

            files_func = os.listdir(tmp_dir_subj_func)
            tmp_ = []
            for i_file in files_func:
                if ".nii" in i_file:
                    tmp_.append(tmp_dir_subj_func + '/' + i_file)
            tmp_ses_info["fname_func"] = tmp_

            subject[i_sess] = tmp_ses_info
        
        # create json object from dictionary
        tmp_json = json.dumps(subject)
        # open file for writing, "w" 
        f = open((fdir_analysis + '/' + i_subj + "/subject_info.json"),"w")
        # write json object to file
        f.write(tmp_json)
        # close file
        f.close()
        
    return print(' [ DONE ]')