#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os
import csv
import numpy as np
import pandas as pd
from dv_code.scripts.misc.feedback_messages import error_message


def update_participants(params):
    
    #tsv_header_name = ['participant_id', 'Session', 'age', 'sex', 'StudyID']
    tsv_header_name = ['participant_id', 'session', 'studyid']
    study_id = params.get('study_id')

    sessions = os.listdir(params.get('fdir_sourcedata'))
    
    dat_w_participants_tsv = []
    for ises in sessions:
        
        if ises[:3] == 'ses':

            subjects_dir_entry = os.listdir(params.get('fdir_sourcedata') + '/' + ises)
            subjects_dir_entry = np.sort( subjects_dir_entry )

            fname_dcm_dat = [
                item
                for item in subjects_dir_entry
                if item[-4:] == '.zip' and item[-9:] != '_logs.zip'
            ]
            
            for i_fname in fname_dcm_dat:
                subj_info = i_fname.split('_')
                participant_id = f'sub-{str(subj_info[0])}'
                ses     = subj_info[1]
                studyID = subj_info[2]
                
                dat_w_participants_tsv.append([str(participant_id), str(ses), str(studyID)])

    fname_participants_tsv = params.get('fdir_data') + 'participants.tsv'
    if os.path.exists(fname_participants_tsv) == False:
        print((f'No participant file present.\n' +
                        f'Creating one at:\n{fname_participants_tsv}'))
        open(fname_participants_tsv, 'x')

        with open(fname_participants_tsv, 'w') as tsvfile:    #csv writer to write in tsv file
            tsv_writer = csv.writer(tsvfile, delimiter='\t')    #write header in tsv file
            tsv_writer.writerow(tsv_header_name)    #write rows
            tsv_writer.writerows(dat_w_participants_tsv)    #close csv file
            tsvfile.close()

    else:
        print((f'Participant file for study "{study_id}" present.\n' +
                f'At: {fname_participants_tsv}'))
        print('Updating it..')

        with open(fname_participants_tsv, 'w') as tsvfile:    #csv writer to write in tsv file
            tsv_writer = csv.writer(tsvfile, delimiter='\t')    #write header in tsv file
            tsv_writer.writerow(tsv_header_name)    #write rows
            tsv_writer.writerows(dat_w_participants_tsv)    #close csv file
            tsvfile.close()

        print("participants.tsv file updated...")


    return print('[ DONE ]')



def get_participants(params):
#
#
#
#
#
#
#
    tsv_header_name = ['participant_id', 'Session', 'StudyID']

    fname_participants_tsv = params.get('fdir_data') + 'participants.tsv'
    if (os.path.exists(fname_participants_tsv) == False):
        error_message('File does not exist')
    else:
        rd = pd.read_csv(fname_participants_tsv, sep='\t')

    return [rd.loc[item].participant_id for item in np.arange(rd.shape[0])]