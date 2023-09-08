#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)

import numpy as np
import pandas as pd
import os
import csv
import glob
import zipfile
from dv_code.scripts.misc import read_analysis_params 

def extract_events_csv(fname_evt=None, subjID=None, ses=None, params=None):

    fdir_evt_src = os.path.dirname(fname_evt)
    fileexts = ['.zip'] # ['.tar.gz', '.bz2', '.zip']
    for ext in fileexts:
        if fname_evt.endswith(ext):
            print ('The file: ', fname_evt, ' has the extension: ', ext)
            
            print(f'# Unpacking: {fname_evt}')
            with zipfile.ZipFile(fname_evt,"r") as zip_ref:
                fname_inside_evt = fdir_evt_src + '/' + zip_ref.filelist[0].filename
                zip_ref.extractall(fdir_evt_src)
            print('done...\n')
                    
            tmp_csv   = []
            tmp_trash = []
            for jtem in zip_ref.filelist:
                tmp_jtem = fdir_evt_src + '/' + jtem.filename
                if '.csv' in tmp_jtem:
                    tmp_csv.append(fdir_evt_src + '/' + jtem.filename)
                else:
                    tmp_trash.append(fdir_evt_src + '/' + jtem.filename)
                    
            for jtem in tmp_trash:
                os.remove(jtem)
            
            fname_evt_ = tmp_csv[0]
            df = pd.read_csv(fname_evt_)
            
            # extract the relevant input
            df_relevant = df[['type','text_welcome.started', 'image_1.started', 'duration_null_event']]
            
    N_stims = df_relevant.shape[0]
    N_runs =  round( N_stims / 24 )
    tsv_header_name = ['onset', 'duration', 'trial_type']
        
    fdir_subjBIDS = params.get('fdir_data') + '/' +  subjID  + '/' + ses + '/func/'
    
    for i_run in np.arange(N_runs):
        
        dat_w_event_tsv = []
        fname_w_tsv = fdir_subjBIDS + subjID + '_' + ses + '_run-0' + (str(i_run+1)) + '_bold.tsv'
        for i_stim in range((i_run) * 24, (i_run + 1) * 24, 1):
            
            if i_run == 0:
                run_offset = df_relevant['text_welcome.started'][i_run * 24]
            else:
                run_offset = df_relevant['image_1.started'][i_run * 24] -10

                
            item_onset = (df_relevant['image_1.started'][i_stim] - run_offset)
            item_onset_ = str(round(item_onset, 4))
            item_trial_type = df_relevant['type'][i_stim]
            dat_w_event_tsv.append([item_onset_, '4', item_trial_type])
            
            item_onset_null = str(round(item_onset+4, 4))
            item_null_dur_ = str( df_relevant['duration_null_event'][i_stim] )
            dat_w_event_tsv.append([item_onset_null, item_null_dur_, 'null_event'])
        
        with open(fname_w_tsv, 'w') as tsvfile:    #csv writer to write in tsv file
            tsv_writer = csv.writer(tsvfile, delimiter='\t')    #write header in tsv file
            tsv_writer.writerow(tsv_header_name)    #write rows
            tsv_writer.writerows(dat_w_event_tsv)    #close csv file
            tsvfile.close()
    
    
    fname_w_tsv     = fdir_subjBIDS + subjID[1] + '_' + ses + '_run-99_bold.tsv'
    dat_w_event_tsv = []
    pass_run_time   = 0
    for i_run in np.arange(N_runs):
        
        for i_stim in range((i_run) * 24, (i_run + 1) * 24, 1):
            
            run_offset = df_relevant['text_welcome.started'][i_run * 24]

            item_onset = (df_relevant['image_1.started'][i_stim] - run_offset + pass_run_time)
            item_onset_ = str(round(item_onset, 4))
            item_trial_type = df_relevant['type'][i_stim]
            dat_w_event_tsv.append([item_onset_, '4', item_trial_type])
            
            item_onset_null = str(round(item_onset+4, 4))
            item_null_dur_ = str( df_relevant['duration_null_event'][i_stim] )
            dat_w_event_tsv.append([item_onset_null, item_null_dur_, 'null_event'])
        
        
        item_onset_welcome = str((float(item_onset_null) + float(item_null_dur_)))
        item_welcome_dur_ = str( 10 )
        dat_w_event_tsv.append([item_onset_welcome, item_welcome_dur_, 'welcome_'])
        pass_run_time = float(item_onset_welcome) + float(item_welcome_dur_)
    
    with open(fname_w_tsv, 'w') as tsvfile:    #csv writer to write in tsv file
        tsv_writer = csv.writer(tsvfile, delimiter='\t')    #write header in tsv file
        tsv_writer.writerow(tsv_header_name)    #write rows
        tsv_writer.writerows(dat_w_event_tsv)    #close csv file
        tsvfile.close()
        
    return