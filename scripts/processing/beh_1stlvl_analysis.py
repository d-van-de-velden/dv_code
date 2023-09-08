#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
#
# import copy
# import math
# import os
# from pathlib import Path
# import shutil
# import stat
# from matplotlib import pyplot as plt
# import pandas as pd
# from dv_code.scripts.misc import check_make_dir
# from dv_code.scripts.misc.analysis_feedback import loadingBar
# from dv_code.scripts.work_participants import get_participants
# import numpy as np
# from scipy import stats  # Visualization
# import seaborn as sns
# from matplotlib import pyplot as plt
# from scipy.stats import norm
# import math


# def do_hddModelling_R(participants=None, params=None):
        
#     if participants is None:
#         participants = get_participants(params)
#     check_make_dir(params.get('fdir_proc_stat'))

#     fdir_lib_bash = params.get('fdir_lib_bash')
#     fdir_results  = params.get('fdir_proc_stat') + '/group/'
#     check_make_dir(fdir_results)      
        
#     print(
#         f'Performing hierachical Drift Diffusion Modelling on behavioral data with R for: {len(participants)} subjects\n'
#     )

#     print('  # Preparing all the behavioral data..\n')

#     list_indicator_session = []
#     for iSubj in range(len(participants)):
#         tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

#         sessions = os.listdir(tmp_fdir_subj)

#         for session in sessions:
#             list_indicator_session.append(session)
#             tmp_fdir_beh_src  = f'{tmp_fdir_subj}/{session}/beh/'
#             tmp_fdir_beh_trgt = f'{fdir_results}/{session}/beh/'
#             check_make_dir(tmp_fdir_beh_trgt)  

#             if os.path.exists(tmp_fdir_beh_src) == True:
#                 fnames = os.listdir(tmp_fdir_beh_src)

#                 if fnames:
#                     for fname in fnames:
#                         # Only use .csv files
#                         if fname[-4:] == '.csv':
#                             loadingBar(iSubj, len(participants), task_part=participants[iSubj])
#                             tmp_fname = fname
#                             tmp_fname_beh = os.path.splitext(os.path.basename(tmp_fname))[0]
#                             fname_beh = os.path.splitext(os.path.basename(tmp_fname_beh))[0]
                            
#                             idents  = fname_beh.split('_')
#                             subjID  = idents[0]
#                             session = idents[1]
#                             studID  = idents[2]
                            
#                             # Copy the csv file
#                             fname_csv_src  = ( tmp_fdir_beh_src  + tmp_fname)
#                             fname_csv_trgt = ( tmp_fdir_beh_trgt + tmp_fname)
                            
#                             shutil.copyfile(fname_csv_src, fname_csv_trgt)

#                             #print(f'  # Copy \n   {fname_csv_src} to \n   {fname_csv_trgt}')
                            

#     run_sessions = list( np.unique( np.array(list_indicator_session)) )
#     for i, session in enumerate(run_sessions):
#         loadingBar(i, len(run_sessions), task_part=session)
#         fdir_bash_script = (params.get('fdir_bash')
#                             + '/' + session
#                             + '/hBayesDDModelling/'
#                             )
#         fdir_beh_analysis = f'{fdir_results}/{session}/beh/'
#         fname_R_script = f'{fdir_lib_bash}/runner_hDDM.r'

#         check_make_dir(fdir_beh_analysis)
#         fname_R_script_ = f'{fdir_beh_analysis}/run_start_R_hBayesDDM.sh'
#         fname_ORIG_R_script = f'{fdir_lib_bash}/runner_hDDM.R'
#         fname_R_script = f'{fdir_beh_analysis}/run_R_hBayesDDM.sh'
#         check_make_dir(fdir_beh_analysis)

#         with open(fname_R_script, 'w') as rsh:
#                             rsh.write('''
#                             #! /bin/sh
#                             echo "I start R now...."
#                             R+ bash ''' + fname_R_script_
#                             )
#         with open(fname_R_script_, 'w') as rsh:
#                             rsh.write('#! /bin/sh\n' +
#                                     'echo I start R script now....\n' + 
#                                     f'cd {fdir_beh_analysis}\n' + 
#                                     f'Rscript  {fname_ORIG_R_script}')
        
#         os.chmod(fname_R_script, stat.S_IRWXU)
#         os.chmod(fname_R_script_, stat.S_IRWXU)
#         use_HPC = 0
#         if use_HPC == 1:
#             print('Use HPC cluster for computation...')
#         else:
#             print('Use local machine for computation...')
#             fname_block_file = f'{fdir_bash_script}/blocked.txt'

#             if os.path.exists(fname_block_file) == True:
#                 print('Job already running....')
#             else:
#                 with open(fname_block_file, 'w') as rsh:
#                     rsh.write("BLOCKED JOB\nfname= " + fname_R_script)

#                     if host == 'mpg':
#                         os.system(f"xterm -e bash {fname_R_script}")
#                     else:
#                         os.system(f"xterm -e bash {fname_R_script}")
#                         os.remove(fname_block_file)
                    
        
        

#     return

# def do_hddModelling(participants=None, params=None, niter=4000,
#                     do_QC_plot=True, forceRedo=True):
        
#     if participants is None:
#         participants = get_participants(params)
#     check_make_dir(params.get('fdir_proc_stat'))

#     fdir_results = params.get('fdir_proc_stat') + '/group/'
#     check_make_dir(fdir_results)      
        
#     print(
#         f'Performing hierachical Drift Diffusion Modelling on behavioral data for: {len(participants)} subjects\n'
#     )

#     print('  # Collecting all the behavioral data..\n')

#     fname_hDDM_results = f'{fdir_results}input_dataset_hDDModelling_N{len(participants)}_niter{niter}.csv'
#     if os.path.exists(fname_hDDM_results) == False or forceRedo == True:
    
#         dataset_csv = []
#         subj_counter = 0
#         subj_data_counter = []
#         NO_beh_data_available = []
#         for iSubj in range(len(participants)):
#             tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

#             sessions = os.listdir(tmp_fdir_subj)

#             for session in sessions:
                
#                 tmp_fdir_beh = f'{tmp_fdir_subj}/{session}/beh/'

#                 if os.path.exists(tmp_fdir_beh) == True:
#                     fnames = os.listdir(tmp_fdir_beh)

#                     if fnames:
#                         for fname in fnames:
#                             # Only use .csv files
#                             if fname[-4:] == '.csv':
#                                 loadingBar(iSubj, len(participants), task_part=participants[iSubj])
#                                 tmp_fname = fname
#                                 tmp_fname_beh = os.path.splitext(os.path.basename(tmp_fname))[0]
#                                 fname_beh = os.path.splitext(os.path.basename(tmp_fname_beh))[0]

#                                 idents  = fname_beh.split('_')
#                                 subjID  = idents[0]
#                                 session = idents[1]
#                                 studID  = idents[2]
                                
#                                 # Load the csv file
#                                 tmp_data= pd.read_csv(( tmp_fdir_beh + tmp_fname))
                                
#                                 # Get data
#                                 tmp_subj_rt       = list( tmp_data['key_resp.rt'] )
#                                 if np.sum( ( np.array(tmp_subj_rt) <= 0.3 ) ) > 0:
#                                     for i, item in enumerate(tmp_subj_rt):
#                                         if item <= 0.3:
#                                             tmp_subj_rt[i] = np.nan
                                
#                                 tmp_subj_response = list( tmp_data['key_resp.corr'] ) # 1= correct ; 0= incorrect
#                                 tmp_subj_stim     = list( tmp_data['type'] ) # 1= correct ; 0= incorrect
                                
#                                 # Remove responses, for which no reaction time is given
#                                 tmp_trash, tmp_subj_stim = filter_corr_nans(tmp_subj_rt, tmp_subj_stim)
#                                 tmp_subj_rt, tmp_subj_response = filter_corr_nans(tmp_subj_rt, tmp_subj_response)
                                
#                                 for i, item in enumerate(tmp_subj_response):
#                                         if item == 0:
#                                             tmp_subj_response[i] = -1
                                
#                                 if len(tmp_subj_response) != len(tmp_subj_rt):
#                                     print('ERROR')
#                                 else:
#                                     subj_data_counter.append( len(tmp_subj_rt) )
                                
#                                 data_single_subject = pd.DataFrame(
#                                     {
#                                         "rt": tmp_subj_rt,
#                                         "response": tmp_subj_response,
#                                         "session": session,
#                                         "type": tmp_subj_stim,
#                                     }
#                                     )
                                
#                                 # Append single subject data to cohort dataset 
#                                 dataset_csv.append(data_single_subject)
                                
#                                 subj_counter = subj_counter + 1
                
#                     else:
#                         NO_beh_data_available.append( participants[iSubj] )
#                 else:
#                         NO_beh_data_available.append( participants[iSubj] )
                                
        
#         # Make single dataframe out of subject-wise datasets
#         dataset = pd.concat(dataset_csv)
        
#         print( (f'\n# Got dataset' +
#                 f'\n      -> Columns      : {list(dataset.head(0))}' +
#                 f'\n      -> Subjects     : {subj_counter}' +
#                 f'\n      -> ' + r'$\SS-Input_max' + f' : {np.max(subj_data_counter)}' +
#                 f'\n      -> ' + r'$\SS-Input_min' + f' : {np.min(subj_data_counter)}' +
#                 f'\n      -> ' + r'$\SS-Input_mean' + f' : {np.mean(subj_data_counter)}'
#                 )
#             )
#         print( (f'\n# Data missing from' +
#                 f'\n      -> Subjects     : {np.unique( NO_beh_data_available )}')
#             )


#         # Save the dataset information before running in model
#         filepath = Path(fname_hDDM_results)  
#         filepath.parent.mkdir(parents=True, exist_ok=True)  
#         dataset.to_csv(filepath)  
            
#     else:
#         dataset = pd.read_csv(fname_hDDM_results, encoding='utf-8')
    
    
    
#     sessions = list( np.unique(dataset.session) )
    
#     HDDM_across_session = []
#     for session in sessions:
        
#         types = list( np.unique(dataset.type) )
#         types.append('total') 
#         for itype in types:
                
#             tmp_ses_dataset = dataset[dataset.session == session]
#             tmp_ses_dataset = tmp_ses_dataset.drop(columns=["session"])
            
#             if itype != 'total':
#                 tmp_ses_dataset = dataset[dataset.type == itype]
#             tmp_ses_dataset = tmp_ses_dataset.drop(columns=["type"])
        
#             ddm_model =  hssm.HSSM(data=tmp_ses_dataset)
            
#             warmup_use = niter * 0.2
#             niter_use  = niter * 0.8
            
#             print('  # Modelling .......')
#             data_ddm_model = ddm_model.sample(
#                 cores   = int( 4 ),             # how many cores to use
#                 chains  = int( 4 ),             # how many chains to run
#                 draws   = int( niter_use ),     # number of draws from the markov chain
#                 tune    = int( warmup_use ),    # number of burn-in samples
#                 )
#             print('  ### Done')
#             az.summary(data_ddm_model)
            
#             # Get the data from the model an PLOT
#             # Response caution (Alpha)
#             fig = plt.figure(dpi=300, facecolor='lightgrey', frameon=True, figsize=[9,8])
#             fig.suptitle(f"Hierarchical Drift Diffusion Modelling\n(N: {len(participants)}, Session: {session}), Stim: {itype}",
#                         fontweight='bold', fontsize=16)
#             Alpha = np.array( data_ddm_model.posterior.data_vars['a'] )
#             plt.subplot(221)
#             plt.title('Response caution (Alpha)', weight='bold')
#             tmp_data = Alpha
#             y, binEdges = np.histogram(tmp_data, bins=100)
#             bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
#             plt.plot(bincenters, y, '-', c='blue', linewidth=3)
#             plt.axvline(tmp_data.mean(), color='k', linestyle='dashed', linewidth=1)
#             plt.xlabel('Alpha [a.U.]')
#             plt.ylabel('Density')
            
#             # Response bias (Beta)
#             Beta  = np.array( data_ddm_model.posterior.data_vars['z'] )
#             plt.subplot(222)
#             plt.title('Response bias (Beta)', weight='bold')
#             tmp_data = Beta
#             y, binEdges = np.histogram(tmp_data, bins=100)
#             bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
#             plt.plot(bincenters, y, '-', c='green', linewidth=3)
#             plt.axvline(tmp_data.mean(), color='k', linestyle='dashed', linewidth=1)
#             plt.xlabel('Beta [a.U.]')
#             plt.ylabel('Density')
#             plt.xlim([0, 1])

#             # DriftRate (Delta)
#             Delta = np.array( data_ddm_model.posterior.data_vars['v'] )
#             plt.subplot(223)
#             plt.title('DriftRate (Delta)', weight='bold')
#             tmp_data = Delta
#             y, binEdges = np.histogram(tmp_data, bins=100)
#             bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
#             plt.plot(bincenters, y, '-', c='red', linewidth=3)
#             plt.axvline(tmp_data.mean(), color='k', linestyle='dashed', linewidth=1)
#             plt.xlabel('Delta [a.U.]')
#             plt.ylabel('Density')
            
#             # Sensory encoding & motor enxecution (Tau)
#             Tau   = np.array( data_ddm_model.posterior.data_vars['t'] )
#             plt.subplot(224)
#             plt.title('Sensory encoding & motor execution (Tau)', weight='bold')
#             tmp_data = Tau
#             y, binEdges = np.histogram(tmp_data, bins=100)
#             bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
#             plt.plot(bincenters, y, '-', c='violet', linewidth=3)
#             plt.axvline(tmp_data.mean(), color='k', linestyle='dashed', linewidth=1)
#             plt.xlabel('time [s]')
#             plt.ylabel('Density')
#             plt.xlim([0.2, 0.7])
            
#             plt.tight_layout()
#             fname_ses_overview_hddm = f'{fdir_results}hDDModelling_results_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#             plt.savefig(fname_ses_overview_hddm)

#             if do_QC_plot:
#                 fig = plt.figure
#                 az.plot_trace(data_ddm_model)
#                 plt.tight_layout()
#                 fname_fig_trace = f'{fdir_results}QC2_trace_hDDModelling_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#                 plt.savefig(fname_fig_trace)
            
#             # Append data from each session
#             data_single_session = pd.DataFrame(
#                                 {
#                                 "v": np.ndarray.flatten(Delta),
#                                 "t": np.ndarray.flatten(Tau),
#                                 "z": np.ndarray.flatten(Beta),
#                                 "a": np.ndarray.flatten(Alpha),
#                                 "session": session,
#                                 "type": itype,
#                                 }
#                                 )
                                    
#             # Append single subject data to cohort dataset 
#             HDDM_across_session.append(data_single_session)
                                
                                
#     dataframe_HDDM_ses = pd.concat(HDDM_across_session)
    
    
#     fig = plt.figure(figsize=[12,12], dpi=300)
#     # Initialize the FacetGrid object
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataframe_HDDM_ses, row="session", hue='type',
#                             aspect=9, height=1.5, legend_out=True)
#     fig_obj.map_dataframe(sns.kdeplot, x="a", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Response caution (Alpha)', y=0.98)
#     legend = plt.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
#     fname_plot1 = fdir_results + f'/hDDModelling_results-Alpha_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#     fig_obj.savefig(fname_plot1)
    
#     fig = plt.figure(figsize=[12,12], dpi=300, facecolor='lightgrey')
#     # Initialize the FacetGrid object
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataframe_HDDM_ses, row="session", hue='type',
#                             aspect=9, height=1.5)
#     fig_obj.map_dataframe(sns.kdeplot, x="z", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Response bias (Beta)', y=0.98)
#     legend = plt.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
#     fname_plot2 = fdir_results + f'/hDDModelling_results-Beta_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#     fig_obj.savefig(fname_plot2)
    
#     fig = plt.figure(figsize=[12,12], dpi=300, facecolor='lightgrey')
#     # Initialize the FacetGrid object
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataframe_HDDM_ses, row="session", hue='type',
#                             aspect=9, height=1.5)
#     fig_obj.map_dataframe(sns.kdeplot, x="v", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Drift Rate (Delta)', y=0.98)
#     legend = plt.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
#     fname_plot3 = fdir_results + f'/hDDModelling_results-Delta_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#     fig_obj.savefig(fname_plot3)
    
#     fig = plt.figure(figsize=[12,12], dpi=300, facecolor='lightgrey')
#     # Initialize the FacetGrid object
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataframe_HDDM_ses, row="session", hue='type',
#                             aspect=9, height=1.5)
#     fig_obj.map_dataframe(sns.kdeplot, x="t", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Response caution (Tau)', y=0.98)
#     legend = plt.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
#     fname_plot4 = fdir_results + f'/hDDModelling_results-Tau_N{len(participants)}_{session}_stim-{itype}_niter{niter}.png'
#     fig_obj.savefig(fname_plot4)
    
#     return


# def do_behave_feature_analysis(participants=None, params=None,
#                     features=None, do_QC_plot=True, forceRedo=True):


#     if participants is None:
#         participants = get_participants(params)
#     check_make_dir(params.get('fdir_proc_stat'))

#     fdir_results = params.get('fdir_proc_stat') + '/group/'
#     check_make_dir(fdir_results)      
        
#     print(
#         f'Performing Feature analysis on behavioral data for: {len(participants)} subjects\n'
#     )

#     print('  # Collecting all the behavioral data..\n')

#     fname_results = f'{fdir_results}input_dataset_Beh-Features_N{len(participants)}.csv'
#     if os.path.exists(fname_results) == False or forceRedo == True:
    
#         dataset_csv = []
#         subj_counter = 0
#         subj_data_counter = []
#         NO_beh_data_available = []
#         for iSubj in range(len(participants)):
#             tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

#             sessions = os.listdir(tmp_fdir_subj)

#             for session in sessions:
                
#                 tmp_fdir_beh = f'{tmp_fdir_subj}/{session}/beh/'

#                 if os.path.exists(tmp_fdir_beh) == True:
#                     fnames = os.listdir(tmp_fdir_beh)

#                     if fnames:
#                         for fname in fnames:
#                             # Only use .csv files
#                             if fname[-4:] == '.csv':
#                                 loadingBar(iSubj+1, len(participants), task_part=participants[iSubj])
#                                 tmp_fname = fname
#                                 tmp_fname_beh = os.path.splitext(os.path.basename(tmp_fname))[0]
#                                 fname_beh = os.path.splitext(os.path.basename(tmp_fname_beh))[0]

#                                 idents  = fname_beh.split('_')
#                                 subjID  = idents[0]
#                                 session = idents[1]
#                                 studID  = idents[2]
                                
#                                 # Load the csv file
#                                 tmp_data= pd.read_csv(( tmp_fdir_beh + tmp_fname))
                                
#                                 # Get data
#                                 tmp_subj_rt       = list( tmp_data['key_resp.rt'] )
#                                 if np.sum( ( np.array(tmp_subj_rt) <= 0.3 ) ) > 0:
#                                     for i, item in enumerate(tmp_subj_rt):
#                                         if item <= 0.3:
#                                             tmp_subj_rt[i] = np.nan
                                
#                                 tmp_subj_response = list( tmp_data['key_resp.corr'] ) # 1= correct ; 0= incorrect
#                                 tmp_subj_stim     = list( tmp_data['type'] ) # 1= correct ; 0= incorrect
                                
#                                 # Remove responses, for which no reaction time is given
#                                 tmp_trash, tmp_subj_stim = filter_corr_nans(tmp_subj_rt, tmp_subj_stim)
#                                 tmp_subj_rt, tmp_subj_response = filter_corr_nans(tmp_subj_rt, tmp_subj_response)
                                
#                                 for i, item in enumerate(tmp_subj_response):
#                                         if item == 0:
#                                             tmp_subj_response[i] = -1
                                
#                                 if len(tmp_subj_response) != len(tmp_subj_rt) or len(tmp_subj_stim) != len(tmp_subj_rt):
#                                     print('ERROR')
#                                 else:
#                                     subj_data_counter.append( len(tmp_subj_rt) )
                                
                                
#                                 data_single_subject_acc = ( sum(tmp_subj_response) / len(tmp_subj_response) ) *100
                                

#                                 hits   = 0
#                                 for i, item in enumerate(tmp_subj_response):
#                                     if item == 1:
#                                         if tmp_subj_stim[i] == 'congruent':
#                                             hits += 1

#                                 misses   = 0
#                                 for i, item in enumerate(tmp_subj_response):
#                                     if item == 0:
#                                         if tmp_subj_stim[i] == 'congruent':
#                                             misses += 1
                                            
#                                 fas   = 0
#                                 for i, item in enumerate(tmp_subj_response):
#                                     if item == 0:
#                                         if tmp_subj_stim[i] == 'incongruent':
#                                             hits += 1
                                            
#                                 crs   = 0
#                                 for i, item in enumerate(tmp_subj_response):
#                                     if item == 1:
#                                         if tmp_subj_stim[i] == 'incongruent':
#                                             crs += 1
                                
#                                 #tmp_dprime = calc_dprime(hits, misses, fas, crs)
#                                 #data_single_subject_dprime = tmp_dprime['d']
#                                 data_single_subject_dprime = calc_simple_dprime(hits, fas, len(tmp_subj_response))
                                
                                
#                                 # TO DO : Append pd Dataframe by "features"
#                                 data_single_subject = pd.DataFrame(
#                                     {
#                                         "rt": tmp_subj_rt,
#                                         "response": tmp_subj_response,
#                                         "participant_id": subjID,
#                                         'stim': tmp_subj_stim,
#                                         'accuracy': data_single_subject_acc,
#                                         'dprime': data_single_subject_dprime,
#                                         "Session": session,
#                                     }
#                                     )
                                
#                                 # Append single subject data to cohort dataset 
#                                 dataset_csv.append(data_single_subject)
                                
#                                 subj_counter = subj_counter + 1
                
#                     else:
#                         NO_beh_data_available.append( participants[iSubj] )
#                 else:
#                         NO_beh_data_available.append( participants[iSubj] )

        
#         # Make single dataframe out of subject-wise datasets
#         dataset = pd.concat(dataset_csv)
        
#         print( (f'\n# Got dataset' +
#                 f'\n      -> Columns      : {list(dataset.head(0))}' +
#                 f'\n      -> Subjects     : {subj_counter}' +
#                 f'\n      -> ' + r'$\SS-Input_max' + f' : {np.max(subj_data_counter)}' +
#                 f'\n      -> ' + r'$\SS-Input_min' + f' : {np.min(subj_data_counter)}' +
#                 f'\n      -> ' + r'$\SS-Input_mean' + f' : {np.mean(subj_data_counter)}'
#                 )
#             )
#         print( (f'\n# Data missing from' +
#                 f'\n      -> Subjects     : {np.unique( NO_beh_data_available )}')
#             )


#         # Save the dataset information before running in model
#         filepath = Path(fname_results)  
#         filepath.parent.mkdir(parents=True, exist_ok=True)  
#         dataset.to_csv(filepath)  
            
#     else:
#         dataset = pd.read_csv(fname_results, encoding='utf-8')
    
#     if len( np.unique( dataset['Session'] )) <= 1:
        
#         tmp_dataset_fake = copy.deepcopy(dataset)
        
#         tmp_dataset_fake['rt']       = tmp_dataset_fake['rt'] * 0.8
#         tmp_dataset_fake['response'] = np.ones(len(tmp_dataset_fake['response']))
#         tmp_dataset_fake["Session"]  = 'ses-2' 
        
#         frames = [dataset, tmp_dataset_fake]
#         dataset = pd.concat(frames)
        
    
    
#     dataset_total = copy.deepcopy(dataset)
#     dataset_total["stim"] = 'total'
    
#     frames = [dataset, dataset_total]
#     dataset = pd.concat(frames)
    
    
#     # using tuple unpacking for multiple Axes
#     fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)
#     fig.set_dpi(300)
#     fig.set_facecolor('lightgrey')
#     fig.set_size_inches(12,6)
    
#     # Do LINE-plot style
#     ax1.set_title('Response Time - Group average')
#     sns.lineplot(data=dataset,
#         x="Session", y="rt", hue="stim",estimator="mean", 
#         palette="flare", legend='brief',
#         markers=True, dashes=False, err_style="bars", errorbar=("se", 2)
#         , ax=ax1)
#     legend = ax1.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
    
#     ax2.set_title('Accuracy - Group average')
#     sns.lineplot(data=dataset,
#         x="Session", y="accuracy", hue="stim",estimator="mean", 
#         palette="flare", legend=False,
#         markers=True, dashes=False, err_style="bars", errorbar=("se", 2)
#         , ax=ax2)
#     ax2.hlines(50, ax1.get_xlim()[0],  ax1.get_xlim()[1],
#             linestyles='dotted', colors='k', label='Chance level')
#     ax2.set_ylim([0,100])
#     legend = ax2.legend(frameon = 1)
#     frame = legend.get_frame()
#     frame.set_color('white')
#     frame.set_facecolor('white')
#     frame.set_linewidth(0)
    
#     ax3.set_title('Sensitivity (dprime) - Group average')
#     sns.lineplot(data=dataset,
#         x="Session", y="dprime", hue="stim",estimator="mean", 
#         palette="flare", legend=False,
#         markers=True, dashes=False, err_style="bars", errorbar=("se", 2)
#         , ax=ax3)
#     ax3.hlines(50, ax1.get_xlim()[0],  ax1.get_xlim()[1],
#             linestyles='dotted', colors='k', label='Chance level')
#     ax3.set_ylim([0,10])
#     fname_plot1 = fdir_results + f'/Group-average_N-{subj_counter}_features.png'
#     fig.savefig(fname_plot1)
    
    
#     # Do RIDGE-plot style
#     # RT
#     fig = plt.figure(figsize=[12,6], dpi=300, facecolor='lightgrey')
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     # Initialize the FacetGrid object
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataset, row="Session", hue='stim',
#                             aspect=9, height=1.5, legend_out=True)
#     fig_obj.map_dataframe(sns.kdeplot, x="rt", fill=False, alpha=1)
#     #fig_obj.map_dataframe(sns.kdeplot, x="rt", color='black')
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Reponse time by Sessions', y=0.98)
#     plt.legend(facecolor="white", frameon=False, loc='lower right')
#     fname_plot2 = fdir_results + f'/Group-density_feature-RT_N-{subj_counter}_features.png'
#     fig_obj.savefig(fname_plot2)
    
    
#     # Sensitivity
#     fig = plt.figure(figsize=[12,6], dpi=300, facecolor='lightgrey')
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     # Initialize the FacetGrid object
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataset, row="Session", hue='stim',
#                             aspect=9, height=1.5,
#                             legend_out=False)
#     fig_obj.map_dataframe(sns.kdeplot, x="dprime", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('D prime by Sessions', y=0.98)
#     plt.legend(facecolor="white", frameon=False, loc='lower right')
#     fname_plot3 = fdir_results + f'/Group-density_feature-dprim_N-{subj_counter}_features.png'
#     fig_obj.savefig(fname_plot3)
    
#     # Accuracy
#     fig = plt.figure(figsize=[12,6], dpi=300, facecolor='lightgrey')
#     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
#     # Initialize the FacetGrid object
#     pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
#     fig_obj = sns.FacetGrid(dataset, row="Session", hue='stim',
#                             aspect=9, height=1.5,
#                             legend_out=False)
#     fig_obj.map_dataframe(sns.kdeplot, x="accuracy", fill=False, alpha=1)
#     # Set the subplots to overlap
#     fig_obj.figure.subplots_adjust(hspace=-.25)
#     # Remove axes details that don't play well with overlap
#     fig_obj.set(yticks=[], ylabel="")
#     fig_obj.despine(bottom=True, left=True)
#     plt.suptitle('Accuracy by Sessions', y=0.98)
#     plt.legend(facecolor="white", frameon=False, loc='lower right')
#     fname_plot4 = fdir_results + f'/Group-density_feature-accuracy_N-{subj_counter}_features.png'
#     fig_obj.savefig(fname_plot4)


    
#     return


# def filter_corr_nans(x, y):
    
#     yn = []
#     for n, ix in enumerate(x):
#         if not math.isnan(ix):
#             yn.append(y[n])
    
#     xn       = [ix for ix in x if str(ix) != 'nan']

#     return xn, yn


# def print_label_ridge(x, color, label):
#     # Define and use a simple function to label the plot in axes coordinates
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color='k',
#             ha="left", va="center", transform=ax.transAxes)
    
#     return


# def calc_dprime(hits, misses, fas, crs):
#     """A central component of Signal Detection Theory is d a measure
#     of the ability to discriminate a signal from noise. The d is flanked
#     by the parameters “beta” and c, which are measures of the criterion 
#     that the observer uses to discriminate between the two.
#     These measures can be calculated in every experiment where there is 
#     a signal (e.g. target trials) and noise (e.g. nontarget trials), and 
#     the observer (e.g. subject) is to indicate whether the signal was 
#     present or not. 
#     d, beta and c are statistical measures in a model where noise and
#     noise+signal are regarded as two “probability-of-detection” 
#     distributions on a “threshold-of-detection” continuum. 
#     -> d is basically a Z-score on how well the observer 
#     discriminates the two distributions, i.e. the number of standard 
#     deviations between the probability-of-response distributions for 
#     signal and noise for this given subject.
#     -> c (the criterion) is the number of standard deviations from
#     the midpoint between these two distributions, i.e. a measure 
#     on a continuum from “conservative” to “liberal”."""
#     Z = norm.ppf
#     # Floors an ceilings are replaced by half hits and half FA's
#     half_hit = 0.5 / (hits + misses)
#     half_fa = 0.5 / (fas + crs)

#     # Calculate hit_rate and avoid d' infinity
#     hit_rate = hits / (hits + misses)
#     if hit_rate == 1: 
#         hit_rate = 1 - half_hit
#     if hit_rate == 0: 
#         hit_rate = half_hit

#     # Calculate false alarm rate and avoid d' infinity
#     fa_rate = fas / (fas + crs)
#     if fa_rate == 1: 
#         fa_rate = 1 - half_fa
#     if fa_rate == 0: 
#         fa_rate = half_fa

#     # Return d', beta, c and Ad'
#     out = {}
#     out['d'] = Z(hit_rate) - Z(fa_rate)
#     out['beta'] = math.exp((Z(fa_rate)**2 - Z(hit_rate)**2) / 2)
#     out['c'] = -(Z(hit_rate) + Z(fa_rate)) / 2
#     out['Ad'] = norm.cdf(out['d'] / math.sqrt(2))
    
#     return(out)

# def calc_simple_dprime(hits, fas, total):

    # hitP = hits / total
    # faP  =  fas / total

    # # z-scores
    # hitZ = stats.norm.ppf(hitP)
    # faZ  = stats.norm.ppf(faP)

    # # d-prime
    # dPrime = hitZ-faZ
    # print(dPrime)
    
    # return(dPrime)