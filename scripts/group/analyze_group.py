#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os
import numpy as np
from seaborn import violinplot
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from dv_code.scripts.misc.select_filetype import get_only_nifti
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.processing.func_MRI_quality_metric import get_tSNR
from dv_code.scripts.viz.combine_plots import viz_combine_plots
from dv_code.scripts.work_participants import get_participants
from seaborn import violinplot, swarmplot

def overview_tSNR(participants=None, params=None, tSNR_threshold=35, N_runsOK=2):
    #
    #
    #    
    if participants is None:
        participants = get_participants(params)
    check_make_dir(params.get('fdir_proc_pre'))

    fdir_fig = params.get('fdir_proc') + '/group/QC/'
    check_make_dir(fdir_fig)      
        
    print(
        f'Providing overview of functional MR data quality in tSNR for: {len(participants)} subjects\n'
    )

    print('  # Collecting all the tSNR data..\n')


    ALL_data_plot = []
    ALL_tSNR_median = []
    ALL_fname_runs = []
    for iSubj in range(len(participants)):
        print(f'    Identifying sessions from "data/{participants[iSubj]}"/ ....')
        tmp_fdir_subj = params.get('fdir_data') + participants[iSubj]

        sessions = os.listdir(tmp_fdir_subj)

        for session in sessions:
            
            tmp_fdir_func = f'{tmp_fdir_subj}/{session}/func/'

            runs = os.listdir(tmp_fdir_func)
            runs = get_only_nifti(runs, 1)
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

            tSNR_median_across_runs_raw     = []
            tSNR_median_across_runs_preproc = []
            
            runs.sort()
            
            for run in runs:
                tmp_fname_funcMR = f'{tmp_fdir_func}/{run}'
                tmp_fname_func = os.path.splitext(os.path.basename(tmp_fname_funcMR))[0]
                fname_func = os.path.splitext(os.path.basename(tmp_fname_func))[0]


                fname_func_tSNR = f'{fdir_derivatives_func}/{fname_func}_tSNR_r.nii.gz'
                if os.path.exists(fname_func_tSNR):
                    _, tSNR_median, _ = get_tSNR(fname_func_tSNR, params)


                fname_func_r_tSNR = f'{fdir_derivatives_func}/{fname_func}_st_mcf_topUP_tSNR_r.nii.gz'
                if os.path.exists(fname_func_r_tSNR):
                    _, tSNR_median1, _ = get_tSNR(fname_func_r_tSNR, params)


                ALL_tSNR_median.extend(
                    (
                        {
                            'Subject': participants[iSubj],
                            'Session': session,
                            'median tSNR (whole)': tSNR_median[0],
                            'median tSNR (gm)': tSNR_median[1],
                            'Processing\nstate': 'raw',
                        },
                        {
                            'Subject': participants[iSubj],
                            'Session': session,
                            'median tSNR (whole)': tSNR_median1[0],
                            'median tSNR (gm)': tSNR_median1[1],
                            'Processing\nstate': 'preproc',
                        },
                    )
                )
                tSNR_median_across_runs_raw.append(tSNR_median[1])
                tSNR_median_across_runs_preproc.append(tSNR_median1[1])
            
            tmp_N_runs = len(tSNR_median_across_runs_raw)
            fig, ax = plt.subplots()
            fig.set_size_inches(8, 3)
            fig.set_dpi(300)
            if ((np.array(tSNR_median_across_runs_preproc) < tSNR_threshold).sum()) < N_runsOK:
                fig.set_facecolor('lightgrey')
            else:
                fig.set_facecolor([1,0.8,0.79])
            fig.set_tight_layout
            
            ax.set_title(f'Subject: {participants[iSubj]} | Session: {session}')
            #create basic scatterplot
            x = np.linspace(1, tmp_N_runs, tmp_N_runs)
            plt.plot(x, tSNR_median_across_runs_raw,
                    'o', color=sns.color_palette("pastel")[0])
            plt.plot(x, tSNR_median_across_runs_preproc,
                    'o', color=sns.color_palette("pastel")[1])

            #obtain m (slope) and b(intercept) of linear regression line           
            coeff = np.polyfit(x, tSNR_median_across_runs_raw, 3)
            yn = np.poly1d(coeff)
            coeff1 = np.polyfit(x, tSNR_median_across_runs_preproc, 3)
            yn1 = np.poly1d(coeff1)
            plt.plot(x, yn(x), 
                    linewidth=4,
                    color=sns.color_palette("pastel")[0])
            plt.plot(x, yn1(x), 
                    linewidth=4,
                    color=sns.color_palette("pastel")[1])
            
            ax.set_xticks(x)
            ax.set_xlim([0.5, 5])
            ax.set_ylim([0, 120])
            
            plot_tSNR_threshold = np.zeros(tmp_N_runs)
            plot_tSNR_threshold1 = np.zeros(tmp_N_runs)
            plot_tSNR_threshold1[:] = tSNR_threshold
            # Fill in area under the curve and the horizontal lines
            plt.axhline(tSNR_threshold, linestyle='--', 
                        color='lightgrey', label="tSNR threshold",
                        alpha=.85)
            plt.fill_between(x=x, 
                            y1=plot_tSNR_threshold, 
                            y2=plot_tSNR_threshold1, 
                            color='red',  interpolate=True, alpha=.15)



            ax_X_labels = []
            for isesLabel in np.arange( tmp_N_runs ):
                ax_X_labels.append('run ' + str(isesLabel+1))
            ax.set_xticklabels(ax_X_labels)
    
            fname_fig = f'overview_tSNR_runs_{participants[iSubj]}_{session}.png'
            
            fig.savefig((fdir_fig + fname_fig), dpi=300)
            ALL_fname_runs.append((fdir_fig + fname_fig))

    viz_combine_plots(ALL_fname_runs, 'pdf', params)

    ALL_data_plot = pd.DataFrame(ALL_tSNR_median)
    fname_data = fdir_fig + 'overview_tSNR.csv'
    ALL_data_plot.to_csv(fname_data, sep=',')

    
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 3)
    fig.set_dpi(300)
    fig.set_facecolor('lightgrey')
    fig.set_tight_layout
    fig.suptitle('Group overview - tSNR', fontsize=16)

    ax = violinplot(
        ALL_data_plot, x='Session', y='median tSNR (gm)', hue='Processing\nstate',
        linewidth=0.5,
        palette=sns.color_palette("pastel")[:2],
        scale_hue=True, cut=0
    )
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax_X_labels = []
    for isesLabel in np.arange( len(set(list(ALL_data_plot.get('Session')))) ):
        ax_X_labels.append('ses-0' + str(isesLabel+1))
    ax.set_xticklabels(ax_X_labels)
    
    plt.axhline(tSNR_threshold, linestyle='--', 
                        color='red', label="tSNR threshold",
                        alpha=.85)
    x = np.linspace(-1, tmp_N_runs+1, tmp_N_runs+2)
    plot_tSNR_threshold = np.zeros(tmp_N_runs+2)
    plot_tSNR_threshold1 = np.zeros(tmp_N_runs+2)
    plot_tSNR_threshold1[:] = tSNR_threshold
    plt.fill_between(x=x,
                    y1=plot_tSNR_threshold,
                    y2=plot_tSNR_threshold1,
                    color='red',  interpolate=True, alpha=.15)
    ax.set_xlim([-1, 4])
    ax.set_ylim([0, 140])

    tMp_dat = ALL_data_plot
    for isesLabel in np.arange( len(set(list(ALL_data_plot.get('Session')))) ):

        a_raw      = tMp_dat.loc[tMp_dat["Processing\nstate"] =='raw' ]
        tot_a_raw  = a_raw.shape[0]
        por_a_raw  = a_raw.loc[a_raw["median tSNR (gm)"] > tSNR_threshold].shape[0]
        perc_a_raw = round(por_a_raw / tot_a_raw * 100)
        str_raw    = f'{perc_a_raw}%'
        annot_raw  = plt.annotate(str_raw, xy=(isesLabel-0.25, tSNR_threshold+5))
        annot_raw.set_fontsize(6)
        
        a_preproc      = tMp_dat.loc[tMp_dat["Processing\nstate"] =='preproc' ]
        tot_a_preproc  = a_preproc.shape[0]
        por_a_preproc  = a_preproc.loc[a_preproc["median tSNR (gm)"] > tSNR_threshold].shape[0]
        perc_a_preproc =  round(por_a_preproc / tot_a_preproc * 100)
        str_preproc    = f'{perc_a_preproc}%'
        annot_preproc  = plt.annotate(str_preproc, xy=(isesLabel+0.25, tSNR_threshold+5))
        annot_preproc.set_fontsize(6)
    
    fname_fig = 'overview_tSNR.png'
    fig.savefig((fdir_fig + fname_fig), dpi=300)
    
    
    
    ALL_data_plot = ALL_data_plot.loc[ALL_data_plot["Processing\nstate"] =='preproc' ]
    
    fig, axis = plt.subplots(len(set(list(ALL_data_plot.get('Session'))))+1, 1)
    fig.set_size_inches(8, 4*len(set(list(ALL_data_plot.get('Session')))))
    fig.set_dpi(300)
    fig.set_facecolor('lightgrey')
    fig.set_tight_layout('tight')
    fig.suptitle('Group overview - tSNR', fontsize=16)
    
    for isesLabel in np.arange( len(set(list(ALL_data_plot.get('Session')))) ):
        ax_idx = int(isesLabel)

        tmp_ses_ALLdata = ALL_data_plot.loc[ALL_data_plot["Session"] ==f'ses-{isesLabel+1}' ]
        
        
        subj_list = np.unique( list(ALL_data_plot.get('Subject')))
        use_idx = np.zeros([len(subj_list)])
        for count, iSubj in enumerate( subj_list ):
            
            tmp_Data_subj = tmp_ses_ALLdata.loc[tmp_ses_ALLdata["Subject"] == iSubj ]
            tmp_tSNR_subj = tmp_Data_subj["median tSNR (gm)"]
            
            if ((np.array(tmp_tSNR_subj) < tSNR_threshold).sum()) < N_runsOK:
                use_idx[count] = 1
            
        # all datasets each session
        if isesLabel == 1:
            swarmplot(tmp_ses_ALLdata, x='Subject', y='median tSNR (gm)', 
                hue='Subject', linewidth=0.5, palette='deep',
                ax=axis[ax_idx], legend="full")
        else:
            swarmplot(tmp_ses_ALLdata, x='Subject', y='median tSNR (gm)', 
                hue='Subject', linewidth=0.5, palette='deep',
                ax=axis[ax_idx], legend=False)
        
        axis[ax_idx].spines['left'].set_color('white')
        axis[ax_idx].spines['right'].set_color('white')
        axis[ax_idx].spines['top'].set_color('white')
        axis[ax_idx].spines['bottom'].set_color('white')
        axis_X_labels = np.unique( list(tmp_ses_ALLdata.get('Subject')) )

        axis[ax_idx].set_xticklabels(axis_X_labels)
        axis[ax_idx].tick_params(axis='x', labelrotation = 90)
        axis[ax_idx].set_title(f'Session {isesLabel+1}  (N={len(subj_list)})')
        
        if isesLabel == 1:
            axis[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        axis[ax_idx].axhline(tSNR_threshold, linestyle='--', 
                            color='red', label="tSNR threshold",
                            alpha=.85)
        N_subs = len(set(list(tmp_ses_ALLdata.get('Subject'))))
        x = np.linspace(-1, N_subs+2, N_subs+2)
        plot_tSNR_threshold = np.zeros(N_subs+2)
        plot_tSNR_threshold1 = np.zeros(N_subs+2)
        plot_tSNR_threshold1[:] = tSNR_threshold
        axis[ax_idx].fill_between(x=x,
                        y1=plot_tSNR_threshold,
                        y2=plot_tSNR_threshold1,
                        color='red',  interpolate=True, alpha=.15)
        
        axis[ax_idx].fill_between(x=np.linspace(0, N_subs, N_subs),
                        y1=np.zeros(N_subs),
                        y2=use_idx,
                        color='red', alpha=.15)
        
        
        axis[ax_idx].set_xlim([-1, N_subs+2])
        axis[ax_idx].set_ylim([0, 120])
    fig.set_tight_layout('tight')
    
    fname_fig = 'overview_tSNR1.png'
    fig.savefig((fdir_fig + fname_fig), dpi=300)

    return
