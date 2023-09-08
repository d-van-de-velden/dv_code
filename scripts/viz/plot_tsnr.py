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
import matplotlib.pyplot as plt
from dv_code.scripts.misc import check_make_dir




def plot_tsnr(fname_tsnr=None, data=None, ROI_label=None, dat_idx=None):
    
    fdir_qc = f'{os.path.dirname(fname_tsnr)}/QC/'
    check_make_dir(fdir_qc)
    
    tmp_fname_tsnr = os.path.splitext(os.path.basename(fname_tsnr))[0]
    tmp_fname_tsnr = os.path.splitext(os.path.basename(tmp_fname_tsnr))[0]

    fname_qc = f'QC_tSNR_{tmp_fname_tsnr}.png'
    
    fname_plot = (fdir_qc + fname_qc)
    
    
    length_data = len(data)
    
    if length_data != len(ROI_label):
        print('ERROR')
    
    figure, axis = plt.subplots(1,3)
    figure.set_size_inches(10, 5)
    figure.set_dpi(100)
    figure.set_facecolor('grey')
    
    # Decide on the data to be plotted
    idat = data[dat_idx]
    
    tra_view = idat[:,:,int(idat.shape[2]/2)]
    axis[0].imshow(tra_view, cmap='gray')

    fro_view = idat[:,int(idat.shape[1]/2), :]
    axis[1].imshow( fro_view , cmap='gray')
    axis[1].set_title(f'File: {os.path.basename(fname_tsnr)}')

    sag_view = idat[int(idat.shape[0]/2), :, :]
    axis[2].imshow(sag_view, cmap='gray')

    y_range = np.arange(-round(length_data/2),round(length_data/2) )
    for i, jdat in enumerate(data):

        pos_x = sag_view.shape[1] * 1.05
        pos_y = sag_view.shape[0] / 2   + (y_range[i] * 50)
        
        s = ( r'$\widetilde{tSNR}$' + ROI_label[i] + ' = ' + str(round( np.nanmedian( jdat ), 2)) )
        
        if i == dat_idx:
            plt.text(pos_x, pos_y, s, bbox=dict(fill=False, edgecolor='red', linewidth=2), weight='bold')
        else:
            plt.text(pos_x, pos_y, s, bbox=dict(fill=False, edgecolor='red', linewidth=2))
        

    figure.set_tight_layout('rect')
    
    plt.savefig(fname_plot)
    
    return