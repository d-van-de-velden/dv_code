#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os

from matplotlib import pyplot as plt

from nilearn import surface

import nilearn
from nilearn import plotting

from PIL import Image

def transform_to_surface(data_in=None):

    fsaverage    = nilearn.datasets.fetch_surf_fsaverage()
    surface_data_RH = surface.vol_to_surf(data_in, fsaverage.pial_right, 
                                    interpolation='nearest', # interpolation : {'linear', 'nearest'}
                                    radius=3, # mm) of the neighbourhood from which samples are drawn
                                    )
    surface_data_LH = surface.vol_to_surf(data_in, fsaverage.pial_left, 
                                    interpolation='nearest', # interpolation : {'linear', 'nearest'}
                                    radius=3, # mm) of the neighbourhood from which samples are drawn
                                    )
    
    return surface_data_LH, surface_data_RH



def viz_surface_plot(surface_data_LH=None, surface_data_RH=None,
                    threshold=0, plotting_title='', fname_out="output.png"):

    if surface_data_LH.shape == surface_data_RH.shape:
        if len( surface_data_LH.shape ) == 1:
            print('   Data is non time-resolved.')
            n_scans = 1
            n_vert  = surface_data_LH.shape[0]
        else:
            print('   Data is time-resolved.')
            n_scans = surface_data_LH.shape[1]
            n_vert  = surface_data_LH.shape[0]

        print('   N scans    : ' + str(n_scans))
        print('   N vertices : ' + str(n_vert) + ' per hemisphere')
    else:
        print('   Left and Right hemisphere are not equal.')
        
    fig_size_y = 5
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4,
                                            sharey=True, sharex=True,
                                            figsize=(fig_size_y*4,fig_size_y),
                                            subplot_kw={'projection': '3d'})
    
    fsaverage    = nilearn.datasets.fetch_surf_fsaverage()
    
    plotting.plot_surf_stat_map(
        fsaverage.infl_right, surface_data_RH, 
        hemi='right', view='lateral',
        threshold=threshold, colorbar=False,
        bg_map=fsaverage.sulc_right, 
        figure=fig , axes=ax1)
    
    plotting.plot_surf_stat_map(
        fsaverage.infl_right, surface_data_RH, 
        hemi='right', view='medial',
        threshold=threshold, colorbar=False,
        bg_map=fsaverage.sulc_right, 
        figure=fig , axes=ax2)
        
    plotting.plot_surf_stat_map(
        fsaverage.infl_left, surface_data_LH, 
        hemi='left', view='medial',
        threshold=threshold, colorbar=False,
        bg_map=fsaverage.sulc_left, 
        figure=fig , axes=ax3)
        
    plotting.plot_surf_stat_map(
        fsaverage.infl_left, surface_data_LH, 
        hemi='left', view='lateral',
        threshold=threshold, colorbar=False,
        bg_map=fsaverage.sulc_left, 
        figure=fig , axes=ax4)
    
    fig.tight_layout()
    fig.suptitle(plotting_title, fontsize=20)
    
    fdir_output = os.path.dirname(fname_out)
    fname_dummy_out_A = fdir_output + 'dummy_A.png'
    plt.savefig(fname_dummy_out_A, dpi=450)
    
    fig2, (ax5) = plt.subplots(1, 1,
                                sharey=True, sharex=True,
                                figsize=(fig_size_y,fig_size_y),
                                subplot_kw={'projection': '3d'})
    plotting.plot_surf_stat_map(
        fsaverage.infl_left, surface_data_LH, 
        hemi='left', view='lateral',
        threshold=threshold, colorbar=True,
        bg_map=fsaverage.sulc_left, 
        figure=fig , axes=ax5)
    
    fname_dummy_out_B = fdir_output + 'dummy_B.png'
    plt.savefig(fname_dummy_out_B, dpi=450)
    
    im = Image.open(fname_dummy_out_B)

    (cord_x, cord_y, width, height) = 1750, 0, im.size[0], im.size[1]
    crop_rectangle = (cord_x, cord_y, width, height)
    cropped_im = im.crop(crop_rectangle)

    fname_dummy_out_C = fdir_output + 'dummy_C.png'
    cropped_im.save(fname_dummy_out_C)
    
    images = [Image.open(x) for x in [fname_dummy_out_A, fname_dummy_out_C]]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]

    new_im.save(fname_out)
    
    
    return