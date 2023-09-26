#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#

from nilearn.glm import threshold_stats_img
import nibabel as nib


def make_overlap_contrastmap_surface(fname_in=None, plotting_title='',
                            alpha=0.05, mcc='fdr', cluster_ext=2, 
                            fname_out="output.png"):
    
    fname_in = ["/data/p_02825/cocoa/data/derivatives/1st_level/Group/ses-1/GroupN20_ses-1_avg_z_unimodal_audios_vs_null_event_tSNR50.nii", 
                "/data/p_02825/cocoa/data/derivatives/1st_level/Group/ses-1/GroupN20_ses-1_avg_z_unimodal_images_vs_null_event_tSNR50.nii",
    ]
    

    if type(fname_in) is list:
        print('Provided multiple filenames:\n' + str(fname_in[0::]))
        vol_data_init   = nib.load(fname_in[0])
    elif type(fname_in) is str:
        print('Provided one filename:\n' + fname_in)
        vol_data_init   = nib.load(fname_in)
    else:
        print('neither a tuple or a list')

    _, threshold = threshold_stats_img(stat_img=vol_data_init, alpha=alpha,
                                        height_control=mcc,
                                        cluster_threshold=cluster_ext)
    print(f"All contrast maps are thresholded by:\n   {mcc} p<{alpha:.3f} threshold: {threshold:.3f}")
    
    img_overlap_mask = vol_data_init.get_fdata() * 0
    for ifname in fname_in:
        
        tmp_vol_data   = nib.load(ifname)

        tmp_vol_data, threshold = threshold_stats_img(stat_img=tmp_vol_data, alpha=alpha,
                                    height_control=mcc,
                                    cluster_threshold=cluster_ext)
        
        img_overlap_mask += tmp_vol_data.get_fdata()
    
    img_overlap_mask[img_overlap_mask > 0] = threshold * 2
    
    img_overlap = nib.Nifti1Image(img_overlap_mask, vol_data_init.affine)

    surface_data_LH, surface_data_RH = transform_to_surface(img_overlap)
    
    viz_surface_plot(surface_data_LH=surface_data_LH,
                    surface_data_RH=surface_data_RH,
                    threshold=threshold, 
                    plotting_title=plotting_title, 
                    fname_out=fname_out)

    
    return