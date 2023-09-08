
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#

from PIL import Image
import os
from dv_code.scripts.misc.analysis_feedback import loadingBar

def viz_concat_h(imgN):
    
    
    img0 = imgN[0]
    init_width = img0.width

    tot_width = 0
    for img in imgN:
        tot_width = tot_width + img.width

    dst = Image.new('RGB', (tot_width, img0.height))

    width_counter = 0
    for img in imgN:
        dst.paste(img, (width_counter, 0))
        width_counter = width_counter + img.width

    return dst

def viz_concat_v(imgN):
    
    img0 = imgN[0]
    init_width = img0.width

    tot_height = 0
    for img in imgN:
        tot_height = tot_height + img.height

    dst = Image.new('RGB', (init_width, tot_height))

    height_counter = 0
    for img in imgN:
        dst.paste(img, (0, height_counter))
        print(height_counter)
        height_counter = height_counter + img.height


    return dst


def viz_combine_plots(fnames, format='pdf', params=None):
    
    fdir_output = params.get('fdir_proc') + '/group/QC/'
    fname_output = 'tSNR_acrossruns_ins_subject'
    iml = []
    print(f"{fnames=}")
    fnames = sorted(fnames)
    for i, fname in enumerate(fnames):
        loadingBar(i, len(fnames), (fname + '\n'))
        imgs = Image.open(fname)
        rgb_im = imgs.convert('RGB') # to prevent errors: cannot save mode RGBA
        iml.append(rgb_im)
    
    print(iml)
    image = iml[0]
    iml.pop(0)

    if format == 'pdf':
        fend_output = ".pdf"
        fname_output = fdir_output + fname_output + fend_output
        image.save(fname_output, "PDF" , resolution=100.0, save_all=True, append_images=iml)
        os.system(fname_output)
    
    if format == 'png':
        fend_output = ".png"
        fname_output = fdir_output + fname_output + fend_output
        image.save(fname_output, dpi=300)
        

    return