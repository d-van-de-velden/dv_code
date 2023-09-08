#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#

def get_only_nifti(fnames=None, compressed=2):
#
# compressed : Find only nifti's (0), 
#              find only packed nifti's (1),
#              find all nifti's (2) [default] 
#
    fnames_nii = list()
    for item in fnames:
        if compressed == 2:
            if item[-4:] == '.nii' or item[-7:] == '.nii.gz':
                fnames_nii.append(item)
        elif compressed == 1:
            if item[-7:] == '.nii.gz':
                fnames_nii.append(item)
        elif compressed == 0:
            if item[-5:] == '.nii':
                fnames_nii.append(item)

    return fnames_nii



def get_only_zip(fnames):   
#
#
#
    fnames_zip = list()
    for item in fnames:
        if item[-4:] == '.zip' and item[-9:] != '_logs.zip':
            fnames_zip.append(item)
            
    return fnames_zip

def get_only_json(fnames):   
#
#
#
    fnames_json = list()
    for item in fnames:
        if item[-5:] == '.json':
            fnames_json.append(item)
    
    return fnames_json