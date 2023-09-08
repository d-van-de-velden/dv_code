#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
# 
import os
import shutil
import tarfile
import zipfile
import numpy as np


def unpack_compressed(finput=None, fdir_out=None):
    
    dir_entry = list()
    if os.path.isdir(finput):
        dir_entry = os.listdir(finput)
        dir_entry = np.sort( dir_entry )
        
    if os.path.isfile(finput):
        dir_entry.append( finput )

    for in_item in dir_entry:
        
        try_tarfile = True
        if in_item[-4:] == '.zip':
            try:
                unpack_zipfiles(in_item, fdir_out)
            except:
                print("zip-file was actually a tar-file :-/")
                unpack_tarfiles(in_item, fdir_out)
                try_tarfile = False


        if in_item[-7:] == '.tar.gz' and try_tarfile == True:
            unpack_tarfiles(in_item, fdir_out)

    return print('done unpacking...\n')
    
    
def unpack_tarfiles(in_item, fdir_out):
    print(f'# Unpacking tar-file: {in_item}')
    with tarfile.open(in_item, "r") as f:
        f.extractall(fdir_out)
        tmp_files = f.getnames()
    
    for jtem in tmp_files:
        shutil.move(( fdir_out + jtem ), fdir_out)

    print(fdir_out)
    for ktem in os.listdir(fdir_out):
        print( fdir_out + ktem )
        if os.path.isdir( ( fdir_out + ktem ) ):
            print(fdir_out + ktem)
            shutil.rmtree(( fdir_out + ktem), ignore_errors=True) 
            
    return

def unpack_zipfiles(in_item, fdir_out):
    print(f'# Unpacking zip-file: {in_item}')
    with zipfile.ZipFile(in_item,"r") as zip_ref:
        zip_ref.extractall(fdir_out)
        
    return
                    