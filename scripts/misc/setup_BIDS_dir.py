#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
import os

def setup_BIDS_dir(subjID=None, ses=None, dir_analysis=None):
#
#
#
#
#
    dir_subjBIDS = dir_analysis + '/data/' + subjID[1] + '/' + ses + '/'

    isExist = os.path.exists(dir_subjBIDS)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS)
        print("New BIDS sub-directory is created!")
    
    isExist = os.path.exists((dir_subjBIDS + 'func'))
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs((dir_subjBIDS + 'func'))
        print("New BIDS sub-directory is created!")
    
    isExist = os.path.exists(dir_subjBIDS + 'anat')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'anat')
        print("New BIDS sub-directory is created!")
    
    isExist = os.path.exists(dir_subjBIDS + 'fmap')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'fmap')
        print("New BIDS sub-directory is created!")
    
    isExist = os.path.exists(dir_subjBIDS + 'dwi')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'dwi')
        print("New BIDS sub-directory is created!")
        
    isExist = os.path.exists(dir_subjBIDS + 'perf')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'perf')
        print("New BIDS sub-directory is created!")
        
    isExist = os.path.exists(dir_subjBIDS + 'meg')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'meg')
        print("New BIDS sub-directory is created!")
        
    isExist = os.path.exists(dir_subjBIDS + 'eeg')
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(dir_subjBIDS + 'eeg')
        print("New BIDS sub-directory is created!")
        
    return dir_subjBIDS
