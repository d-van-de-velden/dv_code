#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#          Alexander Enge          (XX@YY.de)
#          Anne-Sophie Kieslinger  (XX@YY.de)
#
# License: BSD (3-clause)
# 
# 
#
import getpass
import os
from dv_code.scripts.misc.check_make_dir import check_make_dir
from dv_code.scripts.work_participants import update_participants
import keyring
import owncloud
import tkinter as tk
from tkinter import simpledialog



def download_datashare_dcm(params):
    """Downloads new zip files with raw data from MPCDF DataShare.

    Parameters
    ----------
    fdir_datashare : str
        Path of the raw data starting from the DataShare root, like this:
        https://datashare.mpcdf.mpg.de/apps/files/?dir=<fdir_datashare>. Data
        must be organized like <fdir_datashare>/<session>/<participant>_*.zip.
        
    fdir_analysis : datalad.api.Dataset
        The BIDS dataset. New zip files will be downloaded into a 'sourcedata'
        subdataset, with separate subdirectories for each session like on
        DataShare.
    """

    fdir_analysis = params.get('fdir_analysis')
    fdir_datashare = params.get('fdir_datashare')

    # Create subdataset if it doesn't exist
    fdir_source = check_make_dir(fdir_analysis, 'sourcedata')
    fdir_source = params.get('fdir_sourcedata')
    # Get DataShare login credentials
    ROOT = tk.Tk()
    ROOT.withdraw()
    # the input dialog  
    datashare_user = simpledialog.askstring(title="Get DataShare login credentials", prompt="Username:")
    datashare_pass = simpledialog.askstring(title="Get DataShare login credentials", prompt="Password:", show='*')

    if datashare_pass is None:
        datashare_pass = getpass.getpass()
        keyring.set_password('datashare', datashare_user, datashare_pass)

    # Login to DataShare
    domain = 'https://datashare.mpcdf.mpg.de'
    datashare = owncloud.Client(domain)
    datashare.login(datashare_user, datashare_pass)

    # Create empty list / dict to track new data
    new_raw_files = []
    new_participants_sessions = set()

    # Loop over session folders on DataShare

    print('Identifying sessions and datasets from datashare....')
    datashare_sessions = datashare.list(fdir_datashare)

    for datashare_session in datashare_sessions:

        # Loop over files for the current session
        session = datashare_session.name
        if session[:3] == 'ses':

            fdir_session = f'{fdir_source}/{session}'

            files = datashare.list(datashare_session.path)

            for file in files:

                # Explicity exclude certain file names
                if file.name.startswith('_'):
                    continue

                # Download if it doesn't exist
                local_file = f'{fdir_session}/{file.name}'

                if os.path.exists(local_file) == False:

                    # Download zip file
                    print(f'# Downloading {file.path} to {fdir_session}')

                    check_make_dir(fdir_session)
                    datashare.get_file(file, local_file)

                    # Keep track of new data
                    new_raw_files.append(local_file)
                    participant = file.name.split('_')[0]
                    new_participants_sessions.add((participant, session))


    # In the end update/create participant.tsv 
    update_participants(params)

    return print('[DONE]')



def download_datashare_behavioral(params):
    """Downloads new zip files with raw data from MPCDF DataShare.

    Parameters
    ----------
    fdir_datashare : str
        Path of the raw data starting from the DataShare root, like this:
        https://datashare.mpcdf.mpg.de/apps/files/?dir=<fdir_datashare>. Data
        must be organized like <fdir_datashare>/<session>/<participant>_*.zip.
        
    fdir_analysis : datalad.api.Dataset
        The BIDS dataset. New zip files will be downloaded into a 'sourcedata'
        subdataset, with separate subdirectories for each session like on
        DataShare.
    """

    fdir_analysis  = params.get('fdir_analysis')
    fdir_datashare = params.get('fdir_datashare')
    fdir_datashare = fdir_datashare + 'Be_Exp'

    # Create subdataset if it doesn't exist
    fdir_source = check_make_dir(fdir_analysis, 'sourcedata')
    fdir_source = params.get('fdir_sourcedata')
    # Get DataShare login credentials
    ROOT = tk.Tk()
    ROOT.withdraw()
    # the input dialog  
    datashare_user = simpledialog.askstring(title="Get DataShare login credentials", prompt="Username:")
    datashare_pass = simpledialog.askstring(title="Get DataShare login credentials", prompt="Password:", show='*')

    if datashare_pass is None:
        datashare_pass = getpass.getpass()
        keyring.set_password('datashare', datashare_user, datashare_pass)

    # Login to DataShare
    domain = 'https://datashare.mpcdf.mpg.de'
    datashare = owncloud.Client(domain)
    datashare.login(datashare_user, datashare_pass)

    # Create empty list / dict to track new data
    new_raw_files = []
    new_participants_sessions = set()

    # Loop over session folders on DataShare

    print('Identifying sessions and datasets from datashare....')
    datashare_sessions = datashare.list(fdir_datashare)

    for datashare_session in datashare_sessions:

        # Loop over files for the current session
        session = datashare_session.name
        if session[:3] == 'ses':

            fdir_session = f'{fdir_source}/{session}'

            files = datashare.list(datashare_session.path)

            for file in files:

                # Explicity exclude certain file names
                if file.name.startswith('_'):
                    continue

                tmp_fname_idents = os.path.splitext(os.path.basename(file.name))[0]
                tmp_fname_idents = os.path.splitext(os.path.basename(tmp_fname_idents))[0]

                idents  = tmp_fname_idents.split('_')
                subjID  = idents[0]
                tmp_str_BIDS_session = idents[1].replace("-0", "-")
                studyID = idents[2]
                studyID = studyID + '-Behavioral'
                
                # Download if it doesn't exist
                local_file = f'{fdir_session}/{file.name}'

                if os.path.exists(local_file) == False:

                    # Download zip file
                    print(f'# Downloading {file.path} to {fdir_session}')

                    check_make_dir(fdir_session)
                    datashare.get_file(file, local_file)

                    # Keep track of new data
                    new_raw_files.append(local_file)
                    participant = file.name.split('_')[0]
                    new_participants_sessions.add((participant, session))
                
    # In the end update/create participant.tsv 
    update_participants(params)

    return print('[DONE]')
