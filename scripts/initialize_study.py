"""
This function is the very first step in setting up an analysis pipeline for a study.
This function will repeatedly ask the researcher questions in an interactive way.
This function creates the analysis_params.json file based on these answers.
--> This file serves as the foundation of folder/file structure.
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden    (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
# 
# 
#
import os 
import json
import tkinter
from tkinter import simpledialog
from dv_code.scripts.misc import check_make_dir
from dv_code.scripts.misc.feedback_messages import say_welcome

def ___init_stud_dirs(analysis_params_dict, arg1, arg2):
    # Asking for and creating some folders 
  result = analysis_params_dict.get('fdir_study') + arg1
  analysis_params_dict[arg2] = (result)
  check_make_dir(analysis_params_dict.get('fdir_analysis'))

  return result


def initialize_study(RunIt=True):
  #
  #
  #
  #
  say_welcome('0X.X.0X')

  if RunIt == False:
    return print('No study initialization taking place.....')
  
  analysis_params_dict = {}
  
  # Fetch directory where all the data/code/etc is stored
  parent = tkinter.Tk()
  parent.withdraw()
  item = simpledialog.askstring("New Item", 
                                ("Hi and welcome." 
                                + "Please provide the filepath to the "
                                + "study you want to work on:"),
                                parent=parent)
  analysis_params_dict["fdir_study"] = (item)


  fname_analysis_params = analysis_params_dict.get('fdir_study') + "analysis_params.json"

  create_analysis_params = False
  if os.path.exists(fname_analysis_params) == True:
    item = simpledialog.askstring("New Item", 
                                  "Study already initiliazed. Want to restart ?", 
                                  parent=parent)
    if item == 'Yes':
      create_analysis_params = True
  else:
    create_analysis_params = False


  if create_analysis_params:

    # Ask for researchers name
    item = simpledialog.askstring("New Item", 
                                  "Welcome to your study! Please enter your name:", 
                                  parent=parent)
    analysis_params_dict["UserID"] = (item)

    # Fetch study ID
    item = simpledialog.askstring("New Item", 
                                  "Please enter the name of your study:", 
                                  parent=parent)
    analysis_params_dict["study_id"] = (item)

    item = analysis_params_dict.get('fdir_study') + "/dv_code/"
    analysis_params_dict["fdir_analysis"] = (item)  

    item = analysis_params_dict.get('fdir_study') + "/sge"
    analysis_params_dict["fdir_sge"] = (item)  

    item = analysis_params_dict.get('fdir_study') + "/bash"
    analysis_params_dict["fdir_bash"] = (item)  

    item = analysis_params_dict.get('fdir_study') + "/dv_code/scripts/lib_bash"
    analysis_params_dict["fdir_lib_bash"] = (item)  

    item = analysis_params_dict.get('fdir_study') + "/dv_code/dcm2bids_config.json"
    analysis_params_dict["fname_config"] = (item)

    item = ___init_stud_dirs(analysis_params_dict, "/data/", "fdir_data")
    item = ___init_stud_dirs(analysis_params_dict, "/data/FREESURFER", "fdir_fs")
    item = ___init_stud_dirs(analysis_params_dict, "/data/sourcedata", "fdir_sourcedata")
    item = analysis_params_dict.get('fdir_study') + "/data/derivatives"
    analysis_params_dict["fdir_proc"] = (item)  

    item = ___init_stud_dirs(
        analysis_params_dict, "/data/derivatives/preprocessed",
        "fdir_proc_pre")
    item = ___init_stud_dirs(analysis_params_dict, 
                            "/data/derivatives/stats", "fdir_proc_stat")
    item = "/COCOA/raw/german/"
    analysis_params_dict["fdir_datashare"] = (item)


    # Serializing json
    analysis_params_dict_json_object = json.dumps(analysis_params_dict, indent=4)

    # Writing to sample.json
    fname_analysis_params = analysis_params_dict.get('fdir_study') + "analysis_params.json"
    with open(fname_analysis_params, "w") as outfile:
      outfile.write(analysis_params_dict_json_object)


  return



