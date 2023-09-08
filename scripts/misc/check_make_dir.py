#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)

import os

def check_make_dir(fdir, overwrite=0):
    if isExist := os.path.exists(fdir):
        if overwrite == 1:
            print("Directory exists.. Overwriting it.")
            os.makedirs(fdir)

    else:
        # Create a new directory because it does not exist
        os.makedirs(fdir)
        print("The new directory is created!")
    return