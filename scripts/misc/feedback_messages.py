#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Daniel van de Velden (d.vandevelden@yahoo.de)
#
# License: BSD (3-clause)
#
#
#
import warnings

def error_message(input_str):
#   Simple error structure for the toolbox
    raise TypeError(f"DvdV argues: {input_str}")


def warn_message(input_str):
##   Simple warning structure for the toolbox
    return print(f"DvdV warns: {input_str}")

def say_welcome(vers_info='1.52.3'):
    """ Provides user with welcome image

    Parameters
    ----------
    version : str
    Current version of SkeideLab Script repository.

    Returns
    -------
    [] : String Output
    """
    print(
    '       ____\n'+
    '     /     \  __               ___    __\n'+
    '    /       ||  |             /__/   |  |\n'+
    '   /    ___/ |  | __  ______  __     |  | ______\n'+
    '   \   |___  |  |/ / |      ||  | ___|  ||      |\n'+
    '    \___   \ |    /  |  |_| ||  ||      ||  |_| |\n'+
    '       /   \ |    \  |   ___||  || |_|  ||   ___|\n'+
    '    __/    / |  |\ \ |  |___ |  ||      ||  |___\n'+
    '  /_______/  |__| \_\|______||__||______||______|\n'+
    '    ___                   __\n'+
    '   |   |                 |  |\n'+
    '   |   |                 |  |\n'+
    '   |   |                 |  |\n'+
    '   |   |              _  |  |___\n'+
    '   |   |         ____/ \ |      \ \n'+
    '   |   |_____   /      | |  |_| |\n'+
    '   |         \ |  |_|  | |      |\n'+
    '    \________/  \___/|_/ \______/  version: ' + vers_info
    'Written by: Daniel van de Velden                         '
    )
    return


"""
TO-DO: Create class structure for those functionalities, 
before their introduction to the toolbox takes place.
"""