__author__ = 'Matthew Gleeson'

import numpy as np
import pandas as pd
import sys
import os
import warnings
import site
import sysconfig
import importlib
from pathlib import Path
import matplotlib.pyplot as plt
# import subprocess
# import glob
from subprocess import Popen, PIPE
# try:
#     import Thermobar as pt
# except:
#     print('Thermobar import not successful. This is likely due to your numpy version. \n If you are using numpy version > 2 this is currently incompatible with Thermobar.')

from os import path

# # Define the path to your custom environment
# from pathlib import Path
# JULIA_ENV_PATH = Path.home() / ".petthermotools_julia_env"

# if JULIA_ENV_PATH.exists():
#     import juliapkg
    
#     # Lock the Julia version to the 1.11 series
#     # This is lightweight and doesn't download anything yet
#     juliapkg.require_julia("~1.11")

# try:
if os.path.exists('/home/jovyan/shared/Models/alphaMELTS'):
    sys.path.append('/home/jovyan/shared/Models/alphaMELTS')
    # from meltsdynamic import MELTSdynamic
else:
    site_packages_dirs = site.getsitepackages()

    # for dir in site_packages_dirs:
    #     # print(f"Checking {dir} for .pth files")
    #     for file in os.listdir(dir):
    #         if file.endswith(".pth"):
    #             print(f"Found .pth file: {file}")

    # Path to your site-packages directory (adjust as needed)
    site_packages_path = site.getsitepackages()

    fail = True
    test = True
    for i in range(len(site_packages_path)):
        # Path to the .pth file
        pth_file_path = os.path.join(site_packages_path[i], "my_MELTS_path.pth")

        # Check file location the .pth file
        if os.path.exists(pth_file_path):
            print("alphaMELTS for Python files located.")
        else:
            if test:
                print('alphaMELTS for Python files not automatically located on the Python path. \nPlease make sure you have appended the alphaMELTS for Python files to the Python path \nimport sys \nsys.path.append(r"insert_your_path_to_melts_here") \nYou are looking for the location of the meltsdynamic.py file.')
                test = False

# test installation script
from petthermotools.Installation import *
# General functions that are used across multiple types of calculations
from petthermotools.GenFuncs import *
# Plotting functions that are used in multiple types of calculations
from petthermotools.Plotting import *
# This has the main functions that call the barometry calculations
from petthermotools.Barom import *
# This contains the functions required to find liquidis
from petthermotools.Liq import *
# Code to create phase diagrams
from petthermotools.PhaseDiagrams import *
# Functions to carry out specific crystallisation calculations. Simply calls the code in petthermotools.Path
from petthermotools.Path_wrappers import *
# This contains the functions required to perform path calculations
from petthermotools.Path import *
# This contains the functions required to simulate mantle melting
from petthermotools.Melting import *
# This contains the functions required to perform volatile saturation calculations
from petthermotools.Saturation import *
# This contains the functions used to call MELTS calculations
from petthermotools.MELTS import *
# Compsitions saved for use in petthermotools
from petthermotools.Compositions import *
# This contains the functions used to perform calculations using the Holland et al. 2018 thermodynamic dataset - currently disabled
# try:
#     from petthermotools.Holland import *
# except:
#     print('Warning: MAGEMin calculations cannot be performed')

# version
from ._version import __version__
