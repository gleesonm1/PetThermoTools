__author__ = 'Matthew Gleeson'

import numpy as np
import pandas as pd
import sys
import os
import pickle
import matplotlib.pyplot as plt
# import subprocess
# import glob
from subprocess import Popen, PIPE
try:
    import Thermobar as pt
except:
    print('Thermobar import not successful. This is likely due to your numpy version. \n If you are using numpy version > 2 this is currently incompatible with Thermobar.')

from os import path

try:
    if os.path.exists('/home/jovyan/shared/Models/alphaMELTS'):
        sys.path.append('/home/jovyan/shared/Models/alphaMELTS')
    from meltsdynamic import MELTSdynamic
except:
    print('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

# test installation script
from PetThermoTools.Installation import *
# General functions that are used across multiple types of calculations
from PetThermoTools.GenFuncs import *
# Plotting functions that are used in multiple types of calculations
from PetThermoTools.Plotting import *
# This has the main functions that call the barometry calculations
from PetThermoTools.Barom import *
# This contains the functions required to find liquidis
from PetThermoTools.Liq import *
# Code to create phase diagrams
from PetThermoTools.PhaseDiagrams import *
# Functions to carry out specific crystallisation calculations. Simply calls the code in PetThermoTools.Path
from PetThermoTools.Path_wrappers import *
# This contains the functions required to perform path calculations
from PetThermoTools.Path import *
# This contains the functions required to simulate mantle melting
from PetThermoTools.Melting import *
# This contains the functions required to perform volatile saturation calculations
from PetThermoTools.Saturation import *
# This contains the functions used to call MELTS calculations
from PetThermoTools.MELTS import *
# Compsitions saved for use in PetThermoTools
from PetThermoTools.Compositions import *
# This contains the functions used to perform calculations using the Holland et al. 2018 thermodynamic dataset - currently disabled
# try:
#     from PetThermoTools.Holland import *
# except:
#     print('Warning: MAGEMin calculations cannot be performed')

# version
from ._version import __version__
