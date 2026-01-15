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
from os import path

# import importlib.util

# # Check if 'meltsdynamic' exists in the current Python path
# spec = importlib.util.find_spec("meltsdynamic")

# if spec is not None:
#     print("Module found!")
# else:
#     print("Module not found on the current path.")

try:
    if os.path.exists('/home/jovyan/shared/Models/alphaMELTS'):
        sys.path.append('/home/jovyan/shared/Models/alphaMELTS')
    from meltsdynamic import MELTSdynamic
    print("alphaMELTS for Python files successfully located.")
except (ImportError, ModuleNotFoundError) as e:
    warnings.warn(f"Failed to import MELTSdynamic: {e}. \nUse `import sys; sys.path.append(r'path_to_alphaMELTS_4_python')` to add the alphaMELTS for Python files to the Python path.")

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
