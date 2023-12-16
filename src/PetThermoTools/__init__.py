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
import Thermobar as pt
from os import path

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
