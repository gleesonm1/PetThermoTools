__author__ = 'Matthew Gleeson'

import numpy as np
import pandas as pd
import sys
import os
import pickle
import matplotlib.pyplot as plt
import subprocess
import glob
from subprocess import Popen, PIPE
import Thermobar as pt
from os import path

# General functions that are used across multiple types of calculations
from pyMELTScalc.GenFuncs import *
# Plotting functions that are used in multiple types of calculations
from pyMELTScalc.Plotting import *
# This has the main functions that call the barometry calculations
from pyMELTScalc.Barom import *
# This contains the functions required to find liquidi via subprocesses
from pyMELTScalc.Liq import *
# This contains the functions required to perform crystallisation calculations via subprocesses
from pyMELTScalc.Crystallise import *
# This contains the functions used to call MELTS calculations
from pyMELTScalc.MELTS import *
# This contains the functions used to perform calculations using the Holland et al. 2018 thermodynamic dataset
try:
    from pyMELTScalc.Holland import *
except:
    print('Warning: MAGEMin calculations cannot be performed')

# version
from ._version import __version__
