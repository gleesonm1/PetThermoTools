__author__ = 'Matthew Gleeson'

import numpy as np
import pandas as pd
import sys
import os
import warnings
import site
import sysconfig
import importlib
import multiprocessing
from pathlib import Path
import matplotlib.pyplot as plt
# import subprocess
# import glob
from subprocess import Popen, PIPE
from os import path

# import importlib.util

# # Check if 'meltsdynamic' exists in the current Python path
# spec = importlib.util.find_spec("meltsdynamic")

# Child processes spawned by multiprocessing (every multi_* / *_multi worker)
# re-import the whole package from scratch, so these import-time messages
# would otherwise print once per worker (visible in a notebook/IDE kernel as
# the message repeated 100+ times for a batch run). Only the main process
# prints them. Note: `multiprocessing.parent_process()` is NOT a reliable
# signal here — it is only populated by `_bootstrap()`, which runs *after*
# the spawned child has already re-imported this module, so it still reads
# as None during this very import. `current_process().name` is set to
# something other than "MainProcess" (e.g. "Process-1") before the child's
# re-import begins, so it is the check that actually works at import time.
_IS_MAIN_PROCESS = multiprocessing.current_process().name == "MainProcess"

try:
    if os.path.exists('/home/jovyan/shared/Models/alphaMELTS'):
        sys.path.append('/home/jovyan/shared/Models/alphaMELTS')
    from meltsdynamic import MELTSdynamic
    if _IS_MAIN_PROCESS:
        print("alphaMELTS for Python files successfully located.")
except (ImportError, ModuleNotFoundError) as e:
    warnings.warn(f"Failed to import MELTSdynamic: {e}. \nUse `import sys; sys.path.append(r'path_to_alphaMELTS_4_python')` to add the alphaMELTS for Python files to the Python path.")

if _IS_MAIN_PROCESS:
    print("If using the Green et al. (2025) or Weller et al. (2024) thermodynamic models please run `ptt.activate_petthermotools_env()` prior to any calculations.")

# test installation script
from petthermotools.Installation import *
# General functions that are used to assess core numbers
from petthermotools.core_config import *
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
# Compsitions saved for use in petthermotools
from petthermotools.Methods import *
# # import trace element engine
# from petthermotools.trace_engine.partition_standards import *

# ── nGibbs neural-network emulator registry ──────────────────────────────────
# Instantiated once in nGibbs_bridge.py and imported here so users can
# inspect or extend the registry via `ptt._nGibbs_models`.
# Path.py also imports from nGibbs_bridge, so there is no duplication.

from petthermotools.nGibbs_bridge import _nGibbs_models

if not len(_nGibbs_models) and _IS_MAIN_PROCESS:
    print(
        "nGibbs not installed — nGibbs emulators are unavailable. "
        "Install it if you'd like with `ptt.install_nGibbs()` (~300 MB storage required)."
    )

# version
from ._version import __version__


