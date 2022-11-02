import numpy as np
import pandas as pd
from pyMELTScalc.Path import *

def Isobaric_crystallisation(Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None):
    '''
    Simple code to perform isobaric fractional crystallsiation calculations.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each  step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each  step. Default False.

    T_start_C: float or np.ndarray
        Initial temperature used for path calculations.

    T_end_C: float or np.ndarray
        Final temperature in crystallisation calculations.

    dt_C: float or np.ndarray
        Temperature increment during crystallisation calculations.

    P_bar: float or np.ndarray
        Calculation pressure - typically used when calculations are performed at a fixed P (e.g.,isobaric crystallisation).

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, and type(Fe3Fet_Liq) == np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their Fe3Fet_Liq value. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable or an oxygen fugacity buffer must be chosen.

    H2O_Liq: float or np.ndarray
        H2O content of the initial melt phase. If type(comp) == dict, and type(H2O_Liq) = np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their H2O_Liq value. If comp is a pd.DataFrame, a single H2O_Liq value may be passes (float) and will be used as the initial melt H2O content for all starting compositions. Alternatively, if an array of H2O_Liq values is passed, equal to the number of compositions specified in comp, a different initial melt H2O value will be passed for each sample. If None, H2O_Liq must be specified in the comp variable.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float or np.ndarray
        Offset from the buffer spcified in fO2_buffer (log units).

    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.
    '''

    Results = multi_path(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)

    return Results