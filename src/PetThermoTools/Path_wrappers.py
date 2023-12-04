import numpy as np
import pandas as pd
from PetThermoTools.Path import *

def isobaric_crystallisation(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, Crystallinity_limit = None, fluid_sat = None, timeout = None):
    '''
    Computes isobaric crystallization paths using the multi_path function,which
    computes the phase equilibria between solid, liquid and fluid phases.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    bulk: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each step. Default False.

    T_path_C: float or np.ndarray
        Used to specify the temperature path if non constant steps are used.
    
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

    Crystallinity_limit: float
        If not none, the calculation will stop once the weight fraction of crytals in the system exceeds the number specified.

    fluid_sat: True/False
        If True, the calculation will start at the fluid saturated liquidus, but no excess volatile phase will be present at the start of the calculation.

    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.
    '''
    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, Crystallinity_limit = Crystallinity_limit, fluid_sat = fluid_sat, timeout = timeout)

    return Results

def isochoric_crystallisation(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, Crystallinity_limit = None, fluid_sat = None, timeout = None):
    '''
    Computes isochorics crystallization paths using the multi_path function,which
    computes the phase equilibria between solid, liquid and fluid phases.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    bulk: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each step. Default False.

    T_path_C: float or np.ndarray
        Used to specify the temperature path if non constant steps are used.
    
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

    Crystallinity_limit: float
        If not none, the calculation will stop once the weight fraction of crytals in the system exceeds the number specified.

    fluid_sat: True/False
        If True, the calculation will start at the fluid saturated liquidus, but no excess volatile phase will be present at the start of the calculation.

    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.
    '''
    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, Crystallinity_limit = Crystallinity_limit, fluid_sat = fluid_sat, timeout = timeout, isochoric = True)

    return Results

def polybaric_crystallisation_path(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None):
    '''
    Computes polybaric crystallization paths using the multi_path function, which
    computes the phase equilibria between solid, liquid and fluid phases, where
    pressure and temperature both change continuously throughout the crystallization
    model.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    bulk: Dict or pd.DataFrame
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

    P_start_bar: float or np.ndarray
        Initial pressure used for path calculations.

    P_end_bar: float or np.ndarray
        Final pressure in crystallisation calculations. If T_end_C is not None this value will not be used.

    dp_bar: float or np.ndarray
        Pressure increment during crystallisation calculations.

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
    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout)

    return Results

def polybaric_crystallisation_onestep(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, T_step_C = None, P_start_bar = None, P_end_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None):

    comp = bulk.copy()

    if find_liquidus is not None and T_start_C is None:
        T_start_C = 1500

    if type(T_step_C) != np.ndarray and type(P_start_bar) != np.ndarray and type(P_end_bar) != np.ndarray:
        T_path_C = np.linspace(T_start_C, T_end_C, int(1+(T_start_C - T_end_C)/2))
        P_path_bar = np.zeros(len(T_path_C))

        P_path_bar[np.where(T_path_C > T_step_C)] = P_start_bar
        P_path_bar[np.where(T_path_C < T_step_C)] = P_end_bar
    else:
        if type(P_start_bar) == np.ndarray:
            L = len(P_start_bar)
        elif type(P_end_bar) == np.ndarray:
            L = len(P_end_bar)
        else:
            L = len(T_step_C)

        if type(P_start_bar) == float or type(P_start_bar) == int:
            P_start_bar = np.zeros(int(L)) + P_start_bar
        if type(P_end_bar) == float or type(P_end_bar) == int:
            P_end_bar = np.zeros(int(L)) + P_end_bar
        if type(T_step_C) == float or type(T_step_C) == int:
            T_step_C = np.zeros(int(L)) + T_step_C

        T_path_C = np.linspace(T_start_C, T_end_C, int(1+(T_start_C - T_end_C)/2))
        P_path_bar = np.zeros((len(T_step_C),len(T_path_C)))
        for i in range(len(T_step_C)):
            P_path_bar[i,np.where(T_path_C > T_step_C[i])] = P_start_bar[i]
            P_path_bar[i,np.where(T_path_C < T_step_C[i])] = P_end_bar[i]

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, P_path_bar = P_path_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout)

    return Results

def polybaric_crystallisation_multistep(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, T_step_C = None, P_step_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None):

    comp = bulk.copy()

    if find_liquidus is not None and T_start_C is None:
        T_start_C = 1500

    if type(T_step_C) == np.ndarray and type(P_step_bar) == np.ndarray:
        T_path_C = np.linspace(T_start_C, T_end_C, int(1+(T_start_C - T_end_C)/2))
        P_path_bar = np.zeros(len(T_path_C))

        T_sort = np.sort(T_step_C)
        T_step_C = T_sort[::-1]

        P_sort = np.sort(P_step_bar)
        P_step_bar = P_sort[::-1]

        for i in range(len(T_step_C)):
            if i == 0:
                P_path_bar[np.where(T_path_C > T_step_C[i])] = P_step_bar[i]
            else:
                P_path_bar[np.where((T_path_C <= T_step_C[i-1]) & (T_path_C > T_step_C[i]))] = P_step_bar[i]

        P_path_bar[np.where(T_path_C <= T_step_C[-1])] = P_step_bar[-1]
    elif type(T_step_C) == list or type(P_step_bar) == list:
        T_path = np.linspace(T_start_C, T_end_C, int(1+(T_start_C - T_end_C)/2))

        if type(T_step_C) == list:
            L = len(T_step_C)
        else:
            L = len(P_step_bar)

        T_path_C = np.zeros((int(L), len(T_path)))
        P_path_bar = np.zeros((int(L), len(T_path)))
        for i in range(int(L)):
            T_path_C[i, :] = T_path

        if type(T_step_C) == np.ndarray:
            T_sort = np.sort(T_step_C)
            T_step_C = T_sort[::-1]
            for i in range(len(P_step_bar)):
                P_sort = np.sort(P_step_bar[i])
                P_step = P_sort[::1]

                for j in range(len(T_step_C)):
                    if i == 0:
                        P_path_bar[i,np.where(T_path_C[i,:] > T_step_C[j])] = P_step[j]
                    else:
                        P_path_bar[i,np.where((T_path_C[i,:] <= T_step_C[j-1]) & (T_path_C[i,:] > T_step_C[j]))] = P_step[j]

                P_path_bar[i, np.where(T_path_C[i,:] <= T_step_C[-1])] = P_step[-1]

        elif type(P_step_bar) == np.ndarray:
            P_sort = np.sort(P_step_bar)
            P_step_bar = P_sort[::-1]
            for i in range(len(T_step_C)):
                T_sort = np.sort(T_step_C[i])
                T_step = T_sort[::1]

                for j in range(len(P_step_bar)-1):
                    if i == 0:
                        P_path_bar[i,np.where(T_path_C[i,:] > T_step[j])] = P_step_bar[j]
                    else:
                        P_path_bar[i,np.where((T_path_C[i,:] <= T_step[j-1]) & (T_path_C[i,:] > T_step[j]))] = P_step_bar[j]

                P_path_bar[i, np.where(T_path_C[i,:]<= T_step[-1])] = P_step_bar[-1]

        else:
            for i in range(len(T_step_C)):
                T_sort = np.sort(T_step_C[i])
                T_step = T_sort[::1]
                P_sort = np.sort(P_step_bar[i])
                P_step = P_sort[::1]

                for j in range(len(T_step_C[i])):
                    if i == 0:
                        P_path_bar[i,np.where(T_path_C[i,:] > T_step[j])] = P_step[j]
                    else:
                        P_path_bar[i,np.where((T_path_C[i,:] <= T_step[j-1]) & (T_path_C[i,:] > T_step[j]))] = P_step[j]

                P_path_bar[i, np.where(T_path_C[i,:]<= T_step[-1])] = P_step[-1]

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, P_path_bar = P_path_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout)

    return Results

def isothermal_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, fluid_sat = fluid_sat)

    return Results

def isentropic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, isentropic = True, fluid_sat = fluid_sat)

    return Results

def isenthalpic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, isenthalpic = True, fluid_sat = fluid_sat)

    return Results