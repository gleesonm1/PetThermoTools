import numpy as np
import pandas as pd
from PetThermoTools.Path import *

def isobaric_crystallisation(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, Crystallinity_limit = None, fluid_sat = None, timeout = None):
    '''
    Simulates isobaric crystallization paths by calling the `multi_path` function, which computes the phase equilibria between solid, liquid, and fluid phases under specified conditions.

    Parameters
    ----------
    Model: str, optional
        The thermodynamic model to use for calculations. Options include:
        - "Holland" for Holland et al. phase equilibrium models.
        - MELTS versions can be specified as "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default is "MELTSv1.0.2".
    bulk: dict or pandas.DataFrame
        The initial bulk composition for the calculations. If a dictionary is provided, it is used for all calculations. If a DataFrame is provided, it can handle multiple compositions, with rows corresponding to different calculations.
    Frac_solid: bool, optional
        If True, solid phases are fractionated (removed) after each step of crystallization. Default is False.
    Frac_fluid: bool, optional
        If True, fluid phases are fractionated (removed) after each step. Default is False.
    T_path_C: float or numpy.ndarray, optional
        Temperature path for non-linear temperature steps. If specified, overrides the use of `T_start_C`, `T_end_C`, and `dt_C`.
    T_start_C: float, optional
        The starting temperature for crystallization, in degrees Celsius.
    T_end_C: float, optional
        The final temperature for crystallization, in degrees Celsius.
    dt_C: float, optional
        The temperature step size for each iteration, in degrees Celsius.
    P_bar: float or numpy.ndarray, optional
        The pressure at which crystallization occurs, in bars. Default assumes constant pressure (isobaric conditions).
    Fe3Fet_Liq: float or numpy.ndarray, optional
        The initial Fe³⁺/ΣFe ratio in the melt phase. If an array is provided, each entry corresponds to a different bulk composition. If None, this value must be specified in `bulk` or by using an oxygen fugacity buffer.
    H2O_Liq: float or numpy.ndarray, optional
        Initial H₂O content in the melt phase. Similar to `Fe3Fet_Liq`, an array can specify different water contents for each bulk composition.
    CO2_Liq: float or numpy.ndarray, optional
        Initial CO₂ content in the melt phase. Similar in usage to `H2O_Liq`.
    find_liquidus: bool, optional
        If True, the starting temperature is adjusted to the liquidus temperature of the system before crystallization. Default is False.
    fO2_buffer: str, optional
        Specifies an oxygen fugacity buffer to control the system's redox state. Options include "FMQ" (Fayalite-Magnetite-Quartz) or "NNO" (Nickel-Nickel Oxide).
    fO2_offset: float or numpy.ndarray, optional
        Log units of offset from the chosen `fO2_buffer`. For instance, an offset of -1 means 1 log unit below the buffer.
    label: str, optional
        A label or identifier for the simulation run. Useful for tracking multiple results.
    Crystallinity_limit: float, optional
        If specified, the simulation will stop when the total crystallinity (solid fraction) reaches the given limit.
    fluid_sat: bool, optional
        If True, the calculation starts from fluid saturation but no excess volatile phase is present initially. Default is False.
    timeout: float, optional
        Maximum allowed time for the simulation in seconds. If exceeded, the calculation will terminate early.

    Returns
    -------
    Results : dict
        Dictionary containing results of the isothermal decompression simulation, with compositions 
        and thermodynamic properties for each step contained within phase-specific dataframes.
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
    """
    Perform isothermal decompression calculations using the provided model and input parameters.

    This function calls the `multi_path` function to simulate decompression at constant temperature.

    Parameters:
    ----------
    Model : string
        Thermodynamic model to use (e.g., "MELTSv...", "Holland"). MELTS options are "pMELTS", "MELTSv1.0.2", "MELTSv1.1.0", and "MELTSv1.2.0".
    bulk : dict or pd.DataFrame
        Initial bulk composition for calculations. If a dictionary, the same composition is used 
        across all steps.
    Frac_solid : bool
        Whether solid phases are fractionated at each step. Default is None.
    Frac_fluid : bool
        Whether fluid phases are fractionated at each step. Default is None.
    T_C : float
        The constant temperature (in Celsius) for the decompression.
    P_start_bar : float
        Initial pressure in bars for the decompression.
    P_end_bar : float
        Final pressure in bars for the decompression.
    dp_bar : float
        Pressure step size in bars.
    Fe3Fet_Liq : float
        Fe3+/total Fe ratio for the liquid phase.
    H2O_Liq : float
        Initial H2O content in the liquid phase.
    CO2_Liq : float
        Initial CO2 content in the liquid phase.
    find_liquidus : bool
        Whether to search for the liquidus temperature before decompression. Default is None.
    fO2_buffer : string
        Buffer for oxygen fugacity (e.g., "FMQ", "NNO"). Default is None.
    fO2_offset : float
        Offset from the oxygen fugacity buffer in log units.
    label : string
        Optional label for the calculation. Default is None.
    timeout : int
        Time limit for the calculation in seconds. Default is None.
    fluid_sat : bool
        Whether the system is fluid-saturated. Default is None.

    Returns:
    ----------
    Results : dict
        Dictionary containing results of the isothermal decompression simulation, with compositions 
        and thermodynamic properties for each step contained within phase-specific dataframes.
    """
    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, fluid_sat = fluid_sat)

    return Results

def isentropic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):
    """
    Perform isentropic decompression calculations using the provided model and input parameters.

    This function calls the `multi_path` function to simulate decompression at constant entropy.

    Parameters:
    ----------
    Model : string
        Thermodynamic model to use (e.g., "MELTSv...", "Holland"). MELTS options are "pMELTS", "MELTSv1.0.2", "MELTSv1.1.0", and "MELTSv1.2.0".
    bulk : dict or pd.DataFrame
        Initial bulk composition for calculations. If a dictionary, the same composition is used 
        across all steps.
    Frac_solid : bool
        Whether solid phases are fractionated at each step. Default is None.
    Frac_fluid : bool
        Whether fluid phases are fractionated at each step. Default is None.
    T_C : float
        The starting temperature (in Celsius) for the decompression.
    P_start_bar : float
        Initial pressure in bars for the decompression.
    P_end_bar : float
        Final pressure in bars for the decompression.
    dp_bar : float
        Pressure step size in bars.
    Fe3Fet_Liq : float
        Fe3+/total Fe ratio for the liquid phase.
    H2O_Liq : float
        Initial H2O content in the liquid phase.
    CO2_Liq : float
        Initial CO2 content in the liquid phase.
    find_liquidus : bool
        Whether to search for the liquidus temperature before decompression. Default is None.
    fO2_buffer : string
        Buffer for oxygen fugacity (e.g., "FMQ", "NNO"). Default is None.
    fO2_offset : float
        Offset from the oxygen fugacity buffer in log units.
    label : string
        Optional label for the calculation. Default is None.
    timeout : int
        Time limit for the calculation in seconds. Default is None.
    fluid_sat : bool
        Whether the system is fluid-saturated. Default is None.

    Returns:
    ----------
    Results : dict
        Dictionary containing results of the isentropic decompression simulation, with compositions 
        and thermodynamic properties for each step contained within phase-specific dataframes.
    """
    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, isentropic = True, fluid_sat = fluid_sat)

    return Results

def isenthalpic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, isenthalpic = True, fluid_sat = fluid_sat)

    return Results