import numpy as np
import pandas as pd
from petthermotools.Path import *

def isobaric_crystallisation(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                             T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
                             P_bar = None, Fe3Fet_init = None, Fe3Fet_Liq = None, 
                             H2O_init = None, CO2_init = None, H2O_Liq = None, CO2_Liq = None, 
                             find_liquidus = None, fO2_buffer = None, fO2_offset = None, 
                             label = None, Crystallinity_limit = None, fluid_sat = False, timeout = None,
                             Suppress = ['rutile', 'tridymite'], Suppress_except=False,
                             multi_processing = True):
    """
    Simulates isobaric (constant-pressure) crystallization paths using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
            MAGEMin:
                "Green2025"
                "Weller2024"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_path_C : float or ndarray, optional
        Custom temperature path. Overrides T_start_C, T_end_C, dt_C.
    T_start_C, T_end_C, dt_C : float, optional
        Linear temperature path definition.
    P_bar : float or ndarray, optional
        Pressure in bars for isobaric conditions.
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    Crystallinity_limit : float, optional
        Stops if total crystallinity exceeds this value.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                         T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, 
                         P_bar = P_bar, Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init = CO2_init, 
                         find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                         label = label, Crystallinity_limit = Crystallinity_limit, fluid_sat = fluid_sat, 
                         timeout = timeout, Suppress = Suppress, Suppress_except=Suppress_except,
                         multi_processing = multi_processing)

    return Results

def isochoric_crystallisation(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                              T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
                              P_bar = None, Fe3Fet_init = None, Fe3Fet_Liq = None, 
                              H2O_init = None, H2O_Liq = None, CO2_init = None, CO2_Liq = None, 
                              find_liquidus = None, fO2_buffer = None, fO2_offset = None, 
                              label = None, Crystallinity_limit = None, fluid_sat = False, 
                              timeout = None, Suppress = ['rutile', 'tridymite'], Suppress_except=False,
                              multi_processing = True):
    """
    Simulates isochoric (constant-volume) crystallization paths using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_path_C : float or ndarray, optional
        Custom temperature path. Overrides T_start_C, T_end_C, dt_C.
    T_start_C, T_end_C, dt_C : float, optional
        Linear temperature path definition.
    P_bar : float or ndarray, optional
        Pressure in bars for initial condition.
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    Crystallinity_limit : float, optional
        Stops if total crystallinity exceeds this value.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                         T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, 
                         P_bar = P_bar, Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init=CO2_init, 
                         find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                         label = label, Crystallinity_limit = Crystallinity_limit, fluid_sat = fluid_sat, 
                         timeout = timeout, isochoric = True, Suppress = Suppress, Suppress_except=Suppress_except,
                         multi_processing=multi_processing)

    return Results

def polybaric_crystallisation_path(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                                   T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
                                   P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
                                   Fe3Fet_init = None, Fe3Fet_Liq = None, H2O_init = None, H2O_Liq = None, 
                                   CO2_init = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, 
                                   fO2_offset = None, label = None, timeout = None, multi_processing=True, fluid_sat = False,
                                   Suppress = ['rutile', 'tridymite'], Suppress_except=False):
    """
    Simulates polybaric crystallization along a user-specified P-T path using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
            MAGEMin:
                "Green2025"
                "Weller2024"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_path_C : ndarray, optional
        Custom temperature path. Overrides T_start_C, T_end_C, dt_C.
    T_start_C, T_end_C, dt_C : float, optional
        Linear temperature path definition.
    P_path_bar : ndarray, optional
        Pressure path in bars for polybaric conditions.
    P_start_bar, P_end_bar, dp_bar : float, optional
        Linear pressure path definition.
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    Crystallinity_limit : float, optional
        Stops if total crystallinity exceeds this value.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                         T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, 
                         P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, 
                         Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init = CO2_init, find_liquidus = find_liquidus, 
                         fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                         label = label, timeout = timeout, Suppress = Suppress, Suppress_except=Suppress_except,
                         multi_processing=multi_processing, fluid_sat=fluid_sat)

    return Results

def polybaric_crystallisation_onestep(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                                      T_start_C = None, T_end_C = None, dt_C = None, T_step_C = None, 
                                      P_start_bar = None, P_end_bar = None, Fe3Fet_Liq = None, Fe3Fet_init = None,
                                      H2O_Liq = None, H2O_init = None, CO2_Liq = None, CO2_init = None, 
                                      find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, 
                                      timeout = None, multi_processing = True, Suppress = ['rutile', 'tridymite'], 
                                      Suppress_except=False, fluid_sat = False):
    """
    Simulates polybaric crystallization along a user-specified P-T path using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
            MAGEMin:
                "Green2025"
                "Weller2024"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_start_C, T_end_C, dt_C : float
        Linear temperature path definition.
    T_step_C: float
        Temperature at which pressure changes.
    P_start_bar, P_end_bar: float
        Initial and final pressure (bars).
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    Crystallinity_limit : float, optional
        Stops if total crystallinity exceeds this value.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

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

    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                        T_path_C = T_path_C, P_path_bar = P_path_bar, Fe3Fet_init = Fe3Fet_init, 
                        H2O_init = H2O_init, CO2_init = CO2_init, 
                        find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                        label = label, timeout = timeout, Suppress = Suppress, Suppress_except=Suppress_except,
                        multi_processing=multi_processing, fluid_sat=fluid_sat)

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

def isothermal_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                             T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
                             Fe3Fet_Liq = None, Fe3Fet_init = None, H2O_Liq = None, H2O_init = None,
                             CO2_Liq = None, CO2_init = None, find_liquidus = None, fO2_buffer = None, 
                             fO2_offset = None, label = None, timeout = None, fluid_sat = False, 
                             Suppress = ['rutile', 'tridymite'], Suppress_except=False, multi_processing = True):
    """
    Simulates isothermal (constant-temperature) decompression calculations using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
            MAGEMin:
                "Green2025"
                "Weller2024"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_C : float or ndarray
        Fixed temperature for each decompression calculation.
    P_start_bar, P_end_bar, dp_bar : float or ndarray
        Linear pressure path definition.
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                         T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, 
                         Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init = CO2_init, find_liquidus = find_liquidus, 
                         fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, 
                         fluid_sat = fluid_sat, Suppress = Suppress, Suppress_except=Suppress_except,
                         multi_processing=multi_processing)

    return Results

def isentropic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, 
                             T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
                             Fe3Fet_Liq = None, Fe3Fet_init = None, H2O_Liq = None, H2O_init = None,
                             CO2_Liq = None, CO2_init = None, find_liquidus = None, fO2_buffer = None, 
                             fO2_offset = None, label = None, timeout = None, fluid_sat = False, 
                             Suppress = ['rutile', 'tridymite'], Suppress_except=False, multi_processing = True):
    """
    Simulates isentropic (constant-entropy) decompression calculations using the `multi_path` function.

    Parameters
    ----------
    Model : str, optional
        Thermodynamic model, options:
            MELTS:
                "pMELTS"
                "MELTSv1.0.2"
                "MELTSv1.1.0"
                "MELTSv1.2.0"
    bulk : dict or pandas.DataFrame
        Initial bulk composition(s). Dictionary for a single run; DataFrame for multiple runs.
    Frac_solid : bool, optional
        If True, solid phases are fractionated. Default is False.
    Frac_fluid : bool, optional
        If True, fluid phases are fractionated. Default is False.
    T_C : float or ndarray
        Initial temperature for each decompression calculation.
    P_start_bar, P_end_bar, dp_bar : float or ndarray
        Linear pressure path definition.
    Fe3Fet_init, H2O_init, CO2_init : float or ndarray, optional
        Initial redox and volatile content of the system.
    find_liquidus : bool, optional
        If True, starts at the liquidus. Default is False.
    fO2_buffer : str, optional
        Redox buffer (e.g., "FMQ", "NNO").
    fO2_offset : float or ndarray, optional
        Log unit offset from fO2 buffer.
    label : str, optional
        Label for the run.
    fluid_sat : bool, optional
        If True, begins at fluid saturation. Default is False.
    timeout : float, optional
        Maximum runtime in seconds per calculation.
    Suppress : list of str, optional
        Phases to suppress. Default ['rutile', 'tridymite'].
    Suppress_except : bool, optional
        If True, suppresses all phases except those listed.
    multi_processing : bool, optional
        Enables multiprocessing. Default is True.

    Returns
    -------
    dict
        A dictionary of results with DataFrames for each phase.
    """
    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()
    Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                        T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, 
                        Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init = CO2_init, find_liquidus = find_liquidus, 
                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, 
                        isentropic = True, fluid_sat = fluid_sat, Suppress = Suppress, Suppress_except=Suppress_except,
                         multi_processing=multi_processing)

    return Results

# def isenthalpic_decompression(Model = None, bulk = None, Frac_solid = None, Frac_fluid = None, T_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, label = None, timeout = None, fluid_sat = None):

#     comp = bulk.copy()
#     Results = multi_path(Model = Model, bulk = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, label = label, timeout = timeout, isenthalpic = True, fluid_sat = fluid_sat)

#     return Results