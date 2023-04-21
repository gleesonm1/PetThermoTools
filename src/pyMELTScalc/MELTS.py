import numpy as np
import pandas as pd
import sys
import time

def equilibrate_MELTS(Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    Results = {}

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp.copy()
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    from meltsdynamic import MELTSdynamic
    if Model is None or Model == "MELTSv1.0.2":
        melts = MELTSdynamic(1)
    elif Model == "pMELTS":
        melts = MELTSdynamic(2)
    elif Model == "MELTSv1.1.0":
        melts = MELTSdynamic(3)
    elif Model == "MELTSv1.2.0":
        melts = MELTSdynamic(4)

    melts.engine.setBulkComposition(bulk)
    melts.engine.pressure = P_bar
    melts.engine.temperature = T_C

    length = 1

    if fO2_buffer is not None:
        if fO2_offset is None:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
        else:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    # PhaseList = [None]
    # PhaseComp = {}
    # PhaseProp = {}
    # Props = ['mass', 'rho']

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 7)), columns = ['temperature', 'pressure', 'g', 'h', 's', 'v', 'dvdp'])
    # Results['liquid1'] = pd.DataFrame(data = np.zeros((length+1, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    # Results['liquid1_prop'] = pd.DataFrame(data = np.zeros((length+1, 4)), columns = ['h', 'mass', 'v', 'rho'])

    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        return Results

    SolidPhase = melts.engine.solidNames
    # if np.isnan(melts.engine.getProperty('mass', 'liquid1')):
    #     PhaseList = SolidPhase
    # elif len(SolidPhase) > 0:
    #     PhaseList = ['liquid1'] + SolidPhase
    # else:
    #     PhaseList = ['Liquid1']
    if SolidPhase is not None:
        PhaseList = ['liquid1'] + SolidPhase
    else:
        PhaseList = ['liquid1']

    print(PhaseList)

    for R in Results['Conditions']:
        if R == 'temperature':
            Results['Conditions'][R].loc[0] = melts.engine.temperature
        elif R == 'pressure':
            Results['Conditions'][R].loc[0] = melts.engine.pressure
        else:
            Results['Conditions'][R].loc[0] = melts.engine.getProperty(R, 'bulk')

    for phase in PhaseList:
        if phase not in list(Results.keys()):
            Results[phase] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
            Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((length, 5)), columns = ['g','h', 'mass', 'v', 'rho'])

        for el in Results[phase]:
            Results[phase][el].loc[0] = melts.engine.getProperty('dispComposition', phase, el)

        for pr in Results[phase + '_prop']:
            Results[phase + '_prop'][pr].loc[0] = melts.engine.getProperty(pr, phase)

    # if PhaseList is not None:
    #     for p in PhaseList:
    #         PhaseComp[p] = melts.engine.getProperty('dispComposition', p)
    #         PhaseProp[p] = {}
    #         for i in Props:
    #             PhaseProp[p][i] = melts.engine.getProperty(i, p)

    return Results

def findCO2_MELTS(P_bar = None, Model = None, T_C = None, comp = None, melts = None, fO2_buffer = None, fO2_offset = None, Step = None):

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp.copy()
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    from meltsdynamic import MELTSdynamic

    if melts is None:
        if Model is None or Model == "MELTSv1.0.2":
            melts = MELTSdynamic(1)
        elif Model == "pMELTS":
            melts = MELTSdynamic(2)
        elif Model == "MELTSv1.1.0":
            melts = MELTSdynamic(3)
        elif Model == "MELTSv1.2.0":
            melts = MELTSdynamic(4)

    T_Liq = 0
    H2O = 0
    CO2 = 0

    try:
        melts.engine.setBulkComposition(bulk)
        melts.engine.pressure = P_bar
        melts.engine.temperature = T_C
    except:
        return T_Liq, H2O, CO2

    if fO2_buffer is not None:
        if fO2_offset is None:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
        else:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    Liq = ['liquid1','water1', 'fluid1']
    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        return T_Liq, H2O, CO2

    PhaseList = melts.engine.solidNames
    if PhaseList is None:
        PhaseList = ['liquid1']
    else:
        PhaseList = ['liquid1'] + PhaseList

    i = set.intersection(set(Liq),set(PhaseList))

    if Step is None:
        Step = np.array([3,1,0.1])


    for j in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                try:
                    melts = melts.addNodeAfter()
                    melts.engine.temperature = melts.engine.temperature - Step[j]
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    return T_Liq, H2O, CO2

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                try:
                    melts = melts.addNodeAfter()
                    melts.engine.temperature = melts.engine.temperature + Step[j]
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    return T_Liq, H2O, CO2

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))


    if "fluid1" not in PhaseList:
        CO2_step = np.array([0.1, 0.02, 0.005])
        for j in range(len(CO2_step)):
            while "fluid1" not in PhaseList:
                bulk[15] = bulk[15] + CO2_step[j]
                melts.engine.setBulkComposition(bulk)
                if len(i) == len(PhaseList):
                    while len(i) == len(PhaseList):
                        try:
                            melts = melts.addNodeAfter()
                            melts.engine.temperature = melts.engine.temperature - 0.1
                            melts.engine.calcEquilibriumState(1,0)
                        except:
                            return T_Liq, H2O, CO2

                        PhaseList = melts.engine.solidNames
                        if PhaseList is None:
                            PhaseList = ['liquid1']
                        else:
                            PhaseList = ['liquid1'] + PhaseList
                        i = set.intersection(set(Liq),set(PhaseList))

                if len(i) < len(PhaseList):
                    while len(i) < len(PhaseList):
                        try:
                            melts = melts.addNodeAfter()
                            melts.engine.temperature = melts.engine.temperature + 0.1
                            melts.engine.calcEquilibriumState(1,0)
                        except:
                            return T_Liq, H2O, CO2

                        PhaseList = melts.engine.solidNames
                        if PhaseList is None:
                            PhaseList = ['liquid1']
                        else:
                            PhaseList = ['liquid1'] + PhaseList
                        i = set.intersection(set(Liq),set(PhaseList))

            if j != len(CO2_step) - 1:
                bulk[15] = bulk[15] - CO2_step[j]

    T_Liq = melts.engine.getProperty('temperature')
    H2O = melts.engine.getProperty('dispComposition', 'liquid1', 'H2O')
    CO2 = melts.engine.getProperty('dispComposition', 'liquid1', 'CO2')

    return T_Liq, H2O, CO2

def findLiq_MELTS(P_bar = None, Model = None, T_C_init = None, comp = None, melts = None, fO2_buffer = None, fO2_offset = None, Step = None, CO2_return = None, bulk_return = None):
    '''
    Perform a single find liquidus calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.

    Parameters:
    ----------
    P_bar: float
        Specifies the pressure of calculation (bar).

    Model: string
        Dictates the MELTS model to be used: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    T_C_init: float
        Initial 'guess' temperature for findLiq calculations (degrees C).

    comp: list or dict
        Input oxide values required for the calculations.

    melts: melts.engine
        If used as part of a crystallisation or decompression model, an instance of the MELTS C library will have been previously initiated and can be loaded here. When liquidus calculations are run in isolation this should be kept blank and a new instance of the MELTS C library will be initiated.

    Returns:
    ----------
    T_Liq_C: np.ndarray
        Array of liquidus temperatures.

    H2O: np.ndarray
        Array of melt H2O contents at the liquidus.

    '''

    if P_bar is None:
        raise Exception("Please specify a pressure for calculations")

    if Model is None:
        Model = "MELTSv1.0.2"

    if T_C_init is None:
        T_C_init = 1300

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    T_Liq = 0
    H2O_Melt = 0
    CO2_Melt = 0

    from meltsdynamic import MELTSdynamic

    if melts is None:
        if Model is None or Model == "MELTSv1.0.2":
            melts = MELTSdynamic(1)
        elif Model == "pMELTS":
            melts = MELTSdynamic(2)
        elif Model == "MELTSv1.1.0":
            melts = MELTSdynamic(3)
        elif Model == "MELTSv1.2.0":
            melts = MELTSdynamic(4)

    try:
        melts.engine.setBulkComposition(bulk)
        melts.engine.pressure = P_bar
        melts.engine.temperature = T_C_init
    except:
        return T_Liq, H2O_Melt

    if fO2_buffer is not None:
        if fO2_offset is None:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
        else:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    Liq = ['liquid1','water1', 'fluid1']
    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        if CO2_return is not None:
                return T_Liq, H2O_Melt, CO2_Melt

        elif bulk_return is not None:
            return T_Liq, H2O_Melt, melts

        else:
            return T_Liq, H2O_Melt

    PhaseList = melts.engine.solidNames
    if PhaseList is None:
        PhaseList = ['liquid1']
    else:
        PhaseList = ['liquid1'] + PhaseList

    i = set.intersection(set(Liq),set(PhaseList))

    if Step is None:
        Step = np.array([3,1,0.1])

    for j in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                try:
                    melts = melts.addNodeAfter()
                    melts.engine.temperature = melts.engine.temperature - Step[j]
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    if CO2_return is not None:
                            return T_Liq, H2O_Melt, CO2_Melt

                    elif bulk_return is not None:
                        return T_Liq, H2O_Melt, melts

                    else:
                        return T_Liq, H2O_Melt

                    break

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                try:
                    melts = melts.addNodeAfter()
                    melts.engine.temperature = melts.engine.temperature + Step[j]
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    if CO2_return is not None:
                            return T_Liq, H2O_Melt, CO2_Melt

                    elif bulk_return is not None:
                        return T_Liq, H2O_Melt, melts

                    else:
                        return T_Liq, H2O_Melt

                    break

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))


    T_Liq = melts.engine.temperature
    H2O_Melt = melts.engine.getProperty('dispComposition', 'liquid1', 'H2O')

    if CO2_return is not None:
        CO2_Melt = melts.engine.getProperty('dispComposition', 'liquid1', 'CO2')
        return T_Liq, H2O_Melt, CO2_Melt

    elif bulk_return is not None:
        return T_Liq, H2O_Melt, melts

    else:
        return T_Liq, H2O_Melt


def phaseSat_MELTS(Model = None, comp = None, phases = None, T_initial_C = None, T_step_C = None, dt_C = None, P_bar = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
    '''
    Perform a single crystallisation calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.

    Parameters:
    ----------
    Model: string
        Dictates the MELTS model to be used: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: dict
        Input oxide values required for the calculations.

    phases: list
        phases of interest

    T_initial_C: float
        Initial temperature used for liquidus calculations.

    T_step_C: float
        Temperature step at each point of the model.

    dt_C: float
        Total temperature change allowed in the model.

    P_bar: float
        Pressure of the calculation.

    Returns:
    ----------
    Results: Dict
        Dict containing a float for each saturation temperature found and the T_Liq and melt H2O values.

    '''
    Results = {'a_sat': np.nan, 'b_sat': np.nan, 'c_sat': np.nan, 'T_Liq': np.nan, 'H2O_melt': np.nan}
    if len(phases) == 2:
        del Results['c_sat']

    from meltsdynamic import MELTSdynamic

    if Model is None or Model == "MELTSv1.0.2":
        melts = MELTSdynamic(1)
    elif Model == "pMELTS":
        melts = MELTSdynamic(2)
    elif Model == "MELTSv1.1.0":
        melts = MELTSdynamic(3)
    elif Model == "MELTSv1.2.0":
        melts = MELTSdynamic(4)

    melts.engine.setSystemProperties("Suppress", "tridymite")

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    try:
        Results['T_Liq'], Results['H2O_melt'] = findLiq_MELTS(P_bar = P_bar, comp = bulk, T_C_init = T_initial_C, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
    except:
        return Results

    if type(H2O_Liq) == np.ndarray:
        if Results['H2O_melt'] < 0.99*bulk[14]:
            return Results


    T = Results['T_Liq']
    T_final = T - dt_C

    if T>0:
        while T >= T_final:
            try:
                melts = melts.addNodeAfter()
                melts.engine.setBulkComposition(bulk)
                melts.engine.pressure = P_bar
                melts.engine.temperature = T
            except:
                return Results

            try:
                melts.engine.calcEquilibriumState(1,0)
            except:
                return Results

            try:
                PhaseList = melts.engine.solidNames
                print(PhaseList)
                if 'tridymite1' in PhaseList:
                    PhaseList = ['quartz1'] + PhaseList
                if 'clinopyroxene2' in PhaseList:
                    PhaseList = ['orthopyroxene1'] + PhaseList

                if phases[0] in PhaseList and np.isnan(Results['a_sat']):# == 0:
                    Results['a_sat'] = melts.engine.temperature

                if phases[1] in PhaseList and np.isnan(Results['b_sat']):# == 0:
                    Results['b_sat'] = melts.engine.temperature

                if len(phases) == 3:
                    if phases[2] in PhaseList and np.isnan(Results['c_sat']):# == 0:
                        Results['c_sat'] = melts.engine.temperature

                    if ~np.isnan(Results['a_sat']) and ~np.isnan(Results['b_sat']) and ~np.isnan(Results['c_sat']):# > 0:
                        break

                if len(phases) == 2:
                    if ~np.isnan(Results['a_sat']) and ~np.isnan(Results['b_sat']):# > 0:
                        break

                T = T - T_step_C
            except:
                T = T - T_step_C


    return Results

def path_MELTS(Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, fluid_sat = None, Crystallinity_limit = None):
    '''
    Perform a single  calculation in MELTS. WARNING! Running this function directly from the command land/jupyter notebook will initiate the MELTS C library in the main python process. Once this has been initiated the MELTS C library cannot be re-loaded and failures during the calculation will likely cause a terminal error to occur.

    Parameters:
    ----------
    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict
        Initial compositon for calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each crystallisation step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each crystallisation step. Default False.

    T_C: float
        Calculation temperature - typically used when calculations are performed at a fixed T (e.g.,isothermal degassing).

    T_path_C: np.ndarray
        If a specified temperature path is to be used, T_path_C will be used to determine the T at each step of the model. If 2D, this indicates that multiple calculations with different T_path_C arrays are to be performed.

    T_start_C: float
        Initial temperature used for path calculations.

    T_end_C: float
        Final temperature in crystallisation calculations.
melt
    dt_C: float
        Temperature increment during crystallisation calculations.

    P_bar: float
        Calculation pressure - typically used when calculations are performed at a fixed P (e.g.,isobaric crystallisation).

    P_path_bar: np.ndarray
        If a specified pressure path is to be used, P_path_bar will be used to determine the P at each step of the model.

    P_start_bar: float
        Initial pressure used for path calculations.

    P_end_bar: float
        Final pressure in crystallisation calculations.

    dp_bar: float
        Pressure increment during crystallisation calculations.

    isenthalpic: True/False
        If True, calculations will be performed at a constant enthalpy with T treated as a dependent variable.

    isentropic: True/False
        If True, calculations will be performed at a constant entropy with T treated as a dependent variable.

    isochoric: True/False
        If True, the volume of the system will be held constant instead of the pressure. Default is False.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float
        Offset from the buffer spcified in fO2_buffer (log units).

    Crystallinity_limit: float
        Crystallinity of the system (as a fraction) to determine when the calculation will be terminated.    

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    '''
    Results = {}

    if comp is None:
        raise Exception("No composition specified")

    if P_bar is not None and P_path_bar is None:
        P_path_bar = P_bar
    if T_C is not None and T_start_C is None:
        T_start_C = T_C

    if P_path_bar is None and P_start_bar is None:
        raise Exception("Initial P system must be defined")
    if T_path_C is None and T_start_C is None and find_liquidus is None:
        raise Exception("Starting temperature must be specified or the liquidus must be found")

    from meltsdynamic import MELTSdynamic

    if Model is None or Model == "MELTSv1.0.2":
        melts = MELTSdynamic(1)
    elif Model == "pMELTS":
        melts = MELTSdynamic(2)
    elif Model == "MELTSv1.1.0":
        melts = MELTSdynamic(3)
    elif Model == "MELTSv1.2.0":
        melts = MELTSdynamic(4)

    melts.engine.setSystemProperties("Suppress", "rutile")

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    if find_liquidus is not None:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_path_bar[0], comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
                else:
                    T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_path_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
            except:
                return Results
        elif P_start_bar is not None:
            try:
                T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_start_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
            except:
                return Results

        T_start_C = T_Liq

    else:
        if fO2_buffer is not None:
            if fO2_offset is None:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
            else:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    if T_path_C is None:
        if T_end_C is None and dt_C is None:
            T = T_start_C
        elif T_end_C is not None and dt_C is not None:
            T = np.linspace(T_start_C, T_end_C, 1+round((T_start_C-T_end_C)/dt_C))
    elif T_path_C is not None:
        T = T_path_C

    if P_path_bar is None:
        if P_end_bar is None and dp_bar is None:
            P = P_start_bar
        elif P_end_bar is not None and dp_bar is not None:
            P = np.linspace(P_start_bar, P_end_bar, 1+round((P_start_bar-P_end_bar)/dp_bar))
    elif P_path_bar is not None:
        P = P_path_bar

    if type(T) == np.ndarray and P_end_bar is None and dp_bar is not None:
        P = np.linspace(P_start_bar, P_start_bar - dp_bar*(len(T)-1), len(T))
    elif type(P) == np.ndarray and T_end_C is None and dt_C is not None:
        T = np.linspace(T_start_C, T_start_C - dt_C*(len(P)-1), len(P))

    if type(T) == np.ndarray and type(P) == np.ndarray:
        if len(T) != len(P):
            raise Exception("Length of P and T vectors are not the same. Check input parameters")

    if find_liquidus is None:
        if type(T) != np.ndarray:
            melts.engine.temperature = T
        else:
            melts.engine.temperature = T[0]

        if type(P) != np.ndarray:
            melts.engine.pressure = P
        else:
            melts.engine.pressure = P[0]

    if find_liquidus is not None:
        melts.engine.temperature = T_Liq

        if P_path_bar is not None or T_path_C is not None:
            if type(P) == np.ndarray and type(T) == np.ndarray:
                T_Liq_loc = np.abs(T - T_Liq).argmin()
                if T[T_Liq_loc]>T_Liq:
                    T = T[T_Liq_loc:]
                    P = P[T_Liq_loc:]
                else:
                    T = T[T_Liq_loc-1:]
                    P = P[T_Liq_loc-1:]

    melts = melts.addNodeAfter()
    melts.engine.setBulkComposition(bulk)

    if fluid_sat is not None:
        melts.engine.setSystemProperties("Mode", "Fractionate Fluids")
        melts.engine.calcEquilibriumState(1,1)

    if type(T) == np.ndarray:
        length = len(T)
    else:
        length = len(P)

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 7)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'dvdp', 'logfO2'])
    Results['liquid1'] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    Results['liquid1_prop'] = pd.DataFrame(data = np.zeros((length, 4)), columns = ['h', 'mass', 'v', 'rho'])

    for i in range(length):
        if type(T) == np.ndarray:
            melts.engine.temperature = T[i]
        if type(P) == np.ndarray:
            melts.engine.pressure = P[i]

        if Frac_solid is not None:
            melts.engine.setSystemProperties("Mode", "Fractionate Solids")

        if Frac_fluid is not None:
            melts.engine.setSystemProperties("Mode", "Fractionate Fluids")

        if isenthalpic is not None:
            if i == 0:
                melts.engine.setSystemProperties("Mode", "Isenthalpic")
                melts.engine.calcEquilibriumState(1,0)
                melts.engine.setSystemProperties("Mode", "Isenthalpic")

        if isentropic is not None:
            if i == 0:
                melts.engine.setSystemProperties("Mode", "Isentropic")
                melts.engine.calcEquilibriumState(1,0)
                melts.engine.setSystemProperties("Mode", "Isentropic")

        if isochoric is not None:
            if i == 0:
                melts.engine.setSystemProperties("Mode", "Isochoric")
                melts.engine.calcEquilibriumState(1,0)
                melts.engine.setSystemProperties("Mode", "Isochoric")

        if isochoric is None and isenthalpic is None and isentropic is None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(1,1)
                else:
                    melts.engine.calcEquilibriumState(1,0)
            except:
                # return Results
                break

        if isenthalpic is not None:
            try:
                melts.engine.calcEquilibriumState(2,0)
            except:
                # return Results
                break

        if isentropic is not None:
            try:
                melts.engine.calcEquilibriumState(3,0)
            except:
                # return Results
                break

        if isochoric is not None:
            try:
                melts.engine.calcEquilibriumState(4,0)
            except:
                # return Results
                break

        for R in Results['Conditions']:
            if R == 'temperature':
                Results['Conditions'][R].loc[i] = melts.engine.temperature
            elif R == 'pressure':
                Results['Conditions'][R].loc[i] = melts.engine.pressure
            elif R == 'logfO2':
                try:
                    Results['Conditions'][R].loc[i] = melts.engine.getProperty(R)
                except:
                    Results['Conditions'][R].loc[i] = np.nan
            else:
                Results['Conditions'][R].loc[i] = melts.engine.getProperty(R, 'bulk')

        try:
            PhaseList = ['liquid1'] + melts.engine.solidNames
        except:
            # return Results
            break

        for phase in PhaseList:
            if phase not in list(Results.keys()):
                Results[phase] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
                Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((length, 4)), columns = ['h', 'mass', 'v', 'rho'])

            if phase in list(Results.keys()):
                for el in Results[phase]:
                    Results[phase][el].loc[i] = melts.engine.getProperty('dispComposition', phase, el)

                for pr in Results[phase + '_prop']:
                    Results[phase + '_prop'][pr].loc[i] = melts.engine.getProperty(pr, phase)

        if Crystallinity_limit is not None:
            Volume = 0
            Fluid_List = ['liquid1', 'water1', 'fluid1']
            for phase in PhaseList:
                if phase not in Fluid_List:
                    Volume = Volume + float(Results[phase + '_prop']['v'].loc[i])
            
            Total_volume = Volume + float(Results['liquid1_prop']['v'].loc[i])

            if Volume/Total_volume > Crystallinity_limit:
                break

        melts = melts.addNodeAfter()

    return Results

def findSatPressure_MELTS(Model = None, T_C_init = None, P_bar_init = None, comp = None, fO2_buffer = None, fO2_offset = None):
    """
    Calculates the volatile saturation pressure and temperature of a given chemical composition (including volatiles) using the MELTS model.

    Parameters:
        Model (optional, default = None): a string that specifies the version of the MELTS model to be used.
            Options are "MELTSv1.0.2", "pMELTS", "MELTSv1.1.0", and "MELTSv1.2.0".
        T_C_init (optional, default = None): the initial temperature in degrees Celsius for the saturation pressure
            calculation. If not provided, the default value of 1200 °C is used.
        P_bar_init (optional, default = None): the initial pressure in bars for the saturation pressure calculation.
            If not provided, the default value of 10000 bars is used.
        comp: a dictionary of oxide names and their concentrations in the bulk. The oxide names must be in the format
            'OXIDE_Liq'.
        fO2_buffer (optional, default = None): a buffer for the fO2 value used in the MELTS model.
        fO2_offset (optional, default = None): an offset for the fO2 value used in the MELTS model.

    Returns:
        A dictionary with the following keys:
            'SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq',
            'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq': the liquid compositions of the respective elements.
            'P_bar': the saturation pressure in bars.
            'T_Liq': the temperature in degrees Celsius at which the saturation pressure is calculated.
    """

    if P_bar_init is None:
        P_bar_init = 10000

    P_bar = P_bar_init

    if T_C_init is None:
        T_C_init = 1200

    from meltsdynamic import MELTSdynamic

    if Model is None or Model == "MELTSv1.0.2":
        melts = MELTSdynamic(1)
    elif Model == "pMELTS":
        melts = MELTSdynamic(2)
    elif Model == "MELTSv1.1.0":
        melts = MELTSdynamic(3)
    elif Model == "MELTSv1.2.0":
        melts = MELTSdynamic(4)

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    T_Liq, H2O_Melt, melts = findLiq_MELTS(P_bar = P_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_C_init, bulk_return = True, Step = np.array([5,1]))

    if T_Liq == 0.0:
        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
        return out

    PhaseList = melts.engine.solidNames
    if PhaseList is None:
        PhaseList = ['liquid1']
    else:
        PhaseList = ['liquid1'] + PhaseList

    Steps = np.array([500, 50, 10])

    for i in range(len(Steps)):
        if "fluid1" in PhaseList:
            while "fluid1" in PhaseList:
                P_bar = P_bar + Steps[i]
                T_Liq, H2O_Melt, melts = findLiq_MELTS(P_bar = P_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_Liq, bulk_return = True, Step = np.array([5,1]))

                if T_Liq == 0.0:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                    return out
                    break

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList

        if "fluid1" not in PhaseList:
            while "fluid1" not in PhaseList:
                P_bar = P_bar - Steps[i]
                if P_bar < 0:
                    P_bar = P_bar + Steps[i]
                    break

                T_Liq, H2O_Melt, melts = findLiq_MELTS(P_bar = P_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_Liq, bulk_return = True, Step = np.array([5,1]))

                if T_Liq == 0.0:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                    return out
                    break

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList

    bulk_out = melts.engine.getProperty('dispComposition', 'liquid1')
    out = {'SiO2_Liq': bulk_out[0], 'TiO2_Liq': bulk_out[1], 'Al2O3_Liq': bulk_out[2], 'FeOt_Liq': bulk_out[5] + 71.844/(159.69/2)*bulk_out[3], 'MnO_Liq': bulk_out[6], 'MgO_Liq': bulk_out[7], 'CaO_Liq': bulk_out[10], 'Na2O_Liq': bulk_out[11], 'K2O_Liq': bulk_out[12], 'P2O5_Liq': bulk_out[13], 'H2O_Liq': bulk_out[14], 'CO2_Liq': bulk_out[15], 'Fe3Fet_Liq': (71.844/(159.69/2)*bulk_out[3])/(bulk_out[5] + 71.844/(159.69/2)*bulk_out[3]), 'P_bar': P_bar, 'T_Liq': T_Liq}

    return out









