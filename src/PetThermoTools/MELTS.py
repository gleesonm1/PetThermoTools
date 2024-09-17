import numpy as np
import pandas as pd
import sys
import time

def equilibrate_MELTS(Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    Results = {}
    Affinity = {}

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp.copy()
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

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
    melts.engine.setSystemProperties("Suppress", "tridymite")

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
    
    PhaseList = melts.engine.calcSaturationState()
    Affinity_raw = melts.engine.getProperty('affinity', PhaseList)
    # Affinity = {'orthopyroxene': Affinity_raw}
    Affinity = {Phase: float_value for Phase, float_value in zip(PhaseList, Affinity_raw)}

    return Results, Affinity

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

        melts.engine.setSystemProperties("Suppress", "rutile")
        melts.engine.setSystemProperties("Suppress", "tridymite")

    T_Liq = 0
    H2O = 0
    CO2 = 0

    Liq_Results = findLiq_MELTS(P_bar = P_bar, 
                                comp = bulk, 
                                melts = melts, 
                                Model = Model,
                                fO2_buffer = fO2_buffer, 
                                fO2_offset = fO2_offset, 
                                T_C_init = T_C,
                                step = np.array([5,1]))
    
    if Liq_Results['fluid_saturated'] == "No":
        CO2_step = np.array([0.1, 0.02, 0.005])
        for j in range(len(CO2_step)):
            while Liq_Results['fluid_saturated'] == "No":
                bulk[15] = bulk[15] + CO2_step[j]
                Liq_Results = findLiq_MELTS(P_bar = P_bar, 
                                            comp = bulk, 
                                            melts = melts, 
                                            Model = Model,
                                            fO2_buffer = fO2_buffer, 
                                            fO2_offset = fO2_offset, 
                                            T_C_init = Liq_Results['T_Liq'],
                                            step = np.array([0.1]))
                
            if j != len(CO2_step) - 1:
                bulk[15] = bulk[15] - CO2_step[j]
                Liq_Results = findLiq_MELTS(P_bar = P_bar, 
                                            comp = bulk, 
                                            melts = melts, 
                                            Model = Model,
                                            fO2_buffer = fO2_buffer, 
                                            fO2_offset = fO2_offset, 
                                            T_C_init = Liq_Results['T_Liq'],
                                            step = np.array([0.1]))
                
    T_Liq = Liq_Results['T_Liq']
    H2O = bulk[14]
    CO2 = bulk[15]
                

    # try:
    #     melts.engine.setBulkComposition(bulk)
    #     melts.engine.pressure = P_bar
    #     melts.engine.temperature = T_C
    # except:
    #     return T_Liq, H2O, CO2

    # if fO2_buffer is not None:
    #     if fO2_offset is None:
    #         melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
    #     else:
    #         melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    # Liq = ['liquid1','water1', 'fluid1']
    # try:
    #     melts.engine.calcEquilibriumState(1,0)
    # except:
    #     return T_Liq, H2O, CO2

    # PhaseList = melts.engine.solidNames
    # if PhaseList is None:
    #     PhaseList = ['liquid1']
    # else:
    #     PhaseList = ['liquid1'] + PhaseList

    # i = set.intersection(set(Liq),set(PhaseList))

    # if Step is None:
    #     Step = np.array([3,1,0.1])


    # for j in range(len(Step)):
    #     if len(i) == len(PhaseList):
    #         while len(i) == len(PhaseList):
    #             try:
    #                 melts = melts.addNodeAfter()
    #                 melts.engine.temperature = melts.engine.temperature - Step[j]
    #                 melts.engine.calcEquilibriumState(1,0)
    #             except:
    #                 return T_Liq, H2O, CO2

    #             PhaseList = melts.engine.solidNames
    #             if PhaseList is None:
    #                 PhaseList = ['liquid1']
    #             else:
    #                 PhaseList = ['liquid1'] + PhaseList
    #             i = set.intersection(set(Liq),set(PhaseList))

    #     if len(i) < len(PhaseList):
    #         while len(i) < len(PhaseList):
    #             try:
    #                 melts = melts.addNodeAfter()
    #                 melts.engine.temperature = melts.engine.temperature + Step[j]
    #                 melts.engine.calcEquilibriumState(1,0)
    #             except:
    #                 return T_Liq, H2O, CO2

    #             PhaseList = melts.engine.solidNames
    #             if PhaseList is None:
    #                 PhaseList = ['liquid1']
    #             else:
    #                 PhaseList = ['liquid1'] + PhaseList
    #             i = set.intersection(set(Liq),set(PhaseList))


    # if "fluid1" not in PhaseList:
    #     CO2_step = np.array([0.1, 0.02, 0.005])
    #     for j in range(len(CO2_step)):
    #         while "fluid1" not in PhaseList:
    #             bulk[15] = bulk[15] + CO2_step[j]
    #             melts.engine.setBulkComposition(bulk)
    #             if len(i) == len(PhaseList):
    #                 while len(i) == len(PhaseList):
    #                     try:
    #                         melts = melts.addNodeAfter()
    #                         melts.engine.temperature = melts.engine.temperature - 0.1
    #                         melts.engine.calcEquilibriumState(1,0)
    #                     except:
    #                         return T_Liq, H2O, CO2

    #                     PhaseList = melts.engine.solidNames
    #                     if PhaseList is None:
    #                         PhaseList = ['liquid1']
    #                     else:
    #                         PhaseList = ['liquid1'] + PhaseList
    #                     i = set.intersection(set(Liq),set(PhaseList))

    #             if len(i) < len(PhaseList):
    #                 while len(i) < len(PhaseList):
    #                     try:
    #                         melts = melts.addNodeAfter()
    #                         melts.engine.temperature = melts.engine.temperature + 0.1
    #                         melts.engine.calcEquilibriumState(1,0)
    #                     except:
    #                         return T_Liq, H2O, CO2

    #                     PhaseList = melts.engine.solidNames
    #                     if PhaseList is None:
    #                         PhaseList = ['liquid1']
    #                     else:
    #                         PhaseList = ['liquid1'] + PhaseList
    #                     i = set.intersection(set(Liq),set(PhaseList))

    #         if j != len(CO2_step) - 1:
    #             bulk[15] = bulk[15] - CO2_step[j]

    # T_Liq = melts.engine.getProperty('temperature')
    # H2O = melts.engine.getProperty('dispComposition', 'liquid1', 'H2O')
    # CO2 = melts.engine.getProperty('dispComposition', 'liquid1', 'CO2')

    return T_Liq, H2O, CO2

def findLiq_MELTS(P_bar = None, Model = None, T_C_init = None, comp = None, melts = None, fO2_buffer = None, fO2_offset = None, Step = None, fluid_test = None, bulk_return = None, step = None, Affinity = False):
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
    Result: dictionary
        liquidus temperature, phase, and melt composition.

    Affinity: dictionary (optional, select Affinity = True as a kwarg)
        Affinity (J/mol) of phases not present at the liquidus.

    '''
    Results = {}
    Affin = {}
    if P_bar is None:
        raise Exception("Please specify a pressure for calculations")

    if Model is None:
        if melts is None:
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

    # pd.DataFrame(data = np.zeros((1, 17)), columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])

    # T_Liq = 0
    # H2O_Melt = 0
    # CO2_Melt = 0


    if melts is None:
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
        melts.engine.setSystemProperties("Suppress", "tridymite")
        
    if fO2_buffer is not None:
        if fO2_offset is None:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
        else:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    try:
        melts.engine.setBulkComposition(bulk)
        melts.engine.pressure = P_bar
        melts.engine.temperature = T_C_init
    except:
        return Results
    
    if Model == "MELTSv1.2.2":
        try:
            melts.engine.findLiquidus()
        except:
            if Affinity is True:
                if bulk_return is not None:
                    return Results, Affin, melts
                else:
                    return Results, Affin
            else:
                if bulk_return is not None:
                    return Results, melts
                else:
                    return Results
            
        PL = melts.engine.solidNames
        PhaseList = ['liquid1'] + PL

    else:
        Liq = ['liquid1','water1', 'fluid1']
        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            if Affinity is True:
                if bulk_return is not None:
                    return Results, Affin, melts
                else:
                    return Results, Affin
            else:
                if bulk_return is not None:
                    return Results, melts
                else:
                    return Results

        PhaseList = melts.engine.solidNames
        if PhaseList is None:
            PhaseList = ['liquid1']
        else:
            PhaseList = ['liquid1'] + PhaseList

        i = set.intersection(set(Liq),set(PhaseList))

        if Step is None:
            Step = np.array([3,1,0.1])

        for j in range(len(Step)):
            if len(i) < len(PhaseList):
                while len(i) < len(PhaseList):
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.temperature = melts.engine.temperature + Step[j]
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        if Affinity is True:
                            if bulk_return is not None:
                                return Results, Affin, melts
                            else:
                                return Results, Affin
                        else:
                            if bulk_return is not None:
                                return Results, melts
                            else:
                                return Results

                        break

                    PhaseList = melts.engine.solidNames
                    if PhaseList is None:
                        PhaseList = ['liquid1']
                    else:
                        PhaseList = ['liquid1'] + PhaseList
                    i = set.intersection(set(Liq),set(PhaseList))
        
            if len(i) == len(PhaseList):
                while len(i) == len(PhaseList):
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.temperature = melts.engine.temperature - Step[j]
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        if Affinity is True:
                            if bulk_return is not None:
                                return Results, Affin, melts
                            else:
                                return Results, Affin
                        else:
                            if bulk_return is not None:
                                return Results, melts
                            else:
                                return Results

                        break

                    PhaseList = melts.engine.solidNames
                    if PhaseList is None:
                        PhaseList = ['liquid1']
                    else:
                        PhaseList = ['liquid1'] + PhaseList
                    i = set.intersection(set(Liq),set(PhaseList))

        PL = melts.engine.solidNames
        try:
            melts = melts.addNodeAfter()
            melts.engine.temperature = melts.engine.temperature + 4*Step[j]
            melts.engine.calcEquilibriumState(1,0)

            
            melts = melts.addNodeAfter()
            melts.engine.temperature = melts.engine.temperature - 3*Step[j]
            melts.engine.calcEquilibriumState(1,0)
        except:
            if Affinity is True:
                if bulk_return is not None:
                    return Results, Affin, melts
                else:
                    return Results, Affin
            else:
                if bulk_return is not None:
                    return Results, melts
                else:
                    return Results

    columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
    # pd.DataFrame(data = np.zeros((1, 17)), columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    for j in columns:
        if j == "T_Liq":
            Results['T_Liq'] = melts.engine.temperature            
        elif j == "fluid_saturated":
            if type(melts.engine.solidNames) == list:
                if Model == "MELTSv1.0.2":
                    if 'water1' in PhaseList:
                        Results['fluid_saturated'] = "Yes"
                    else:
                        Results['fluid_saturated'] = "No"
                if 'water1' not in PhaseList:
                    if 'fluid1' in PhaseList:
                        Results['fluid_saturated'] = "Yes"
                    else:
                        Results['fluid_saturated'] = "No"
        elif j == "liquidus_phase":
            if type(PL) == list:
                if 'fluid1' in PL:
                    PL.remove('fluid1')
                    Results['liquidus_phase'] = PL[0]
                    
                elif 'water1' in PL:
                    PL.remove('water1')
                    Results['liquidus_phase'] = PL[0]

                else:
                    Results['liquidus_phase'] = PL[0]  
            else:
                Results['liquidus_phase'] = PL[0]
        else:
            Results[j] = melts.engine.getProperty('dispComposition', 'liquid1', j)

    if Affinity is True:
        PhaseList = melts.engine.calcSaturationState()
        Affinity_raw = melts.engine.getProperty('affinity', PhaseList)
        Affin = {Phase: float_value for Phase, float_value in zip(PhaseList, Affinity_raw)}

        if bulk_return is not None:
            return Results, Affin, melts
        else:
            return Results, Affin
    else:
        if bulk_return is not None:
            return Results, melts
        else:
            return Results

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
    Results = {phases[0]: np.nan, phases[1]: np.nan, phases[2]: np.nan, 'T_Liq': np.nan, 'H2O_melt': np.nan}
    if len(phases) == 2:
        del Results[phases[2]]

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
    melts.engine.setSystemProperties("Suppress", "tridymite")

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    try:
        Liq_Results = findLiq_MELTS(P_bar = P_bar, comp = bulk, T_C_init = T_initial_C, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        Results['T_Liq'] = Liq_Results['T_Liq']
        Results['H2O_melt'] = Liq_Results['H2O']
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

                if phases[0] in PhaseList and np.isnan(Results[phases[0]]):# == 0:
                    Results[phases[0]] = melts.engine.temperature

                if phases[1] in PhaseList and np.isnan(Results[phases[1]]):# == 0:
                    Results[phases[1]] = melts.engine.temperature

                if len(phases) == 3:
                    if phases[2] in PhaseList and np.isnan(Results[phases[2]]):# == 0:
                        Results[phases[2]] = melts.engine.temperature

                    if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]) and ~np.isnan(Results[phases[2]]):# > 0:
                        break

                if len(phases) == 2:
                    if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]):# > 0:
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
    melts.engine.setSystemProperties("Suppress", "tridymite")

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    if find_liquidus is not None:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    Liq_Results = findLiq_MELTS(P_bar = P_path_bar[0], comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
                else:
                    Liq_Results = findLiq_MELTS(P_bar = P_path_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
            except:
                return Results
        elif P_start_bar is not None:
            try:
                Liq_Results = findLiq_MELTS(P_bar = P_start_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = 1400)
            except:
                return Results

        T_start_C = Liq_Results['T_Liq'] + 0.1

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
    if type(T) == np.ndarray and P_end_bar is not None and dp_bar is None:
        P = np.linspace(P_start_bar, P_end_bar, len(T))
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
        melts.engine.temperature = T_start_C

        if P_path_bar is not None or T_path_C is not None:
            if type(P) == np.ndarray and type(T) == np.ndarray:
                T_Liq_loc = np.abs(T - T_start_C).argmin()
                if T[T_Liq_loc]>T_start_C:
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

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 8)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'mass', 'dvdp', 'logfO2'])
    Results['liquid1'] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    Results['liquid1_prop'] = pd.DataFrame(data = np.zeros((length, 4)), columns = ['h', 'mass', 'v', 'rho'])

    for i in range(length):
        if type(T) == np.ndarray:
            melts.engine.temperature = T[i]
        if type(P) == np.ndarray:
            if i == 0:
                melts.engine.pressure = P[i]
            if isochoric is None:
                melts.engine.pressure = P[i]

        if Frac_solid is not None:
            melts.engine.setSystemProperties("Mode", "Fractionate Solids")

        if Frac_fluid is not None:
            melts.engine.setSystemProperties("Mode", "Fractionate Fluids")

        if isenthalpic is not None:
            if i == 0:
                try:
                    melts.engine.setSystemProperties("Mode", "Isenthalpic")
                    melts.engine.calcEquilibriumState(1,0)
                    melts.engine.setSystemProperties("Mode", "Isenthalpic")
                except:
                    break

        if isentropic is not None:
            if i == 0:
                try:
                    melts.engine.setSystemProperties("Mode", "Isentropic")
                    melts.engine.calcEquilibriumState(1,0)
                    melts.engine.setSystemProperties("Mode", "Isentropic")
                except:
                    break

        if isochoric is not None:
            if i == 0:
                try:
                    melts.engine.setSystemProperties("Mode", "Isochoric")
                    melts.engine.calcEquilibriumState(1,0)
                    melts.engine.setSystemProperties("Mode", "Isochoric")
                except:
                    break

        if isochoric is None and isenthalpic is None and isentropic is None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(1,1)
                else:
                    melts.engine.calcEquilibriumState(1,0)
            except:
                return Results
                break

        if isenthalpic is not None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(2,1)
                else:
                    melts.engine.calcEquilibriumState(2,0)
            except:
                # return Results
                break

        if isentropic is not None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(3,1)
                else:
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
            return Results
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

def findSatPressure_MELTS_multi(Model = None, comp = None, fO2_buffer = None, fO2_offset = None, P_bar = None, T_fixed_C = None):
    out = pd.DataFrame(columns = ['Sample_ID_Liq', 'SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq', 'P_bar', 'T_Liq'])

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
    melts.engine.setSystemProperties("Suppress", "tridymite")

    if "Sample_ID_Liq" not in comp.keys():
        comp.loc[:,"Sample_ID_Liq"] = np.linspace(0, len(comp['SiO2_Liq'])-1, len(comp['SiO2_Liq']))

    for i in range(len(T_fixed_C)):
        bulk = [comp.loc[i,'SiO2_Liq'], comp.loc[i,'TiO2_Liq'], comp.loc[i,'Al2O3_Liq'], 
                comp.loc[i,'Fe3Fet_Liq']*((159.59/2)/71.844)*comp.loc[i,'FeOt_Liq'], 
                comp.loc[i,'Cr2O3_Liq'], (1- comp.loc[i,'Fe3Fet_Liq'])*comp.loc[i,'FeOt_Liq'], 
                comp.loc[i,'MnO_Liq'], comp.loc[i,'MgO_Liq'], 0.0, 0.0, comp.loc[i,'CaO_Liq'], 
                comp.loc[i,'Na2O_Liq'], comp.loc[i,'K2O_Liq'], comp.loc[i,'P2O5_Liq'], 
                comp.loc[i,'H2O_Liq'], comp.loc[i,'CO2_Liq'], 0.0, 0.0, 0.0]
        bulk = list(100*np.array(bulk)/np.sum(bulk))
        if i == 0:
            melts.engine.pressure = np.random.normal(P_bar[i], P_bar[i]/10)
            melts.engine.temperature = T_fixed_C[i] + 500
            melts.engine.setBulkComposition(bulk)
            PL = melts.engine.calcSaturationState()
            for p in PL:
                if p != "fluid":
                    if p != "water":
                        melts.engine.setSystemProperties("Suppress", p)

        if fO2_buffer is not None:
            if fO2_offset[i] is None:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
            else:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset[i])])

        melts = melts.addNodeAfter()
        melts.engine.pressure = P_bar[i]
        melts.engine.temperature = T_fixed_C[i]
        melts.engine.setBulkComposition(bulk)
        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            return out

        for j in range(2):
            if len(melts.engine.solidNames) > 0:
                while len(melts.engine.solidNames) > 0:
                    P_low = P_bar[i]
                    if P_bar[i] > 5000:
                        P_bar[i] = P_bar[i]*1.1
                    else:
                        P_bar[i] = P_bar[i]*1.5
                    melts.engine.pressure = P_bar[i]
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        return out
                    
                    if len(melts.engine.solidNames) == 0:
                        P_high = P_bar[i]
                
            else:
                while len(melts.engine.solidNames) == 0:
                    P_high = P_bar[i]
                    P_bar[i] = P_bar[i]/2
                    melts.engine.pressure = P_bar[i]
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        return out
                    
                    if len(melts.engine.solidNames) == 0:
                        P_low = P_bar[i]

        for j in range(5):
            P_bar[i] = P_low + (P_high - P_low)/2
            melts.engine.pressure = P_bar[i]
            try:
                melts = melts.addNodeAfter()
                melts.engine.calcEquilibriumState(1,0)
            except:
                return out
            
            if len(melts.engine.solidNames) > 0:
                P_low = P_bar[i]
            else:
                P_high = P_bar[i]

            if P_high/P_low < 1.002:
                break

        res = {'SiO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'SiO2'), 
            'TiO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'TiO2'),
            'Al2O3_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'Al2O3'), 
            'Cr2O3_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'Cr2O3'),
            'FeOt_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'FeO') + 71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3'), 
            'MnO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'MnO'), 
            'MgO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'MgO'), 
            'CaO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'CaO'),
            'Na2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'Na2O'), 
            'K2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'K2O'), 
            'P2O5_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'P2O5'), 
            'H2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'H2O'), 
            'CO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'CO2'), 
            'Fe3Fet_Liq': (71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3'))/(melts.engine.getProperty('dispComposition', 'liquid1', 'FeO') + 71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3')), 
            'P_bar': P_bar[i], 
            'T_Liq': T_fixed_C[i],
            'Sample_ID_Liq': comp.loc[i, "Sample_ID_Liq"]}
        
        out.loc[i,:] = res

    return out

def findSatPressure_MELTS(Model = None, T_C_init = None, T_fixed_C = None, P_bar_init = None, comp = None, fO2_buffer = None, fO2_offset = None):
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
        P_bar_init = 5000

    P_bar = np.random.normal(P_bar_init, P_bar_init/10)

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

    melts.engine.setSystemProperties("Suppress", "rutile")
    melts.engine.setSystemProperties("Suppress", "tridymite")

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    if T_fixed_C is not None:
        melts.engine.pressure = P_bar_init
        melts.engine.temperature = T_fixed_C + 500
        melts.engine.setBulkComposition(bulk)
        PL = melts.engine.calcSaturationState()
        for p in PL:
            if p != "fluid":
                if p != "water":
                    melts.engine.setSystemProperties("Suppress", p)

        melts.engine.pressure = P_bar_init
        melts.engine.temperature = T_fixed_C
        melts.engine.setBulkComposition(bulk)
        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
            return out
        
        for i in range(2):
            if len(melts.engine.solidNames) > 0:
                while len(melts.engine.solidNames) > 0:
                    P_low = P_bar
                    if P_bar > 5000:
                        P_bar = P_bar*1.1
                    else:
                        P_bar = P_bar*1.5
                    melts.engine.pressure = P_bar
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                        return out
                    
                    if len(melts.engine.solidNames) == 0:
                        P_high = P_bar
                
            else:
                while len(melts.engine.solidNames) == 0:
                    P_high = P_bar
                    P_bar = P_bar/2
                    melts.engine.pressure = P_bar
                    try:
                        melts = melts.addNodeAfter()
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                        return out
                    
                    if len(melts.engine.solidNames) == 0:
                        P_low = P_bar

        for i in range(5):
            P_bar = P_low + (P_high - P_low)/2
            melts.engine.pressure = P_bar
            try:
                melts = melts.addNodeAfter()
                melts.engine.calcEquilibriumState(1,0)
            except:
                out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                return out
            
            if len(melts.engine.solidNames) > 0:
                P_low = P_bar
            else:
                P_high = P_bar

            if P_high/P_low < 1.002:
                break

        out = {'SiO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'SiO2'), 
            'TiO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'TiO2'),
            'Al2O3_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'Al2O3'), 
            'FeOt_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'FeO') + 71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3'), 
            'MnO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'MnO'), 
            'MgO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'MgO'), 
            'CaO_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'CaO'),
            'Na2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'Na2O'), 
            'K2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'K2O'), 
            'P2O5_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'P2O5'), 
            'H2O_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'H2O'), 
            'CO2_Liq': melts.engine.getProperty('dispComposition', 'liquid1', 'CO2'), 
            'Fe3Fet_Liq': (71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3'))/(melts.engine.getProperty('dispComposition', 'liquid1', 'FeO') + 71.844/(159.69/2)*melts.engine.getProperty('dispComposition', 'liquid1', 'Fe2O3')), 
            'P_bar': P_bar, 
            'T_Liq': T_fixed_C}

    else:    
        Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                    fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_C_init, Step = np.array([5,1]))
        
        if Liq_Results['T_Liq'] == 0.0:
            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
            return out
        
        if Liq_Results['fluid_saturated'] == "Yes":
            while Liq_Results['fluid_saturated'] == "Yes":
                P_low = P_bar
                P_bar = P_bar*2
                Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq'], Step = np.array([5,1]))

                if Liq_Results['T_Liq'] == 0.0:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                    return out
                
                if Liq_Results['fluid_saturated'] == "No":
                    P_high = P_bar
                    break
        
        else:
            while Liq_Results['fluid_saturated'] == "No":
                P_high = P_bar
                P_bar = P_bar/2
                Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq'], Step = np.array([5,1]))

                if Liq_Results['T_Liq'] == 0.0:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                    return out
                
                if Liq_Results['fluid_saturated'] == "Yes":
                    P_low = P_bar
                    break

        while P_high/P_low > 1.01:
            P_new = P_low + (P_high - P_low)/2
            Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_new, comp = bulk, melts = melts, 
                                    fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq'], Step = np.array([3,1]))
            
            if Liq_Results['T_Liq'] == 0.0:
                out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq': 0.0}
                return out
                break

            if Liq_Results['fluid_saturated'] == "Yes":
                P_low = P_new
            else:
                P_high = P_new

            if P_high/P_low < 1.01:
                break

        P_bar = round((P_high+P_low)/2)
        Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                                        T_C_init = Liq_Results['T_Liq'], Step = np.array([1]))
        
        out = {'SiO2_Liq': Liq_Results['SiO2'], 
            'TiO2_Liq': Liq_Results['TiO2'],
            'Al2O3_Liq': Liq_Results['Al2O3'], 
            'FeOt_Liq': Liq_Results['FeO'] + 71.844/(159.69/2)*Liq_Results['Fe2O3'], 
            'MnO_Liq': Liq_Results['MnO'], 
            'MgO_Liq': Liq_Results['MgO'], 
            'CaO_Liq': Liq_Results['CaO'],
            'Na2O_Liq': Liq_Results['Na2O'], 
            'K2O_Liq': Liq_Results['K2O'], 
            'P2O5_Liq': Liq_Results['P2O5'], 
            'H2O_Liq': Liq_Results['H2O'], 
            'CO2_Liq': Liq_Results['CO2'], 
            'Fe3Fet_Liq': (71.844/(159.69/2)*Liq_Results['Fe2O3'])/(Liq_Results['FeO'] + 71.844/(159.69/2)*Liq_Results['Fe2O3']), 
            'P_bar': P_bar, 
            'T_Liq': Liq_Results['T_Liq']}

    return out

def AdiabaticDecompressionMelting_MELTS(Model = None, comp = None, Tp_C = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Frac = False, fO2_buffer = None, fO2_offset = None):
    try:
        import pyMelt as m     
        Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
                    'KG1': m.lithologies.matthews.kg1(),
                    'G2': m.lithologies.matthews.eclogite(),
                    'hz': m.lithologies.shorttle.harzburgite()}
    except ImportError:
        raise RuntimeError('You havent installed pyMelt or there is an error when importing pyMelt. pyMelt is currently required to estimate the starting point for the melting calculations.')

    Results = {}

    if comp is None:
        raise Exception("No composition specified")

    if P_path_bar is None and P_start_bar is None:
        raise Exception("Initial P system must be defined")

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
    melts.engine.setSystemProperties("Suppress", "tridymite")
    
    lz = m.lithologies.matthews.klb1()
    mantle = m.mantle([lz], [1], ['Lz'])
    if P_start_bar is not None:
        T_start_C = mantle.adiabat(P_start_bar/10000, Tp_C)
    else:
        T_start_C = mantle.adiabat(P_path_bar[0]/10000, Tp_C)

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    melts.engine.setBulkComposition(bulk)
    melts.engine.temperature = T_start_C

    if P_path_bar is not None:
        P = P_path_bar
    else:
        P = np.linspace(P_start_bar, P_end_bar, 1+int(round((P_start_bar - P_end_bar)/dp_bar)))

    Results['Conditions'] = pd.DataFrame(data = np.zeros((len(P), 7)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'dvdp', 'logfO2'])

    for k in range(len(P)):
        if k != 0:
            melts = melts.addNodeAfter()

        melts.engine.pressure = P[k]
        if k == 0:
            T = T_start_C

            try:
                melts.engine.calcEquilibriumState(1,0)
                s = melts.engine.getProperty('s','bulk')
            except:
                return Results
            # try:
            #     melts.engine.setSystemProperties("Mode", "Isentropic")
            #     melts.engine.calcEquilibriumState(1,0)
            #     melts.engine.setSystemProperties("Mode", "Isentropic")
            # except:
            #     return Results
        if k > 0:
            s_check = 0
            n = 4
            while np.abs(s_check-s)/s > 0.0002:
                T_save = np.linspace(T, T - 3.0, n)
                s_save = np.zeros(len(T_save))
                for i in range(len(T_save)):
                    melts.engine.temperature = T_save[i]
                    try:
                        melts.engine.calcEquilibriumState(1,0)
                    except:
                        return Results

                    s_save[i] = melts.engine.getProperty('s','bulk')
                    melts = melts.addNodeAfter()

                    # print(s_save[i])

                p = np.polyfit(s_save, T_save, 3)
                T_next = p[0]*s**3 + p[1]*s**2 + p[2]*s + p[3]

                melts.engine.temperature = T_next

                try:
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    return Results
                
                s_check = melts.engine.getProperty('s','bulk')

                n = n*2
                if n > 32:
                    break

            T = T_next
            
        for R in Results['Conditions']:
            if R == 'temperature':
                Results['Conditions'][R].loc[k] = melts.engine.temperature
            elif R == 'pressure':
                Results['Conditions'][R].loc[k] = melts.engine.pressure
            elif R == 'logfO2':
                try:
                    Results['Conditions'][R].loc[k] = melts.engine.getProperty(R)
                except:
                    Results['Conditions'][R].loc[k] = np.nan
            else:
                Results['Conditions'][R].loc[k] = melts.engine.getProperty(R, 'bulk')

        try:
            if ~np.isnan(melts.engine.getProperty('mass', 'liquid1')):
                PhaseList = ['liquid1'] + melts.engine.solidNames
            else:
                PhaseList = melts.engine.solidNames
        except:
            return Results

        for phase in PhaseList:
            if phase not in list(Results.keys()):
                Results[phase] = pd.DataFrame(data = np.zeros((len(P), 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
                Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((len(P), 4)), columns = ['h', 'mass', 'v', 'rho'])

            if phase in list(Results.keys()):
                for el in Results[phase]:
                    Results[phase][el].loc[k] = melts.engine.getProperty('dispComposition', phase, el)

                for pr in Results[phase + '_prop']:
                    Results[phase + '_prop'][pr].loc[k] = melts.engine.getProperty(pr, phase)
        
        # if k == 0:
        #     s = melts.engine.getProperty('s','bulk')

    return Results







