import numpy as np
import pandas as pd
import sys
import time
import copy
from petthermotools.GenFuncs import estimate_t

def equilibrate_MELTS(Model = None, P_bar = None, T_C = None, comp = None, 
                      fO2_buffer = None, fO2_offset = None, Suppress = None, Suppress_except=None,
                      trail = None, melts = None):
    Results = {}
    Affinity = {}



    if trail is not None:
        trail = False

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp.copy()
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    bulk = list(100*np.array(bulk)/np.sum(bulk))

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

    melts.engine.setBulkComposition(bulk)
    melts.engine.pressure = P_bar
    melts.engine.temperature = T_C

    if Suppress is None:
        melts.engine.setSystemProperties("Suppress", "rutile")
        melts.engine.setSystemProperties("Suppress", "tridymite")
    else:
        if Suppress == "All":
            melts.engine.pressure = P_bar # np.random.normal(P_bar, P_bar/10)
            melts.engine.temperature = T_C + 500
            PL = melts.engine.calcSaturationState()
            for p in PL:
                if p != "fluid":
                    if p != "water":
                        melts.engine.setSystemProperties("Suppress", p)
        elif Suppress_except is not None:
            if type(Suppress_except) == str:
                Suppress_except = [Suppress_except]
            melts.engine.pressure = P_bar # np.random.normal(P_bar, P_bar/10)
            melts.engine.temperature = T_C + 500
            PL = melts.engine.calcSaturationState()
            for p in PL:
                if p not in Suppress_except:
                    melts.engine.setSystemProperties("Suppress", p)   
        else:
            for p in Suppress:
                melts.engine.setSystemProperties("Suppress", p)

        melts = melts.addNodeAfter()

        melts.engine.pressure = P_bar
        melts.engine.temperature = T_C

    length = 1

    if fO2_buffer is not None:
        if fO2_offset is None:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
        else:
            melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 8)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'mass', 'dvdp', 'logfO2'])
    properties = ['g', 'h', 's', 'v', 'cp', 'dcpdt', 'dvdt', 'dpdt', 'd2vdt2', 'd2vdtdp', 'd2vdp2', 'molwt', 'rho', 'mass']

    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        if trail is not None:
            return Results, Affinity, trail
        else:
            return Results, Affinity
    
    if melts.engine.getProperty('mass', 'liquid1') > 0.0:
        PhaseList = ['liquid1'] + melts.engine.solidNames
    else:
        PhaseList = melts.engine.solidNames

    for R in Results['Conditions']:
        if R == 'temperature':
            Results['Conditions'].loc[0,R] = melts.engine.temperature
        elif R == 'pressure':
            Results['Conditions'].loc[0,R] = melts.engine.pressure
        elif R == 'logfO2':
            try:
                Results['Conditions'].loc[0,R] = melts.engine.getProperty(R)
            except:
                Results['Conditions'].loc[0,R] = np.nan
        else:
            Results['Conditions'].loc[0,R] = melts.engine.getProperty(R, 'bulk')

    for phase in PhaseList:
        if phase not in list(Results.keys()):
            thermo_properties = []
            endmembers = melts.endMemberFormulas[phase[:-1]]
            thermo = ['activity', 'activity0', 'mu','mu0', 'X']
            for iii in thermo:
                for jjj in endmembers:
                    thermo_properties.append(iii+'_'+jjj)

            Results[phase] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
            Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((length, len(properties) + len(thermo_properties))), columns = properties + thermo_properties) #['h', 'mass', 'v', 'rho'])

        if phase in list(Results.keys()):
            for el in Results[phase]:
                Results[phase].loc[0,el] = melts.engine.getProperty('dispComposition', phase, el)

            melts.engine.calcEndMemberProperties(phase, melts.engine.getProperty('dispComposition', phase))

            for pr in Results[phase + '_prop']:
                if pr in properties:
                    Results[phase + '_prop'].loc[0,pr] = melts.engine.getProperty(pr, phase)
                else:
                    Results[phase + '_prop'].loc[0, pr] = melts.engine.getProperty(pr.split('_')[0], phase, pr.split('_')[1])
        # melts.engine.calcPhaseProperties(phase, melts.engine.dispComposition[phase])
        # melts.engine.calcEndMemberProperties(phase, melts.engine.dispComposition[phase])

        # endmembers = melts.endMemberFormulas[phase]
        # thermo = ['activity', 'activity0', 'mu','mu0', 'X']
        # thermo_properties = []
        # thermodynamics = np.empty(len(endmembers)*len(thermo))
        # idx = 0
        # for i in thermo:
        #     for j in endmembers:
        #         thermo_properties.append(i+'_'+j)
        #         thermodynamics[idx] = melts.engine.getProperty(i, phase, j)
        #         idx = idx + 1
    
    PhaseList = melts.engine.calcSaturationState()
    Affinity_raw = melts.engine.getProperty('affinity', PhaseList)
    Affinity = {Phase: float_value for Phase, float_value in zip(PhaseList, Affinity_raw)}
    # Results['Affinity'] = pd.DataFrame(Affinity)

    if trail is not None:
        return Results, Affinity, trail
    else:
        return Results, Affinity

def findCO2_MELTS(P_bar = None, Model = None, T_C = None, comp = None, melts = None, fO2_buffer = None, fO2_offset = None, Step = None):

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp.copy()
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    bulk = list(100*np.array(bulk)/np.sum(bulk))

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

    T_Liq_C = 0
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
                                            T_C_init = Liq_Results['T_Liq_C'],
                                            step = np.array([0.1]))
                
            if j != len(CO2_step) - 1:
                bulk[15] = bulk[15] - CO2_step[j]
                Liq_Results = findLiq_MELTS(P_bar = P_bar, 
                                            comp = bulk, 
                                            melts = melts, 
                                            Model = Model,
                                            fO2_buffer = fO2_buffer, 
                                            fO2_offset = fO2_offset, 
                                            T_C_init = Liq_Results['T_Liq_C'],
                                            step = np.array([0.1]))
                
    T_Liq_C = Liq_Results['T_Liq_C']
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

    return T_Liq_C, H2O, CO2

def findLiq_MELTS(P_bar = None, Model = None, T_C_init = None, comp = None, melts = None, 
                fO2_buffer = None, fO2_offset = None, Step = None, fluid_test = None, 
                bulk_return = None, step = None, Affinity = False):
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

    # T_C_init = np.round(np.random.normal(T_C_init, T_C_init/100))

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    
    bulk = list(100*np.array(bulk)/np.sum(bulk))
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
        
        melts.engine.setSystemProperties("Output",  "none")

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
        # try:
        #     melts = melts.addNodeAfter()
        #     melts.engine.temperature = melts.engine.temperature + 4*Step[j]
        #     melts.engine.calcEquilibriumState(1,0)

            
        #     melts = melts.addNodeAfter()
        #     melts.engine.temperature = melts.engine.temperature - 3*Step[j]
        #     melts.engine.calcEquilibriumState(1,0)
        # except:
        #     if Affinity is True:
        #         if bulk_return is not None:
        #             return Results, Affin, melts
        #         else:
        #             return Results, Affin
        #     else:
        #         if bulk_return is not None:
        #             return Results, melts
        #         else:
        #             return Results

    columns = ['T_Liq_C', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
    # pd.DataFrame(data = np.zeros((1, 17)), columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    for j in columns:
        if j == "T_Liq_C":
            Results['T_Liq_C'] = melts.engine.temperature            
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

        
def supCalc_MELTS(Model = "MELTSv1.0.2", comp = None, phase = None, T_C = None, P_bar = None,
             fO2_buffer = None, fO2_offset = None, 
             melts = None):

    if type(comp) == pd.core.frame.DataFrame: # simplest scenario - one calculation per bulk composition imported
        L = len(comp['SiO2_Liq'])
    else:
        if type(P_bar) == np.ndarray: #one calculation per P loaded
            L = len(P_bar)
        elif type(T_C) == np.ndarray: # one calculation per T loaded.
            L = len(T_C)
        elif type(fO2_offset) == np.ndarray or type(fO2_offset) == list: # one calculation per fO2 loaded.
            L = len(fO2_offset)
        else:
            L = -1

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

    if L == -1:
        melts.engine.temperature = T_C
        melts.engine.pressure = P_bar

        if fO2_buffer is not None:
            if fO2_offset is None:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
            else:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

        bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

        melts.engine.calcPhaseProperties(phase, bulk)
        melts.engine.calcEndMemberProperties(phase, melts.engine.dispComposition[phase])

        properties = ['g', 'h', 's', 'v', 'cp', 'dcpdt',
              'dvdt', 'dpdt', 'd2vdt2', 'd2vdtdp',
              'd2vdp2', 'molwt', 'rho', 'mass']
        results = np.empty(len(properties))  # Creates an array of the same length as the list
        # Populate the NumPy array
        for idx, prop in enumerate(properties):
            results[idx] = melts.engine.getProperty(prop, phase)

        oxides = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
        conc = np.empty(len(oxides))
        for idx, i in enumerate(oxides):
            conc[idx] = melts.engine.getProperty('dispComposition', phase, i)

        endmembers = melts.endMemberFormulas[phase]
        thermo = ['activity', 'activity0', 'mu','mu0', 'X']
        thermo_properties = []
        thermodynamics = np.empty(len(endmembers)*len(thermo))
        idx = 0
        for i in thermo:
            for j in endmembers:
                thermo_properties.append(i+'_'+j)
                thermodynamics[idx] = melts.engine.getProperty(i, phase, j)
                idx = idx + 1

        Output = {'Properties': pd.Series(index = properties, data = results), 
                'Composition': pd.Series(index = oxides, data = conc),
                'Thermodynamics': pd.Series(index = thermo_properties, data = thermodynamics)}
    else:
        if type(T_C) == float or type(T_C) == int:
            T_C = np.zeros(L) + T_C
        if type(P_bar) == float or type(P_bar) == int:
            P_bar = np.zeros(L) + P_bar
        if fO2_offset is None:
            fO2_offset = [None]*L
        
        properties = ['g', 'h', 's', 'v', 'cp', 'dcpdt',
                'dvdt', 'dpdt', 'd2vdt2', 'd2vdtdp',
                'd2vdp2', 'molwt', 'rho', 'mass']
        results = np.empty((L,len(properties)))
        
        oxides = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
        conc = np.empty((L,len(oxides)))

        endmembers = None

        for l in range(L):
            melts.engine.temperature = T_C[l]
            melts.engine.pressure = P_bar[l]

            if fO2_buffer is not None:
                if fO2_offset[l] is None:
                    melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
                else:
                    melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset[l])])

            if type(comp) == dict:
                bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1 - comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
            else:
                bulk = [comp.loc[l,'SiO2_Liq'], comp.loc[l,'TiO2_Liq'], comp.loc[l,'Al2O3_Liq'], comp.loc[l,'Fe3Fet_Liq']*((159.59/2)/71.844)*comp.loc[l,'FeOt_Liq'], comp.loc[l,'Cr2O3_Liq'], (1 - comp.loc[l,'Fe3Fet_Liq'])*comp.loc[l,'FeOt_Liq'], comp.loc[l,'MnO_Liq'], comp.loc[l,'MgO_Liq'], 0.0, 0.0, comp.loc[l,'CaO_Liq'], comp.loc[l,'Na2O_Liq'], comp.loc[l,'K2O_Liq'], comp.loc[l,'P2O5_Liq'], comp.loc[l,'H2O_Liq'], comp.loc[l,'CO2_Liq'], 0.0, 0.0, 0.0]

            melts.engine.calcPhaseProperties(phase, bulk)
            melts.engine.calcEndMemberProperties(phase, melts.engine.dispComposition[phase])

            # Populate the NumPy array
            for idx, prop in enumerate(properties):
                results[l,idx] = melts.engine.getProperty(prop, phase)

            for idx, i in enumerate(oxides):
                conc[l,idx] = melts.engine.getProperty('dispComposition', phase, i)

            if endmembers is None:
                endmembers = melts.endMemberFormulas[phase]
                thermo_properties = []
                thermo = ['activity', 'activity0', 'mu','mu0', 'X']
                thermodynamics = np.empty((L,len(endmembers)*len(thermo)))

            idx = 0
            for i in thermo:
                for j in endmembers:
                    if l == 0:
                        thermo_properties.append(i+'_'+j)
                    thermodynamics[l,idx] = melts.engine.getProperty(i, phase, j)
                    idx = idx + 1

        Output = {'Properties': pd.DataFrame(columns = properties, data = results), 
                'Composition': pd.DataFrame(columns = oxides, data = conc),
                'Thermodynamics': pd.DataFrame(columns = thermo_properties, data = thermodynamics)}
    
    return Output

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
    Results = {phases[0]: np.nan, phases[1]: np.nan, phases[2]: np.nan, 'T_Liq_C': np.nan, 'H2O_melt': np.nan}
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

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    try:
        Liq_Results = findLiq_MELTS(P_bar = P_bar, comp = bulk, T_C_init = T_initial_C, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        Results['T_Liq_C'] = Liq_Results['T_Liq_C']
        Results['H2O_melt'] = Liq_Results['H2O']
    except:
        return Results

    if type(H2O_Liq) == np.ndarray:
        if Results['H2O_melt'] < 0.99*bulk[14]:
            return Results


    T = Results['T_Liq_C']
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

def path_MELTS(Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_initial_C = None,
               T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, T_maxdrop_C = None,
               P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
               isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = False, 
               fO2_buffer = None, fO2_offset = None, fluid_sat = False, Crystallinity_limit = None, 
               Suppress = ['rutile', 'tridymite'], Suppress_except=False, phases=None, trail = None, melts = None,
               thermo_prop = False, shared_dict = None, index = None, q = None):
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

    # if shared_dict is not None:
    #     shared_dict[index] = copy.deepcopy(Results)

    if trail is not None:
        trail = False

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

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    melts.engine.setSystemProperties("Output",  "none")

    if Suppress_except is False:
        if Suppress is not None:
            if Suppress == "All":
                melts.engine.pressure = 500 #np.random.normal(500, 500/10)
                melts.engine.temperature = 1200 + 200
                melts.engine.setBulkComposition(bulk)
                PL = melts.engine.calcSaturationState()
                for p in PL:
                    if p != "fluid":
                        if p != "water":
                            melts.engine.setSystemProperties("Suppress", p)
            else:
                if type(Suppress) == list:
                    for p in Suppress:
                        melts.engine.setSystemProperties("Suppress", p)
                else:
                    melts.engine.setSystemProperties("Suppress", Suppress)
    else:
        melts.engine.pressure = 500 #np.random.normal(500, 500/10)
        melts.engine.temperature = 1200 + 200
        melts.engine.setBulkComposition(bulk)
        PL = melts.engine.calcSaturationState()
        for p in PL:
            if p != "fluid":
                if p != "water":
                    if type(Suppress_except) == list:
                        if p not in Suppress_except:
                            melts.engine.setSystemProperties("Suppress", p)
                    else:
                        if p != Suppress_except:
                            melts.engine.setSystemProperties("Suppress", p)

    if find_liquidus:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    T_initial_C = estimate_t(comp = comp, P_bar = P_path_bar[0]) + 60.0
                    Liq_Results = findLiq_MELTS(P_bar = P_path_bar[0], comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_initial_C)
                else:
                    T_initial_C = estimate_t(comp = comp, P_bar = P_path_bar) + 60.0
                    Liq_Results = findLiq_MELTS(P_bar = P_path_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_initial_C)
            except:
                print('Liquidus calculation failed.')
                if trail is not None:
                    return Results, trail
                else:
                    return Results
        elif P_start_bar is not None:
            try:
                T_initial_C = estimate_t(comp = comp, P_bar = P_start_bar) + 60.0
                Liq_Results = findLiq_MELTS(P_bar = P_start_bar, comp = bulk, melts = melts, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_initial_C)
            except:
                print('Liquidus calculation failed')
                if trail is not None:
                    return Results, trail
                else:
                    return Results

        T_start_C = Liq_Results['T_Liq_C'] + 0.1
        if T_end_C is None and T_maxdrop_C is not None:
            T_end_C = T_start_C - T_maxdrop_C

    else:
        if fO2_buffer is not None:
            if fO2_offset is None:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer])
            else:
                melts.engine.setSystemProperties(["Log fO2 Path: " + fO2_buffer, "Log fO2 Offset: " + str(fO2_offset)])

    T = None
    P = None

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

    if type(P) == np.ndarray and T_end_C is not None and dt_C is None:
        T = np.linspace(T_start_C, T_end_C, len(P))

    if type(T) == np.ndarray and P_end_bar is None and dp_bar is not None:
        P = np.linspace(P_start_bar, P_start_bar - dp_bar*(len(T)-1), len(T))
    if type(T) == np.ndarray and P_end_bar is not None and dp_bar is None:
        P = np.linspace(P_start_bar, P_end_bar, len(T))
    elif type(P) == np.ndarray and T_end_C is None and dt_C is not None:
        T = np.linspace(T_start_C, T_start_C - dt_C*(len(P)-1), len(P))
   

    if type(T) == np.ndarray and type(P) == np.ndarray:
        if len(T) != len(P):
            raise Exception("Length of P and T vectors are not the same. Check input parameters")

    if find_liquidus is False:
        if type(T) != np.ndarray:
            melts.engine.temperature = T
        else:
            melts.engine.temperature = T[0]

        if type(P) != np.ndarray:
            melts.engine.pressure = P
        else:
            melts.engine.pressure = P[0]

    else:
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

    if fluid_sat:
        melts.engine.setSystemProperties("Mode", "Fractionate Fluids")
        melts.engine.calcEquilibriumState(1,1)

    if type(T) == np.ndarray:
        length = len(T)
    else:
        length = len(P)

    properties = ['g', 'h', 's', 'v', 'cp', 'dcpdt', 'dvdt', 'dpdt', 'd2vdt2', 'd2vdtdp', 'd2vdp2', 'molwt', 'rho', 'mass']

    phase_prop = {}
    for i in range(length):
        if i == 1:
            if trail is not None:
                trail = True

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
                if trail is not None:
                    trail = False
                # if trail is not None:
                    return Results, trail
                else:
                    return Results
                # break

        if isenthalpic is not None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(2,1)
                else:
                    melts.engine.calcEquilibriumState(2,0)
            except:
                if trail is not None:
                    trail = False
                # if trail is not None:
                #     return Results, trail
                # else:
                #     return Results
                break

        if isentropic is not None:
            try:
                if Frac_solid is not None or Frac_fluid is not None:
                    melts.engine.calcEquilibriumState(3,1)
                else:
                    melts.engine.calcEquilibriumState(3,0)
            except:
                if trail is not None:
                    trail = False
                # return Results
                break

        if isochoric is not None:
            try:
                melts.engine.calcEquilibriumState(4,0)
            except:
                if trail is not None:
                    trail = False
                # return Results
                break
            
        try:
            step_results = {}
                
            if "Quadratic iterations" in melts.engine.status.message:
                if trail is not None:
                    trail = False
                    return Results, trail
                else:
                    return Results
            if "Convergence error" in melts.engine.status.message:
                if trail is not None:
                    trail = False
                    return Results, trail
                else:
                    return Results

            step_results['Conditions'] = np.array([melts.engine.temperature, melts.engine.pressure, melts.engine.getProperty('h', 'bulk'),
                                              melts.engine.getProperty('s', 'bulk'),melts.engine.getProperty('v', 'bulk'),melts.engine.getProperty('mass', 'bulk'),
                                              melts.engine.getProperty('dvdp', 'bulk'),melts.engine.getProperty('logfO2')])

            if melts.engine.getProperty('mass', 'liquid1') > 0.0:
                PhaseList = ['liquid1'] + melts.engine.solidNames
            else:
                PhaseList = melts.engine.solidNames

            for phase in PhaseList:
                compo = np.array(melts.engine.getProperty('dispComposition', phase))
                step_results[phase] = compo[[0,1,2,3,4,5,6,7,10,11,12,13,14,15]]
                if phase not in list(phase_prop.keys()):
                    if thermo_prop:
                        thermo_properties = []
                        endmembers = melts.endMemberFormulas[phase[:-1]]
                        thermo = ['activity', 'activity0', 'mu','mu0', 'X']
                        for iii in thermo:
                            for jjj in endmembers:
                                thermo_properties.append(iii+'_'+jjj)
                            # melts.engine.calcEndMemberProperties(phase, melts.engine.getProperty('dispComposition', phase))
                        full_prop = properties + thermo_properties
                    else:
                        full_prop = properties

                    phase_prop[phase] = full_prop

                step_results[phase+'_prop_keys'] = phase_prop[phase]
                prop_phase = []
                for prop_name in phase_prop[phase]:
                    if prop_name in properties:
                        prop_phase.append(melts.engine.getProperty(prop_name, phase))
                    else:
                        if thermo_prop:
                            melts.engine.calcEndMemberProperties(phase, melts.engine.getProperty('dispComposition', phase))
                            prop_phase.append(melts.engine.getProperty(prop_name.split('_')[0], phase, prop_name.split('_')[1]))
                
                step_results[phase+'_prop'] = np.array(prop_phase)

            if q is not None:
                q.put([index, i, step_results])

            Results[i] = step_results
        except:
            if trail is not None:
                trail = False
            break

        

        if Crystallinity_limit is not None:
            Volume = 0
            Fluid_List = ['liquid1', 'water1', 'fluid1']
            for phase in PhaseList:
                if phase not in Fluid_List:
                    Volume = Volume + float(Results[phase + '_prop']['v'].loc[i])

            Total_volume = Volume + float(Results['liquid1_prop']['v'].loc[i])

            if Volume/Total_volume > Crystallinity_limit:
                break

        if phases is not None:
            ll = 0
            for p in phases:
                if p in Results.keys():
                    ll = ll + 1
            
            if ll == len(phases):
                break

        if i != length - 1:
            try:
                melts = melts.addNodeAfter()
            except:
                if trail is not None:
                    return Results, False
                else:
                    return Results
    
    if trail is not None:
        return Results, trail
    else:
        return Results

def findSatPressure_MELTS_multi(Model = None, comp = None, fO2_buffer = None, fO2_offset = None, P_bar = None, T_fixed_C = None, index = None):
    out = pd.DataFrame(columns = ['Sample_ID_Liq', 'SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq', 'P_bar', 'T_Liq_C'])

    from meltsdynamic import MELTSdynamic

    if Model is None or Model == "MELTSv1.0.2":
        melts = MELTSdynamic(1)
    elif Model == "pMELTS":
        melts = MELTSdynamic(2)
    elif Model == "MELTSv1.1.0":
        melts = MELTSdynamic(3)
    elif Model == "MELTSv1.2.0":
        melts = MELTSdynamic(4)

    # melts.engine.setSystemProperties("Suppress", "rutile")
    # melts.engine.setSystemProperties("Suppress", "tridymite")

    if "Sample_ID_Liq" not in comp.keys():
        if index is not None:
            comp.loc[:, "Sample_ID_Liq"] = index
        else:
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
            melts.engine.pressure = P_bar[i] # np.round(np.random.normal(P_bar[i], P_bar[i]/10))
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
            'T_Liq_C': T_fixed_C[i],
            'Sample_ID_Liq': comp.loc[i, "Sample_ID_Liq"]}
        
        out.loc[i,:] = res

    return out

def findSatPressure_MELTS(Model = None, T_C_init = None, T_fixed_C = None, P_bar_init = None, 
                          comp = None, fO2_buffer = None, fO2_offset = None, trial = None, melts = None, suppressed = None):
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
            'T_Liq_C': the temperature in degrees Celsius at which the saturation pressure is calculated.
    """
    if trial is not None:
        trial = False
    
    if P_bar_init is None:
        P_bar_init = 5000

    P_bar = P_bar_init # np.random.normal(P_bar_init, P_bar_init/20)

    if T_C_init is None:
        T_C_init = 1200

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

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    if T_fixed_C is not None:
        if suppressed is None:
            melts.engine.pressure = P_bar
            melts.engine.temperature = T_fixed_C + 500
            melts.engine.setBulkComposition(bulk)
            PL = melts.engine.calcSaturationState()
            for p in PL:
                if p != "fluid":
                    if p != "water":
                        melts.engine.setSystemProperties("Suppress", p)

        melts = melts.addNodeAfter()
        melts.engine.pressure = P_bar
        melts.engine.temperature = T_fixed_C
        melts.engine.setBulkComposition(bulk)
        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            if trial is not None:
                out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                return out, trial
            else:
                out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
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
                        if trial is not None:
                            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                            return out, trial
                        else:
                            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
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
                        if trial is not None:
                            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                            return out, trial
                        else:
                            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
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
                if trial is not None:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                    return out, trial
                else:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
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
            'T_Liq_C': T_fixed_C}

    else:    
        try:
            T_C_init = estimate_t(comp = comp, P_bar = P_bar) + 60.0
            Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = T_C_init, Step = np.array([5,1]))
        except:
            if trial is not None:
                out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                return out, trial
            else:
                return {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
        
        if Liq_Results['T_Liq_C'] == 0.0:
            out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
            return out
        
        if Liq_Results['fluid_saturated'] == "Yes":
            while Liq_Results['fluid_saturated'] == "Yes":
                P_low = P_bar
                if P_bar > 5000:
                    P_bar = P_bar*1.2
                else:
                    P_bar = P_bar*1.6

                try:
                    Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq_C'], Step = np.array([5,1]))
                except:
                    if trial is not None:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out, trial
                    else:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out
                
                if Liq_Results['T_Liq_C'] == 0.0:
                    if trial is not None:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out, trial
                    else:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out
                
                if Liq_Results['fluid_saturated'] == "No":
                    P_high = P_bar
                    break
        
        else:
            while Liq_Results['fluid_saturated'] == "No":
                P_high = P_bar
                P_bar = P_bar/2
                try:
                    Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq_C'], Step = np.array([5,1]))
                except:
                    if trial is not None:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out, trial
                    else:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out
                
                if Liq_Results['T_Liq_C'] == 0.0:
                    if trial is not None:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out, trial
                    else:
                        out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                        return out
                
                if Liq_Results['fluid_saturated'] == "Yes":
                    P_low = P_bar
                    break

        for jj in range(5):
            P_new = P_low + (P_high - P_low)/2

            try:
                Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_new, comp = bulk, melts = melts, 
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, T_C_init = Liq_Results['T_Liq_C'], Step = np.array([3,1]))
            except:
                if trial is not None:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                    return out, trial
                else:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                    return out
            
            if Liq_Results['T_Liq_C'] == 0.0:
                if trial is not None:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                    return out, trial
                else:
                    out = {'SiO2_Liq': 0.0, 'TiO2_Liq': 0.0, 'Al2O3_Liq': 0.0, 'FeOt_Liq': 0.0, 'MnO_Liq': 0.0, 'MgO_Liq': 0.0, 'CaO_Liq': 0.0, 'Na2O_Liq': 0.0, 'K2O_Liq': 0.0, 'P2O5_Liq': 0.0, 'H2O_Liq': 0.0, 'CO2_Liq': 0.0, 'Fe3Fet_Liq': 0.0, 'P_bar': 0.0, 'T_Liq_C': 0.0}
                    return out

            if Liq_Results['fluid_saturated'] == "Yes":
                P_low = P_new
            else:
                P_high = P_new

            if P_high/P_low < 1.01:
                break

        P_bar = round((P_high+P_low)/2)
        try:
            Liq_Results = findLiq_MELTS(Model = Model, P_bar = P_bar, comp = bulk, melts = melts, 
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, 
                                            T_C_init = Liq_Results['T_Liq_C'], Step = np.array([1]))
        except:
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
            'T_Liq_C': Liq_Results['T_Liq_C']}
              
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
            'T_Liq_C': Liq_Results['T_Liq_C']}

    if trial is not None:
        trial = True
        return out, trial
    else:
        return out

def AdiabaticDecompressionMelting_MELTS(Model = None, comp = None, Tp_C = None, Tp_Method = None, T_start_C = None,
                                        P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
                                        Frac = False, fO2_buffer = None, fO2_offset = None):
    if T_start_C is None:
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
    
    
    if T_start_C is None:
        lz = m.lithologies.matthews.klb1()
        mantle = m.mantle([lz], [1], ['Lz'])
        if P_start_bar is not None:
            T_start_C = mantle.adiabat(P_start_bar/10000, Tp_C)
        else:
            T_start_C = mantle.adiabat(P_path_bar[0]/10000, Tp_C)

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    melts.engine.setBulkComposition(bulk)
    melts.engine.temperature = T_start_C

    if P_path_bar is not None:
        P = P_path_bar
    else:
        P = np.linspace(P_start_bar, P_end_bar, 1+int(round((P_start_bar - P_end_bar)/dp_bar)))

    # Results['Conditions'] = pd.DataFrame(data = np.zeros((len(P), 7)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'dvdp', 'logfO2'])
    Results['Conditions'] = pd.DataFrame(data = np.zeros((len(P), 8)), columns = ['temperature', 'pressure', 'h', 's', 'v', 'mass', 'dvdp', 'logfO2'])
    properties = ['g', 'h', 's', 'v', 'cp', 'dcpdt', 'dvdt', 'dpdt', 'd2vdt2', 'd2vdtdp', 'd2vdp2', 'molwt', 'rho', 'mass']
    
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
                Results['Conditions'].loc[k,R] = melts.engine.temperature
            elif R == 'pressure':
                Results['Conditions'].loc[k,R] = melts.engine.pressure
            elif R == 'logfO2':
                try:
                    Results['Conditions'].loc[k,R] = melts.engine.getProperty(R)
                except:
                    Results['Conditions'].loc[k,R] = np.nan
            else:
                Results['Conditions'].loc[k,R] = melts.engine.getProperty(R, 'bulk')

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
                # Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((len(P), 4)), columns = ['h', 'mass', 'v', 'rho'])
                Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((len(P), len(properties))), columns = properties)
                                                        
            if phase in list(Results.keys()):
                for el in Results[phase]:
                    Results[phase].loc[k,el] = melts.engine.getProperty('dispComposition', phase, el)

                for pr in Results[phase + '_prop']:
                    Results[phase + '_prop'].loc[k,pr] = melts.engine.getProperty(pr, phase)
        
        # if k == 0:
        #     s = melts.engine.getProperty('s','bulk')

    return Results


# ─────────────────────────────────────────────────────────────────────────────
# Worker-process globals (one MELTS instance per worker, initialised once)
# ─────────────────────────────────────────────────────────────────────────────
_PARALLEL_MELTS = None

# Property fields returned by calcPhaseProperties in this order (from MELTSstatus.fields)
_MELTS_FIELDS = (
    'g', 'h', 's', 'v', 'cp', 'dcpdt',
    'dvdt', 'dvdp', 'd2vdt2', 'd2vdtdp', 'd2vdp2',
    'molwt', 'rho', 'mass',
)

_MODEL_CODES = {
    "MELTSv1.0.2":  1,
    "pMELTS":       2,
    "MELTSv1.1.0":  3,
    "MELTSv1.2.0":  4,
}


def _init_melts_worker(model_code: int, extra_paths: list):
    """ProcessPoolExecutor initializer – run once per worker process."""
    global _PARALLEL_MELTS
    import warnings
    for p in extra_paths:
        if p not in sys.path:
            sys.path.insert(0, p)
    from meltsdynamic import MELTSdynamic
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _PARALLEL_MELTS = MELTSdynamic(model_code)


def _worker_phase_props(phase_name: str,
                        T_arr,     # ndarray (n,)
                        P_arr,     # ndarray (n,)
                        bulk_mat,  # ndarray (n, 19)
                        ):
    """Evaluate MELTS supplemental calculator for a batch of rows.

    Returns ndarray shape (n, 14) in _MELTS_FIELDS order. NaN on failure.
    """
    global _PARALLEL_MELTS
    melts = _PARALLEL_MELTS
    n = len(T_arr)
    out = np.full((n, 14), np.nan, dtype=np.float64)

    for i in range(n):
        try:
            melts.engine.temperature = float(T_arr[i])
            melts.engine.pressure = float(P_arr[i])
            melts.engine.calcPhaseProperties(phase_name, bulk_mat[i].tolist())
            if not melts.engine.status.failed:
                for j, field in enumerate(_MELTS_FIELDS):
                    d = melts.engine.__dict__.get(field)
                    if isinstance(d, dict):
                        val = d.get(phase_name)
                        if val is not None:
                            try:
                                out[i, j] = float(val)
                            except (TypeError, ValueError):
                                pass
        except Exception:
            pass
    return out


def _build_bulk_matrix(comp_liq: pd.DataFrame) -> np.ndarray:
    """Convert a comp_liq DataFrame to an (n, 19) MELTS-order oxide mass matrix.

    MELTS oxide order for MELTSv1.0.2:
    0:SiO2  1:TiO2  2:Al2O3  3:Fe2O3  4:Cr2O3  5:FeO  6:MnO  7:MgO
    8:NiO(0)  9:CoO(0)  10:CaO  11:Na2O  12:K2O  13:P2O5  14:H2O  15:CO2
    16:SO3(0)  17:Cl2O-1(0)  18:F2O-1(0)
    """
    n = len(comp_liq)
    mat = np.zeros((n, 19), dtype=np.float64)
    _LIQ_TO_IDX = {
        'SiO2_Liq': 0, 'TiO2_Liq': 1, 'Al2O3_Liq': 2, 'Cr2O3_Liq': 4,
        'MnO_Liq':  6, 'MgO_Liq':  7, 'CaO_Liq':   10, 'Na2O_Liq': 11,
        'K2O_Liq': 12, 'P2O5_Liq': 13, 'H2O_Liq':   14, 'CO2_Liq':  15,
    }
    for col, idx in _LIQ_TO_IDX.items():
        if col in comp_liq.columns:
            mat[:, idx] = comp_liq[col].values

    fe3fet = comp_liq['Fe3Fet_Liq'].values if 'Fe3Fet_Liq' in comp_liq.columns else np.zeros(n)
    feot   = comp_liq['FeOt_Liq'].values   if 'FeOt_Liq'   in comp_liq.columns else np.zeros(n)
    mat[:, 3] = fe3fet * (159.59 / 2.0 / 71.844) * feot   # Fe2O3
    mat[:, 5] = (1.0 - fe3fet) * feot                       # FeO
    return mat

def calc_phase_props_MELTS(Results, Model="MELTSv1.0.2", fO2_buffer=None, fO2_offset=None, melts=None):
    """
    Calculates per-phase thermodynamic properties for all phases in a PTT/nGibbs
    output dictionary using supCalc_MELTS, and computes bulk thermodynamic properties.

    Parameters
    ----------
    Results : dict
        PTT/nGibbs-format dict containing at minimum:
        - 'Conditions': DataFrame with 'T_C' and 'P_bar' columns
        - '<phase>': DataFrames with phase oxide compositions
          (columns like 'SiO2_<suffix>', 'FeOt_<suffix>', 'Fe3Fet_<suffix>')
        - 'mass_g': DataFrame with phase masses (g per 100 g total system mass)
    Model : str
        MELTS model: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS".
    fO2_buffer : str, optional
    fO2_offset : float, optional
    melts : MELTSdynamic, optional
        Pre-initialised MELTS instance; one will be created if not provided.

    Returns
    -------
    Results : dict
        Input dict with the following additions/modifications:
        - '<phase>_prop' DataFrames with columns:
            mass_g, rho_kg/m^3, V_cm^3, G_J, H_J, S_J/K, Cp_J/(kg.K^2),
            dCpdT_J/(kg.K^2), dVdT_cm^3/K, dPdT_bar/K, d2VdT2_cm^3/K^2,
            d2VdTdP_cm^3/(bar.K), d2VdP2_cm^3/bar^2, molwt
          NaN where the phase has zero mass. Extensive properties (V, G, H, S,
          Cp, etc.) are scaled to actual phase mass (nGibbs_mass / 100 * supCalc
          value) so that Conditions bulk sums are straightforward. Intensive
          properties (rho, molwt, dPdT) are stored unscaled.
        - Updated 'Conditions' DataFrame with columns:
            T_C, P_bar, mass_g, H_J, S_J/K, V_cm^3, rho_kg/m^3, log10(fO2),
            dVdP_cm^3/bar (NaN: dvdp is not available per-phase from supCalc_MELTS)
    """
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

    T_C = Results['Conditions']['T_C'].values
    P_bar_vals = Results['Conditions']['P_bar'].values
    n = len(T_C)

    _SKIP_KEYS = {'Conditions', 'mass_g', 'All'}
    phases = [k for k in Results if k not in _SKIP_KEYS and not k.endswith('_prop')]

    prop_col_names = [
        'mass_g', 'rho_kg/m^3', 'V_cm^3', 'G_J', 'H_J', 'S_J/K',
        'Cp_J/(kg.K^2)', 'dCpdT_J/(kg.K^2)', 'dVdT_cm^3/K', 'dPdT_bar/K',
        'd2VdT2_cm^3/K^2', 'd2VdTdP_cm^3/(bar.K)', 'd2VdP2_cm^3/bar^2', 'molwt',
    ]
    melts_prop_names = [
        'mass', 'rho', 'v', 'g', 'h', 's',
        'cp', 'dcpdt', 'dvdt', 'dpdt',
        'd2vdt2', 'd2vdtdp', 'd2vdp2', 'molwt',
    ]
    # Intensive properties are stored as-is; all others are scaled by phase mass fraction
    _INTENSIVE = {'rho', 'molwt', 'dpdt'}

    _BASE_OXIDES = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'MnO', 'MgO',
                    'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']

    has_mass_table = 'mass_g' in Results

    for phase in phases:
        # reset_index so active_idx (0-based positions from np.where) matches labels.
        # divide_ptt_tables preserves original row labels, so sub-tables may start at
        # non-zero values (e.g. 16, 32, ...) causing .loc[active_idx] to fail.
        phase_df = Results[phase].reset_index(drop=True)

        sio2_cols = [c for c in phase_df.columns if c.startswith('SiO2_')]
        if not sio2_cols:
            continue
        suffix = sio2_cols[0].split('SiO2_', 1)[1]

        if has_mass_table and phase in Results['mass_g'].columns:
            active = Results['mass_g'][phase].values > 0
        else:
            active = np.ones(n, dtype=bool)
        active_idx = np.where(active)[0]

        prop_df = pd.DataFrame(np.nan, index=range(n), columns=prop_col_names)

        if len(active_idx) == 0:
            Results[phase + '_prop'] = prop_df
            continue

        # Build comp DataFrame in _Liq format expected by supCalc_MELTS
        na = len(active_idx)
        comp_liq = pd.DataFrame(index=range(na), dtype=float)

        for ox in _BASE_OXIDES:
            src = f'{ox}_{suffix}'
            comp_liq[f'{ox}_Liq'] = (
                phase_df.loc[active_idx, src].fillna(0).values
                if src in phase_df.columns else np.zeros(na)
            )

        feot_src = f'FeOt_{suffix}'
        fe3_src = f'Fe3Fet_{suffix}'

        if feot_src in phase_df.columns:
            comp_liq['FeOt_Liq'] = phase_df.loc[active_idx, feot_src].fillna(0).values
        else:
            feo = (phase_df.loc[active_idx, f'FeO_{suffix}'].fillna(0).values
                   if f'FeO_{suffix}' in phase_df.columns else np.zeros(na))
            fe2o3 = (phase_df.loc[active_idx, f'Fe2O3_{suffix}'].fillna(0).values
                     if f'Fe2O3_{suffix}' in phase_df.columns else np.zeros(na))
            comp_liq['FeOt_Liq'] = feo + (71.844 / (159.69 / 2)) * fe2o3

        if fe3_src in phase_df.columns:
            comp_liq['Fe3Fet_Liq'] = phase_df.loc[active_idx, fe3_src].fillna(0).values
        else:
            fe2o3 = (phase_df.loc[active_idx, f'Fe2O3_{suffix}'].fillna(0).values
                     if f'Fe2O3_{suffix}' in phase_df.columns else np.zeros(na))
            feot = comp_liq['FeOt_Liq'].values
            comp_liq['Fe3Fet_Liq'] = np.where(feot > 0, (71.844 / (159.69 / 2)) * fe2o3 / feot, 0.0)

        if has_mass_table and phase in Results['mass_g'].columns:
            phase_mass = Results['mass_g'][phase].values[active_idx]
        else:
            phase_mass = np.full(na, 100.0)
        scale = phase_mass / 100.0

        bulk_mat = _build_bulk_matrix(comp_liq)

        for ii, row_idx in enumerate(active_idx):
            try:
                melts.engine.temperature = float(T_C[row_idx])
                melts.engine.pressure = float(P_bar_vals[row_idx])
                melts.engine.calcPhaseProperties(phase, bulk_mat[ii].tolist())
                if melts.engine.status.failed:
                    continue
                prop_df.loc[row_idx, 'mass_g'] = phase_mass[ii]
                for melts_field, user_col in zip(melts_prop_names, prop_col_names):
                    if user_col == 'mass_g':
                        continue
                    d = melts.engine.__dict__.get(melts_field)
                    val = d.get(phase) if isinstance(d, dict) else None
                    if val is not None:
                        try:
                            v = float(val)
                            prop_df.loc[row_idx, user_col] = (
                                v if melts_field in _INTENSIVE else v * scale[ii]
                            )
                        except (TypeError, ValueError):
                            pass
            except Exception:
                pass

        Results[phase + '_prop'] = prop_df

    # ---- Bulk Conditions (simple sums of already-scaled extensive properties) ----
    total_mass = np.zeros(n)
    total_H = np.zeros(n)
    total_S = np.zeros(n)
    total_V = np.zeros(n)

    for phase in phases:
        prop_key = phase + '_prop'
        if prop_key not in Results:
            continue
        prop_df = Results[prop_key]
        total_mass += np.nan_to_num(prop_df['mass_g'].values)
        total_H += np.nan_to_num(prop_df['H_J'].values)
        total_S += np.nan_to_num(prop_df['S_J/K'].values)
        total_V += np.nan_to_num(prop_df['V_cm^3'].values)

    with np.errstate(divide='ignore', invalid='ignore'):
        rho_bulk = np.where(total_V > 0, (total_mass / total_V) * 1000.0, np.nan)

    existing_cond = Results['Conditions']
    if 'log10(fO2)' in existing_cond.columns:
        fO2_vals = existing_cond['log10(fO2)'].values
    elif 'logfO2' in existing_cond.columns:
        fO2_vals = existing_cond['logfO2'].values
    else:
        fO2_vals = np.full(n, np.nan)

    Results['Conditions'] = pd.DataFrame({
        'T_C': T_C,
        'P_bar': P_bar_vals,
        'mass_g': total_mass,
        'H_J': total_H,
        'S_J/K': total_S,
        'V_cm^3': total_V,
        'rho_kg/m^3': rho_bulk,
        'log10(fO2)': fO2_vals,
        'dVdP_cm^3/bar': np.full(n, np.nan),
    })

    return Results


def calc_phase_props_MELTS_parallel(
    Results,
    Model="MELTSv1.0.2",
    fO2_buffer=None,
    fO2_offset=None,
    n_workers=None,
    _extra_sys_paths=None,
):
    """Parallel version of calc_phase_props_MELTS using ProcessPoolExecutor.

    Spawns n_workers subprocesses each with their own MELTS library instance,
    distributing the row-wise property evaluations across them.  Because the
    alphaMELTS C library uses process-global state, multiple threads inside
    one process cannot safely share an instance; separate processes are
    required.

    The supplemental-calculator calls (getMeltsPhaseProperties) are pure
    analytic EOS evaluations with no iterative minimisation, so they are
    embarrassingly parallel and scale linearly with worker count.

    Parameters
    ----------
    Results : dict
        nGibbs-style PTT dict.  Must contain 'Conditions' (with 'T_C' and
        'P_bar'), one DataFrame per phase, and optionally 'mass_g'.
    Model : str
        MELTS model string: "MELTSv1.0.2", "pMELTS", "MELTSv1.1.0", or
        "MELTSv1.2.0".
    fO2_buffer : str, optional
        fO2 buffer name (passed to setSystemProperties).  Currently not
        forwarded to workers; if supplied a ValueError is raised.
    n_workers : int or None
        Number of worker processes.  None → core_config.MAX_WORKERS.
        Pass 1 to disable parallelism (useful for profiling comparisons).
    _extra_sys_paths : list or None
        Extra sys.path entries to propagate to workers so they can import
        meltsdynamic.  Defaults to the current sys.path.

    Returns
    -------
    Results : dict
        Input dict augmented with '<phase>_prop' DataFrames and an updated
        'Conditions' DataFrame.  Column layout matches calc_phase_props_MELTS,
        with the following differences:
          - 'dVdP_cm^3/bar' is now populated (was NaN in serial version).
          - 'dPdT_bar/K' is derived as -dVdT / dVdP (= 0 when dVdP == 0).
    """
    if fO2_buffer is not None:
        raise ValueError(
            "fO2_buffer is not yet supported in calc_phase_props_MELTS_parallel. "
            "Use calc_phase_props_MELTS for fO2-buffered calculations."
        )

    import concurrent.futures as _cf
    import warnings

    if n_workers is None:
        from petthermotools.core_config import MAX_WORKERS
        n_workers = MAX_WORKERS

    if n_workers < 1:
        n_workers = 1

    model_code = _MODEL_CODES.get(Model, 1)

    if _extra_sys_paths is None:
        _extra_sys_paths = [p for p in sys.path if p]

    T_C_all   = Results['Conditions']['T_C'].values
    P_bar_all = Results['Conditions']['P_bar'].values
    n = len(T_C_all)

    _SKIP_KEYS = {'Conditions', 'mass_g', 'All'}
    phases = [k for k in Results if k not in _SKIP_KEYS and not k.endswith('_prop')]

    # Output column names (must match existing calc_phase_props_MELTS layout)
    prop_col_names = [
        'mass_g', 'rho_kg/m^3', 'V_cm^3', 'G_J', 'H_J', 'S_J/K',
        'Cp_J/(kg.K^2)', 'dCpdT_J/(kg.K^2)', 'dVdT_cm^3/K', 'dPdT_bar/K',
        'dVdP_cm^3/bar', 'd2VdT2_cm^3/K^2', 'd2VdTdP_cm^3/(bar.K)',
        'd2VdP2_cm^3/bar^2', 'molwt',
    ]
    # Map each output column to the index in _MELTS_FIELDS
    # _MELTS_FIELDS = g(0) h(1) s(2) v(3) cp(4) dcpdt(5) dvdt(6) dvdp(7)
    #                 d2vdt2(8) d2vdtdp(9) d2vdp2(10) molwt(11) rho(12) mass(13)
    _PROP_FIELD_IDX = {
        'mass_g':                13,   # mass
        'rho_kg/m^3':            12,   # rho  (intensive)
        'V_cm^3':                 3,   # v
        'G_J':                    0,   # g
        'H_J':                    1,   # h
        'S_J/K':                  2,   # s
        'Cp_J/(kg.K^2)':          4,   # cp
        'dCpdT_J/(kg.K^2)':       5,   # dcpdt
        'dVdT_cm^3/K':            6,   # dvdt
        # dPdT_bar/K is derived: -dvdt / dvdp
        'dVdP_cm^3/bar':          7,   # dvdp
        'd2VdT2_cm^3/K^2':        8,   # d2vdt2
        'd2VdTdP_cm^3/(bar.K)':   9,   # d2vdtdp
        'd2VdP2_cm^3/bar^2':     10,   # d2vdp2
        'molwt':                 11,   # molwt (intensive)
    }
    # Intensive properties are stored unscaled (no × phase_mass/100)
    _INTENSIVE = {'rho_kg/m^3', 'molwt'}

    _BASE_OXIDES = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'MnO', 'MgO',
                    'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2']
    has_mass_table = 'mass_g' in Results

    with _cf.ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_melts_worker,
        initargs=(model_code, _extra_sys_paths),
    ) as executor:

        for phase in phases:
            phase_df = Results[phase].reset_index(drop=True)
            sio2_cols = [c for c in phase_df.columns if c.startswith('SiO2_')]
            if not sio2_cols:
                continue
            suffix = sio2_cols[0].split('SiO2_', 1)[1]

            if has_mass_table and phase in Results['mass_g'].columns:
                active = Results['mass_g'][phase].values > 0
            else:
                active = np.ones(n, dtype=bool)
            active_idx = np.where(active)[0]

            prop_df = pd.DataFrame(np.nan, index=range(n), columns=prop_col_names)

            if len(active_idx) == 0:
                Results[phase + '_prop'] = prop_df
                continue

            na = len(active_idx)

            # ---- Build comp_liq DataFrame (same logic as calc_phase_props_MELTS) ----
            comp_liq = pd.DataFrame(index=range(na), dtype=float)
            for ox in _BASE_OXIDES:
                src = f'{ox}_{suffix}'
                comp_liq[f'{ox}_Liq'] = (
                    phase_df.loc[active_idx, src].fillna(0).values
                    if src in phase_df.columns else np.zeros(na)
                )
            feot_src = f'FeOt_{suffix}'
            fe3_src  = f'Fe3Fet_{suffix}'
            if feot_src in phase_df.columns:
                comp_liq['FeOt_Liq'] = phase_df.loc[active_idx, feot_src].fillna(0).values
            else:
                feo   = phase_df.loc[active_idx, f'FeO_{suffix}'].fillna(0).values   if f'FeO_{suffix}'   in phase_df.columns else np.zeros(na)
                fe2o3 = phase_df.loc[active_idx, f'Fe2O3_{suffix}'].fillna(0).values if f'Fe2O3_{suffix}' in phase_df.columns else np.zeros(na)
                comp_liq['FeOt_Liq'] = feo + (71.844 / (159.69 / 2)) * fe2o3
            if fe3_src in phase_df.columns:
                comp_liq['Fe3Fet_Liq'] = phase_df.loc[active_idx, fe3_src].fillna(0).values
            else:
                fe2o3    = phase_df.loc[active_idx, f'Fe2O3_{suffix}'].fillna(0).values if f'Fe2O3_{suffix}' in phase_df.columns else np.zeros(na)
                feot_v   = comp_liq['FeOt_Liq'].values
                comp_liq['Fe3Fet_Liq'] = np.where(feot_v > 0, (71.844 / (159.69 / 2)) * fe2o3 / feot_v, 0.0)

            if has_mass_table and phase in Results['mass_g'].columns:
                phase_mass = Results['mass_g'][phase].values[active_idx]
            else:
                phase_mass = np.full(na, 100.0)
            scale = phase_mass / 100.0

            # ---- Build MELTS oxide matrix (vectorised) ----
            bulk_mat  = _build_bulk_matrix(comp_liq)
            T_active  = T_C_all[active_idx]
            P_active  = P_bar_all[active_idx]

            # ---- Dispatch to workers ----
            chunk_size  = max(1, (na + n_workers - 1) // n_workers)
            chunk_starts = list(range(0, na, chunk_size))
            futures = []
            for start in chunk_starts:
                end = min(start + chunk_size, na)
                fut = executor.submit(
                    _worker_phase_props,
                    phase,
                    T_active[start:end],
                    P_active[start:end],
                    bulk_mat[start:end],
                )
                futures.append((start, end, fut))

            # ---- Collect results ----
            try:
                raw = np.full((na, 14), np.nan)
                for start, end, fut in futures:
                    raw[start:end] = fut.result()

                # dPdT (thermal pressure coefficient) derived from dvdt and dvdp
                dvdt_col = raw[:, 6]
                dvdp_col = raw[:, 7]
                with np.errstate(divide='ignore', invalid='ignore'):
                    dpdt_derived = np.where(dvdp_col != 0.0, -dvdt_col / dvdp_col, 0.0)

                prop_df.loc[active_idx, 'mass_g'] = phase_mass

                for col_name in prop_col_names:
                    if col_name == 'mass_g':
                        continue
                    elif col_name == 'dPdT_bar/K':
                        # Intensive derived property: -dVdT/dVdP  (stored unscaled)
                        prop_df.loc[active_idx, col_name] = dpdt_derived
                    elif col_name in _PROP_FIELD_IDX:
                        vals = raw[:, _PROP_FIELD_IDX[col_name]]
                        if col_name in _INTENSIVE:
                            prop_df.loc[active_idx, col_name] = vals
                        else:
                            prop_df.loc[active_idx, col_name] = vals * scale

            except Exception:
                pass

            Results[phase + '_prop'] = prop_df

    # ---- Bulk Conditions (sum of scaled extensive properties) ----
    total_mass = np.zeros(n)
    total_H    = np.zeros(n)
    total_S    = np.zeros(n)
    total_V    = np.zeros(n)
    total_dVdP = np.zeros(n)

    for phase in phases:
        prop_key = phase + '_prop'
        if prop_key not in Results:
            continue
        pdf = Results[prop_key]
        total_mass += np.nan_to_num(pdf['mass_g'].values)
        total_H    += np.nan_to_num(pdf['H_J'].values)
        total_S    += np.nan_to_num(pdf['S_J/K'].values)
        total_V    += np.nan_to_num(pdf['V_cm^3'].values)
        if 'dVdP_cm^3/bar' in pdf.columns:
            total_dVdP += np.nan_to_num(pdf['dVdP_cm^3/bar'].values)

    with np.errstate(divide='ignore', invalid='ignore'):
        rho_bulk = np.where(total_V > 0, (total_mass / total_V) * 1000.0, np.nan)

    existing_cond = Results['Conditions']
    if 'log10(fO2)' in existing_cond.columns:
        fO2_vals = existing_cond['log10(fO2)'].values
    elif 'logfO2' in existing_cond.columns:
        fO2_vals = existing_cond['logfO2'].values
    else:
        fO2_vals = np.full(n, np.nan)

    Results['Conditions'] = pd.DataFrame({
        'T_C':           T_C_all,
        'P_bar':         P_bar_all,
        'mass_g':        total_mass,
        'H_J':           total_H,
        'S_J/K':         total_S,
        'V_cm^3':        total_V,
        'rho_kg/m^3':    rho_bulk,
        'log10(fO2)':    fO2_vals,
        'dVdP_cm^3/bar': total_dVdP,
    })

    return Results


