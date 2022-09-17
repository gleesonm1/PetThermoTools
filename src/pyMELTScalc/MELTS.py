import numpy as np
import pandas as pd
import sys
from pyMELTScalc.Barom import *
from pyMELTScalc.Holland import *
from pyMELTScalc.Liq import *

def findLiq_MELTS(P_bar = None, Model = None, T_C_init = None, comp = None, melts = None):

    if P_bar is None:
        raise Exception("Please specify a pressure for calculations")

    if Model is None:
        Model = "MELTS"

    if T_C_init is None:
        T_C_init = 1300

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp
        else:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1 - comp['Fe3Fet'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]

    T_Liq = 0
    H2O_Melt = 0

    from meltsdynamic import MELTSdynamic

    if melts is None:
        melts = MELTSdynamic(1)

    melts.engine.setBulkComposition(bulk)
    melts.engine.pressure = P_bar
    melts.engine.temperature = T_C_init

    Liq = ['liquid1','water1', 'fluid1']
    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        return T_Liq, H2O_Melt

    PhaseList = melts.engine.solidNames
    if PhaseList is None:
        PhaseList = ['liquid1']
    else:
        PhaseList = ['liquid1'] + PhaseList

    i = set.intersection(set(Liq),set(PhaseList))
    Step = np.array([3,1,0.1])

    for j in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                melts = melts.addNodeAfter()
                melts.engine.temperature = melts.engine.temperature - Step[j]
                try:
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    return T_Liq, H2O_Melt

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                melts = melts.addNodeAfter()
                melts.engine.temperature = melts.engine.temperature + Step[j]
                try:
                    melts.engine.calcEquilibriumState(1,0)
                except:
                    return T_Liq, H2O_Melt

                PhaseList = melts.engine.solidNames
                if PhaseList is None:
                    PhaseList = ['liquid1']
                else:
                    PhaseList = ['liquid1'] + PhaseList
                i = set.intersection(set(Liq),set(PhaseList))


    T_Liq = melts.engine.temperature
    H2O_Melt = melts.engine.getProperty('dispComposition', 'liquid1', 'H2O')

    return T_Liq, H2O_Melt

def satTemperature_MELTS(q, Model, phases, P, T_initial, bulk, dt, T_step):

    from meltsdynamic import MELTSdynamic

    melts = MELTSdynamic(1)

    # melts.engine.setBulkComposition(bulk)
    # melts.engine.pressure = P
    # melts.engine.temperature = T_initial

    a_sat = 0
    b_sat = 0
    if len(phases) == 3:
        c_sat = 0

    T_Liq = 0
    H2O_Melt = 0

    try:
        T_Liq, H2O_Melt = findLiq_MELTS(P, Model, T_initial, bulk, melts = melts)
    except:
        q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
        return

    # Liq = ['liquid1','water1', 'fluid1']
    # try:
    #     melts.engine.calcEquilibriumState(1,0)
    # except:
    #     if len(phases) == 2:
    #         q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
    #         return
    #     else:
    #         q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
    #         return
    #
    # PhaseList = melts.engine.solidNames
    # if PhaseList is None:
    #     PhaseList = ['liquid1']
    # else:
    #     PhaseList = ['liquid1'] + PhaseList
    #
    # i = set.intersection(set(Liq),set(PhaseList))
    # Step = np.array([3,1,0.1])
    #
    # for j in range(len(Step)):
    #     if len(i) == len(PhaseList):
    #         while len(i) == len(PhaseList):
    #             melts = melts.addNodeAfter()
    #             melts.engine.temperature = melts.engine.temperature - Step[j]
    #             try:
    #                 melts.engine.calcEquilibriumState(1,0)
    #             except:
    #                 if len(phases) == 2:
    #                     q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
    #                     return
    #                 else:
    #                     q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
    #                     return
    #
    #             PhaseList = melts.engine.solidNames
    #             if PhaseList is None:
    #                 PhaseList = ['liquid1']
    #             else:
    #                 PhaseList = ['liquid1'] + PhaseList
    #             i = set.intersection(set(Liq),set(PhaseList))
    #
    #     if len(i) < len(PhaseList):
    #         while len(i) < len(PhaseList):
    #             melts = melts.addNodeAfter()
    #             melts.engine.temperature = melts.engine.temperature + Step[j]
    #             try:
    #                 melts.engine.calcEquilibriumState(1,0)
    #             except:
    #                 if len(phases) == 2:
    #                     q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
    #                     return
    #                 else:
    #                     q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
    #                     return
    #
    #             PhaseList = melts.engine.solidNames
    #             if PhaseList is None:
    #                 PhaseList = ['liquid1']
    #             else:
    #                 PhaseList = ['liquid1'] + PhaseList
    #             i = set.intersection(set(Liq),set(PhaseList))
    #
    #
    # T_Liq = melts.engine.temperature
    # H2O_Melt = melts.engine.getProperty('dispComposition', 'liquid1', 'H2O')

    T_fin = T_Liq - dt

    melts = melts.addNodeAfter()
    melts.engine.temperature = T_Liq
    melts.engine.pressure = P

    try:
        melts.engine.calcEquilibriumState(1,0)
    except:
        if len(phases) == 2:
            q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
            return
        else:
            q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
            return

    j = 0
    while melts.engine.temperature>T_fin:
        if j == 0:
            j = 1
            melts = melts.addNodeAfter()
            melts.engine.temperature = T_Liq-0.1

        else:
            melts = melts.addNodeAfter()
            melts.engine.temperature = melts.engine.temperature - T_step

        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            if len(phases) == 2:
                q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
                return
            else:
                q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
                return

        PhaseList = melts.engine.solidNames

        if phases[0] in PhaseList and a_sat == 0:
            a_sat = melts.engine.temperature

        if phases[1] in PhaseList and b_sat == 0:
            b_sat = melts.engine.temperature

        if len(phases) == 3:
            if phases[2] == 'orthopyroxene1':
                if 'orthopyroxene1' in PhaseList or 'clinopyroxene2' in PhaseList and c_sat == 0:
                    c_sat = melts.engine.temperature
            else:
                if phases[2] in PhaseList and c_sat == 0:
                    c_sat = melts.engine.temperature

            if a_sat and b_sat and c_sat > 0:
                break

        if len(phases) == 2:
            if a_sat and b_sat > 0:
                break

    if len(phases) == 2:
        q.put([a_sat, b_sat, T_Liq, H2O_Melt, P])
        return
    else:
        q.put([a_sat, b_sat, c_sat, T_Liq, H2O_Melt, P])
        return

def crystallise_MELTS(Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, isochoric = None, find_liquidus = None):

    Results = {}

    if comp is None:
        raise Exception("No composition specified")

    if P_path_bar is None and P_start_bar is None:
        raise Exception("Initial P system must be defined")
    if T_path_C is None and T_start_C is None and find_liquidus is None:
        raise Exception("Starting temperature must be specified or the liquidus must be found")

    from meltsdynamic import MELTSdynamic

    if Model is None or Model == "MELTS":
        melts = MELTSdynamic(1)

    bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], comp['Fe3Fet']*((159.59/2)/71.844)*comp['FeOt_Liq'], 0.0, (1- comp['Fe3Fet'])*comp['FeOt_Liq'], comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    if find_liquidus is not None:
        if P_path_bar is not None:
            try:
                if type(P_path_bar) == np.ndarray:
                    T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_path_bar[0], comp = bulk, melts = melts)
                else:
                    T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_path_bar, comp = bulk, melts = melts)
            except:
                return Results
        elif P_start_bar is not None:
            try:
                T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_start_bar, comp = bulk, melts = melts)
            except:
                return Results

        T_start_C = T_Liq

    if T_path_C is None:
        if T_end_C is None and dt is None:
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

    if type(T) == np.ndarray and type(P) == np.ndarray:
        if len(T) != len(P):
            raise Exception("Length of P and T vectors are not the same. Check input parameters")

    if type(T) == np.ndarray and P_end_bar is None and dp_bar is not None:
        P = np.linspace(P_start_bar, P_start_bar - dp_bar*len(T), len(T))
    elif type(P) == np.ndarray and T_end_C is None and dt_C is not None:
        T = np.linspace(T_start_C, T_start_C - dt_C*len(P), len(P))

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

    if type(T) == np.ndarray:
        length = len(T)
    else:
        length =len(P)

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 5)), columns = ['temperature', 'pressure', 'h', 's', 'v'])
    Results['liquid1'] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
    Results['liquid1_prop'] = pd.DataFrame(data = np.zeros((length, 4)), columns = ['h', 'mass', 'v', 'rho'])

    for i in range(length):
        if type(T) == np.ndarray:
            melts.engine.temperature = T[i]
        if type(P) == np.ndarray:
            melts.engine.pressure = P[i]

        try:
            melts.engine.calcEquilibriumState(1,0)
        except:
            return Results

        if isochoric is not None:
            if i == 0:
                v = melts.engine.getProperty('v', 'bulk')
            else:
                v_new = melts.engine.getProperty('v', 'bulk')
                step = np.array([5, 2, 0.1])
                for j in range(step):
                    if v_new < v:
                        while v_new < v:
                            melts = melts.addNodeAfter()
                            melts.engine.pressure = melts.engine.pressure - step[j]
                            try:
                                melts.engine.calcEquilibriumState(1,0)
                            except:
                                return Results

                            v_new = melts.engine.getProperty('v', 'bulk')

                    if v_new > v:
                        while v_new > v:
                            melts = melts.addNodeAfter()
                            melts.engine.pressure = melts.engine.pressure + step[j]
                            try:
                                melts.engine.calcEquilibriumState(1,0)
                            except:
                                return Results

                            v_new = melts.engine.getProperty('v', 'bulk')

        for R in Results['Conditions']:
            if R == 'temperature':
                Results['Conditions'][R].loc[i] = melts.engine.temperature
            elif R == 'pressure':
                Results['Conditions'][R].loc[i] = melts.engine.pressure
            else:
                Results['Conditions'][R].loc[i] = melts.engine.getProperty(R, 'bulk')

        PhaseList = ['liquid1'] + melts.engine.solidNames

        for phase in PhaseList:
            if phase not in list(Results.keys()):
                Results[phase] = pd.DataFrame(data = np.zeros((length, 14)), columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
                Results[phase + '_prop'] = pd.DataFrame(data = np.zeros((length, 4)), columns = ['h', 'mass', 'v', 'rho'])

            if phase in list(Results.keys()):
                for el in Results[phase]:
                    Results[phase][el].loc[i] = melts.engine.getProperty('dispComposition', phase, el)

                for pr in Results[phase + '_prop']:
                    Results[phase + '_prop'][pr].loc[i] = melts.engine.getProperty(pr, phase)

        if Frac_solid is not None and Frac_fluid is not None:
            bulk = np.array(melts.engine.getProperty('dispComposition', 'liquid1'))
            bulk = list(melts.engine.getProperty('mass', 'liquid1')*bulk/np.sum(bulk))

        melts = melts.addNodeAfter()

        if Frac_solid is not None and Frac_fluid is not None:
            try:
                melts.engine.setBulkComposition(bulk)
            except:
                return Results

    return Results









