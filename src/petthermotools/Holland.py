import numpy as np
import pandas as pd
import julia
import time
from julia import MAGEMinCalc

def findLiq_holland(P_bar = None, T_C_init = None, comp = None):

    if P_bar is None:
        raise Exception("Please specify a pressure for calculations")

    if T_C_init is None:
        T_C_init = 1300

    if comp is None:
        raise Exception("No composition specified")
    else:
        if type(comp) == list:
            bulk = comp
        else:
            bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]

    T_Liq = 0

    T = T_C_init

    bulk_in = bulk.copy()

    Liq = ["liq","fl"]

    start = time.time()

    PhaseList = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)

    i = set.intersection(set(Liq),set(PhaseList))

    Step = np.array([3,1,0.1])
    for k in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T - Step[k]
                    Ret = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T + Step[k]
                    Ret = MAGEMinCalc.satPhase(P_bar/1000, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

    if "liq" in Ret['Phase']:
        if Ret['Liq_Frac'] < 0.9:
            return T_Liq
        else:
            T_Liq = T
            return T_Liq
    else:
        return T_Liq

def crystallise_holland(comp = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, isochoric = None, find_liquidus = None):

    Results = {}

    if comp is None:
        raise Exception("No composition specified")

    if P_path_bar is None and P_start_bar is None:
        raise Exception("Initial P system must be defined")
    if T_path_C is None and T_start_C is None and find_liquidus is None:
        raise Exception("Starting temperature must be specified or the liquidus must be found")

    bulk = [comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['Fe3Fet']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq']), comp['Cr2O3_Liq'], comp['H2O_Liq']]
    bulk = list(100*np.array(bulk)/np.sum(bulk))

    # find liquidus if specified
    if find_liquidus is not None:
        if P_path_bar is not None:
            #try:
            if T_start_C is None and T_path_C is None:
                if type(P_path_bar) == np.ndarray:
                    T_Liq = findLiq_holland(P_bar = P_path_bar[0], comp = bulk)
                else:
                    T_Liq = findLiq_holland(P_bar = P_path_bar, comp = bulk)
            elif T_start_C is not None and T_path_C is None:
                if type(P_path_bar) == np.ndarray:
                    T_Liq = findLiq_holland(P_bar = P_path_bar[0], T_C_init = T_start_C, comp = bulk)
                else:
                    T_Liq = findLiq_holland(P_bar = P_path_bar, T_C_init = T_start_C, comp = bulk)
            else:
                if type(P_path_bar) == np.ndarray:
                    T_Liq = findLiq_holland(P_bar = P_path_bar[0], T_C_init = T_path_C[0], comp = bulk)
                else:
                    T_Liq = findLiq_holland(P_bar = P_path_bar, T_C_init = T_path_C[0], comp = bulk)
            #except:
            #    return Results
        elif P_start_bar is not None:
            #try:
            if T_start_C is None and T_path_C is None:
                T_Liq = findLiq_holland(P_bar = P_start_bar, comp = bulk)
            elif T_start_C is not None and T_path_C is None:
                T_Liq = findLiq_holland(P_bar = P_start_bar, T_C_init = T_start_C, comp = bulk)
            else:
                T_Liq = findLiq_holland(P_bar = P_start_bar, T_C_init = T_path_C[0], comp = bulk)
            #except:
            #    return Results

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

    if find_liquidus is not None:
        if type(P) == np.ndarray and type(T) == np.ndarray:
            T_Liq_loc = np.abs(T - T_Liq).argmin()
            if T[T_Liq_loc]>T_Liq:
                T = T[T_Liq_loc:]
                P = P[T_Liq_loc:]
            else:
                T = T[T_Liq_loc-1:]
                P = P[T_Liq_loc-1:]

    if type(T) == np.ndarray:
        length = len(T)
    else:
        length =len(P)

    Results['Conditions'] = pd.DataFrame(data = np.zeros((length, 2)), columns = ['temperature', 'pressure'])
    Results['liq'] = pd.DataFrame(data = np.zeros((length, 11)), columns = ['SiO2', 'Al2O3', 'CaO', 'MgO', 'FeO', 'K2O', 'Na2O', 'TiO2', 'O', 'Cr2O3', 'H2O'])

    if type(T) != np.ndarray:
        T_in = T
    if type(P) != np.ndarray:
        P_in = P/1000

    for i in range(length):
        if type(T) == np.ndarray:
            T_in = T[i]
        if type(P) == np.ndarray:
            P_in = P[i]/1000

        #try:
        Ret = MAGEMinCalc.satPhase(P_in, T_in, bulk)
        #except:
        #    return Results

        for R in Results['Conditions']:
            if R == 'temperature':
                Results['Conditions'][R].loc[i] = T_in
            elif R == 'pressure':
                Results['Conditions'][R].loc[i] = P_in*1000

        A = dict(zip(Ret['Oxides'], Ret['Liq_Comp']))
        for el in Results['liq']:
            Results['liq'][el].loc[i] = A[el]

        if Frac_solid is not None and Frac_fluid is not None:
            bulk = np.array(Res['Liq_Comp'])

    return Results