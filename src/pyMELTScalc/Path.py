import numpy as np
import pandas as pd
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.notebook import tqdm, trange

def multi_path(cores = None, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, Print_suppress = None, fluid_sat = None, Crystallinity_limit = None):
    '''
    Carry out multiple calculations in parallel. Allows isobaric, polybaric and isochoric crystallisation to be performed as well as isothermal, isenthalpic or isentropic decompression. All temperature inputs/outputs are reported in degrees celcius and pressure is reported in bars.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default will be determined using Multiprocessing.cpu_count().

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each  step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each  step. Default False.

    T_C: float or np.ndarray
        Calculation temperature - typically used when calculations are performed at a fixed T (e.g.,isothermal degassing).

    T_path_C: np.ndarray
        If a specified temperature path is to be used, T_path_C will be used to determine the T at each step of the model. If 2D, this indicates that multiple calculations with different T_path_C arrays are to be performed.

    T_start_C: float or np.ndarray
        Initial temperature used for path calculations.

    T_end_C: float or np.ndarray
        Final temperature in crystallisation calculations.

    dt_C: float or np.ndarray
        Temperature increment during crystallisation calculations.

    P_bar: float or np.ndarray
        Calculation pressure - typically used when calculations are performed at a fixed P (e.g.,isobaric crystallisation).

    P_path_bar: np.ndarray
        If a specified pressure path is to be used, P_path_bar will be used to determine the P at each step of the model. If 2D, this indicates that multiple calculations with different P_path_bar arrays are to be performed.

    P_start_bar: float or np.ndarray
        Initial pressure used for path calculations.

    P_end_bar: float or np.ndarray
        Final pressure in crystallisation calculations.

    dp_bar: float or np.ndarray
        Pressure increment during crystallisation calculations..

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, and type(Fe3Fet_Liq) == np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their Fe3Fet_Liq value. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable or an oxygen fugacity buffer must be chosen.

    H2O_Liq: float or np.ndarray
        H2O content of the initial melt phase. If type(comp) == dict, and type(H2O_Liq) = np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their H2O_Liq value. If comp is a pd.DataFrame, a single H2O_Liq value may be passes (float) and will be used as the initial melt H2O content for all starting compositions. Alternatively, if an array of H2O_Liq values is passed, equal to the number of compositions specified in comp, a different initial melt H2O value will be passed for each sample. If None, H2O_Liq must be specified in the comp variable.

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

    fO2_offset: float or np.ndarray
        Offset from the buffer spcified in fO2_buffer (log units).

    Print_suppress: True/False
        If True, print messages concerning the status of the thermodynamic calculations will not be displayed.

    Crystallinity_limit: float
        If value given, calculation will stop when the volume mass fraction of the system (excludin fluids) exceeds this value.

    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.

    '''

    # set default values if required
    if Model is None:
        Model == "MELTSv1.0.2"

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("Adding small amounts of H$_{2}$O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H$_{2}$O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    if cores is None:
        cores = multiprocessing.cpu_count()

    # specify the number of calculations to be performed in each sequence
    One = 0
    if type(comp) == pd.core.frame.DataFrame: # simplest scenario - one calculation per bulk composition imported
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq'])%cores
    else:
        if P_bar is not None and type(P_bar) == np.ndarray: #one calculation per P loaded
            A = len(P_bar)//cores
            B = len(P_bar)%cores
        elif T_C is not None and type(T_C) == np.ndarray: # one calculation per T loaded.
            A = len(T_C)//cores
            B = len(T_C)%cores
        elif P_start_bar is not None and type(P_start_bar) == np.ndarray: # one calculation per starting P
            A = len(P_start_bar)//cores
            B = len(P_start_bar)%cores
        elif T_start_C is not None and type(T_start_C) == np.ndarray: # one calculation per starting T
            A = len(T_start_C)//cores
            B = len(T_start_C)%cores
        elif P_path_bar is not None and len(np.shape(P_path_bar)) > 1: # one calculation per P path
            A = np.shape(P_path_bar)[0]//cores
            B = np.shape(P_path_bar)[0]%cores
        elif T_path_C is not None and len(np.shape(T_path_C)) > 1: # one calculation per T path
            A = np.shape(T_path_C)[0]//cores
            B = np.shape(T_path_C)[0]%cores
        elif fO2_offset is not None and type(fO2_offset) == np.ndarray: # one calculation per offset
            A = len(fO2_offset)//cores
            B = len(fO2_offset)%cores

        else: # just one calculation
            One = 1
            A = 1
            B = 0

    Group = np.zeros(A) + cores
    if B > 0:
        Group = np.append(Group, B)

    qs = []
    q = Queue()

    # perform calculation if only 1 calculation is specified
    if One == 1:
        if Print_suppress is None:
            print("Running " + Model + " calculation...", end = "", flush = True)
            s = time.time()

        p = Process(target = path, args = (q, 1),
                    kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                            'T_C': T_C, 'T_path_C': T_path_C, 'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                            'P_bar': P_bar, 'P_path_bar': P_path_bar, 'P_start_bar': P_start_bar, 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                            'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit})

        p.start()
        try:
            ret = q.get(timeout = 180)
        except:
            ret = []

        TIMEOUT = 5
        start = time.time()
        if p.is_alive():
            while time.time() <= TIMEOUT:
                if not p.is_alive():
                    p.join()
                    p.terminate()
                    break
                time.sleep(.1)
            else:
                p.terminate()
                p.join(5)
        else:
            p.join()
            p.terminate()

        if Print_suppress is None:
            print(" Complete (time taken = " + str(round(time.time() - s,2)) + " seconds)", end = "", flush = True)

        if len(ret) > 0:
            Results, index = ret
            Results = stich(Results, Model = Model)
            return Results
        else:
            Results = {}
            return Results

    else: # perform multiple crystallisation calculations
        # first make sure everything is the right length
        L = np.sum(Group)
        if P_bar is None:
            P_bar = [None] * int(L)
        elif type(P_bar) == float or type(P_bar) == int:
            P_bar = np.zeros(int(L)) + P_bar
        if T_C is None:
            T_C = [None] * int(L)
        elif type(T_C) == float or type(T_C) == int:
            T_C = np.zeros(int(L)) + T_C

        if P_start_bar is None:
            P_start_bar = [None] * int(L)
        elif type(P_start_bar) == float or type(P_start_bar) == int:
            P_start_bar = np.zeros(int(L)) + P_start_bar
        if T_start_C is None:
            T_start_C = [None] * int(L)
        elif type(T_start_C) == float or type(T_start_C) == int:
            T_start_C = np.zeros(int(L)) + T_start_C

        if P_path_bar is None:
            P_path_bar = [None] * int(L)
        elif len(np.shape(P_path_bar))  == 1:
            P_path_bar = np.vstack([P_path_bar] * int(L))
        if T_path_C is None:
            T_path_C = [None] * int(L)
        elif len(np.shape(T_path_C))  == 1:
            T_path_C = np.vstack([T_path_C] * int(L))

        if fO2_offset is None:
            fO2_offset = [None] * int(L)
        elif type(fO2_offset) == float or type(fO2_offset) == int:
            fO2_offset = np.zeros(int(L)) + fO2_offset

        for j in tqdm(range(len(Group))):
            ps = []

            if Print_suppress is None:
                print("Running " + Model + " calculations " + str(int(cores*j)) + " to " + str(int(cores*j) + Group[j] - 1) + " ...", end = "", flush = True)
                s = time.time()

            for i in range(int(cores*j), int(cores*j + Group[j])):
                if type(comp) == dict:
                    p = Process(target = path, args = (q, i),
                                kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit})
                else:
                    p = Process(target = path, args = (q, i),
                                kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit})

                ps.append(p)
                p.start()

            for p in ps:
                try:
                    ret = q.get(timeout = 180)
                except:
                    ret = []

                qs.append(ret)

            TIMEOUT = 5
            start = time.time()
            for p in ps:
                if p.is_alive():
                    while time.time() - start <= TIMEOUT:
                        if not p.is_alive():
                            p.join()
                            p.terminate()
                            break
                        time.sleep(.1)
                    else:
                        p.terminate()
                        p.join(5)
                else:
                    p.join()
                    p.terminate()

            if Print_suppress is None:
                print(" Complete (time taken = " + str(round(time.time() - s,2)) + " seconds)", end = "\n", flush = True)

        Results = {}
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Res, index = qs[i]
                Results['index = ' + str(index)] = Res

        if "MELTS" in Model:
            Results = stich(Results, multi = True, Model = Model)

        return Results



def path(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, fO2_buffer = None, fO2_offset = None, fluid_sat = None, Crystallinity_limit = None):
    '''
    Crystallisation calculations to be performed in parallel. Calculations may be either isobaric or isochoric.

    Parameters:
    ----------
    q: Multiprocessing Queue instance
        Queue instance to record the output variables

    index: int
        index of the calculation in the master code (e.g., position within a for loop) to aid indexing results after calculations are complete.

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
        If a specified pressure path is to be used, P_path_bar will be used to determine the P at each step of the model. If 2D, this indicates that multiple calculations with different P_path_bar arrays are to be performed.

    P_start_bar: float
        Initial pressure used for path calculations.

    P_end_bar: float
        Final pressure in crystallisation calculations.

    dp_bar: float
        Pressure increment during crystallisation calculations..

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio.

    H2O_Liq: float or np.ndarray
        H2O content of the initial bulk composition.

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

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    index: int
        index of the calculation

    '''

    Results = {}
    if "MELTS" in Model:
        try:
            Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit)
            q.put([Results, index])
        except:
            q.put([])

        return

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        Results = MM.path(comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index])
        return

