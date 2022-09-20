import numpy as np
import pandas as pd
# import Thermobar as pt
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.MELTS import *
try:
    from pyMELTScalc.Holland import *
except:
    pass
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.notebook import tqdm, trange

def multi_iso_crystallise(cores = None, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, isochoric = None, find_liquidus = None, Print_suppress = None):
    '''
    Carry out multiple crystallisation calculations in parallel. Allows isobaric or isochoric calculations to be performed. All temperature inputs/outputs are reported in degrees celcius and pressure is reported in bars.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default is 4.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid:
        If True, solid phases will be removed from the system at the end of each crystallisation step. Default False.

    Frac_fluid:
        If True, fluid phases will be removed from the system at the end of each crystallisation step. Default False.

    T_start_C: float
        Initial temperature used for crystallisation calculations.

    T_end_C: float
        Final temperature in crystallisation calculations. Default = 750.

    dt_C: float
        Temperature increment during crystallisation calculations. Default = 2.

    P_path_bar: float or np.ndarray
        Single pressure or an array of pressures equal to the number of calculations to be performed. Specifies the pressure of calculation (bar).

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, Fe3Fet_Liq must be a float and will set the Fe redox state in the initial composition. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable.

    H2O_Liq: float or np.ndarray
        H2O content of the initial melt phase. If type(comp) == dict, H2O_Liq must be a float. If comp is a pd.DataFrame, a single H2O_Liq value may be passes (float) and will be used as the initial melt H2O content for all starting compositions. Alternatively, if an array of H2O_Liq values is passed, equal to the number of compositions specified in comp, a different initial melt H2O value will be passed for each sample. If None, H2O_Liq must be specified in the comp variable.

    isochoric: True/False
        If True, the volume of the system will be held constant instead of the pressure. Default is False.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    Print_suppress: True/False
        If True, print messages concerning the status of the thermodynamic calculations will not be displayed.
    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.

    '''

    P_bar = P_path_bar

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

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition.")

    if cores is None:
        cores = 4

    if dt_C is None and T_path_C is None:
        dt_C = 2

    # specify the number of calculations to be performed in each sequence
    if type(P_bar) == np.ndarray:
        A = len(P_bar)//cores
        B = len(P_bar) % cores
    elif type(P_bar) != np.ndarray and type(comp) == dict:
        A = 1
        B = 0
    else:
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq']) % cores

    Group = np.zeros(A) + cores
    if B > 0:
        Group = np.append(Group, B)

    qs = []
    q = Queue()

    # perform calculation if only 1 calculation is specified
    if type(comp) == dict and type(P_bar) != np.ndarray:
        if Print_suppress is None:
            print("Running " + Model + " calculation...", end = "", flush = True)
            s = time.time()

        p = Process(target = iso_crystallise, args = (q, 1),
                    kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                            'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                            'P_path_bar': P_bar, 'isochoric': isochoric, 'find_liquidus': find_liquidus})

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
        for j in tqdm(range(len(Group))):
            ps = []

            if Print_suppress is None:
                print("Running " + Model + " calculations " + str(int(cores*j)) + " to " + str(int(cores*j) + Group[j] - 1) + " ...", end = "", flush = True)
                s = time.time()

            for i in range(int(cores*j), int(cores*j + Group[j])):
                if type(comp) == dict and type(P_bar) == np.ndarray:
                    p = Process(target = iso_crystallise, args = (q, i),
                                kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_path_bar': P_bar[i], 'isochoric': isochoric, 'find_liquidus': find_liquidus})
                elif type(comp) != dict and type(P_bar) != np.ndarray:
                    p = Process(target = iso_crystallise, args = (q, i),
                                kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_path_bar': P_bar, 'isochoric': isochoric, 'find_liquidus': find_liquidus})
                else:
                    p = Process(target = iso_crystallise, args = (q, i),
                                kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_path_bar': P_bar[i], 'isochoric': isochoric, 'find_liquidus': find_liquidus})

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

        if type(P_bar) == np.ndarray:
            Results = {}
            for i in range(len(qs)):
                if len(qs[i]) > 0:
                    Res, index = qs[i]
                    Results['P = ' + str(round(P_bar[index],2)) + ' bars'] = Res
        else:
            Results = {}
            for i in range(len(qs)):
                if len(qs[i]) > 0:
                    Res, index = qs[i]
                    Results['index = ' + str(index)] = Res

        Results = stich(Results, multi = True, Model = Model)

        return Results



def iso_crystallise(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, isochoric = None, find_liquidus = None):
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

    T_start_C: float
        Initial temperature used for crystallisation calculations.

    T_end_C: float
        Final temperature in crystallisation calculations. Default = 750.

    dt_C: float
        Temperature increment during crystallisation calculations. Default = 2.

    P_path_bar: float
         Specifies the pressure of calculation (bar).

    isochoric: True/False
        If True, the volume of the system will be held constant instead of the pressure. Default is False.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    index: int
        index of the calculation

    '''

    Results = {}
    if "MELTS" in Model:
        Results = crystallise_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, isochoric = isochoric, find_liquidus = find_liquidus)
        q.put([Results, index])
        return

    if Model == "Holland":
        Results = crystallise_holland(comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, isochoric = isochoric, find_liquidus = find_liquidus)
        q.put([Results, index])
        return

