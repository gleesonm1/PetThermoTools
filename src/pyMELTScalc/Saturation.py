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

def findSatPressure(cores = None, Model = None, bulk = None, T_C_init = None, P_bar_init = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
    """
    Calculates the saturation pressure of a specified composition in the liquid or melt phase. The function will return the saturation pressure of the liquid or melt as a pandas dataframe. If the saturation pressure cannot be calculated, the function will return an empty dataframe.

    Parameters:
    cores (int): The number of CPU cores to use for parallel processing. If left as None, all available CPU cores will be used.
    Model (str): The geochemical model to be used. Options include "MELTSv1.2.0" (default) or "MELTSv1.1.0".
    bulk (dict or pandas.core.frame.DataFrame): A dictionary or pandas dataframe containing the starting composition of the liquid or melt. Elements should be specified in weight percent.
    T_C_init (float, int, or np.ndarray): The starting temperature for the calculation in degrees Celsius. If left as None, the default value of 1200 C will be used.
    P_bar_init (float, int, or np.ndarray): The starting pressure for the calculation in bar. If left as None, the default value of 10000 bar will be used.
    Fe3Fet_Liq (float): The initial concentration of ferric iron in the liquid or melt in weight percent. If left as None, the value will be taken from the bulk composition.
    H2O_Liq (float): The initial concentration of water in the liquid or melt in weight percent. If left as None, the value will be taken from the bulk composition.
    fO2_buffer (float): The  oxygen fugacity buffer to be used in the calculation. If left as None, no buffer will be applied
    fO2_offset (float, int, or np.ndarray): The oxygen fugacity offset (log units) from the specified buffer to be applied to the calculation. If left as None, no offset will be applied. If a list is provided, each offset will be applied to the corresponding element in the bulk composition.

    Returns:
    Res (dict or pandas.core.frame.DataFrame): A dict or dataframe containing the chemical composition of the melt phase and the saturation pressure and temperature of the melt.

    Examples:

    Calculate the saturation pressure and temperature of a basaltic liquid, with initial conditions set at 1200 C and 10000 bar using all available CPU cores and the MELTSv1.2.0 model.
    >>> findSatPressure(bulk = {'SiO2_Liq': 50, 'TiO2_Liq': 1, 'Al2O3_Liq': 15, 'FeOt_Liq': 10, 'MnO_Liq': 0.5, 'MgO_Liq': 8, 'CaO_Liq': 10, 'Na2O_Liq': 2, 'K2O_Liq': 1, 'P2O5_Liq': 0.5, 'H2O_Liq': 0.5, 'Fe3Fet_Liq': 0.15}, T_C_init = 1200, P_bar_init = 10000)

    """

    # set default values if required
    if Model is None:
        Model == "MELTSv1.2.0"

    comp = bulk.copy()

    if T_C_init is None:
        T_C_init = 1200

    if P_bar_init is None:
        P_bar_init = 10000

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

    qs = []
    q = Queue()

    # specify the number of calculations to be performed in each sequence
    if type(comp) == dict:
        p = Process(target = satP, args = (q, 0),
                    kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init, 'P_bar_init': P_bar_init,
                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

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

        if len(ret) > 0:
            Results, index = ret

        Res = Results.copy()

        return Res

    else:
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq'])%cores

        if A > 0:
            Group = np.zeros(A) + cores
            if B > 0:
                Group = np.append(Group, B)
        else:
            Group = np.array([B])

        L = np.sum(Group)
        if type(P_bar_init) == float or type(P_bar_init) == int:
            P_bar_init = np.zeros(int(L)) + P_bar_init
        if type(T_C_init) == float or type(T_C_init) == int:
            T_C_init = np.zeros(int(L)) + T_C_init
        if fO2_offset is None:
            fO2_offset = [None] * int(L)
        else:
            if type(fO2_offset) == float or type(fO2_offset) == int:
                fO2_offset = np.zeros(int(L)) + fO2_offset

        for j in tqdm(range(len(Group))):
            ps = []

            for i in range(int(cores*j), int(cores*j + Group[j])):
                p = Process(target = satP, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'T_C_init': T_C_init[i], 'P_bar_init': P_bar_init[i],
                                    'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i]})

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

        Res = pd.DataFrame(data = np.zeros((int(L),15)), columns = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq', 'P_bar', 'T_Liq'])
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Results, index = qs[i]
                Res.loc[index] = Results

        return Res



def satP(q, index, *, Model = None, comp = None, T_C_init = None, P_bar_init = None, fO2_buffer = None, fO2_offset = None):
    """Find the saturation pressure for a given composition and temperature.

    This function calculates the volatile saturation pressure for a given composition using the MELTS models. The results are returned in the form of a tuple and added to the specified queue.

    Args:
        q (Queue): The queue to which the results should be added.
        index (int): The index of the calculation.
        Model (str, optional): The model to use for the calculation. Must include "MELTS".
        comp (dict): A dictionary of oxide names and their concentration (in wt%) in the composition.
        T_C_init (float, optional): The initial temperature in degrees Celsius.
        P_bar_init (float, optional): The initial pressure in bar.
        fO2_buffer (float, optional): The fO2 buffer.
        fO2_offset (float, optional): The fO2 offset.

    Returns:
        None
    """

    if "MELTS" in Model:
        Results = findSatPressure_MELTS(Model = Model, T_C_init = T_C_init, P_bar_init = P_bar_init, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index])
        return