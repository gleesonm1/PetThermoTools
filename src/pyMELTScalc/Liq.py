import numpy as np
import pandas as pd
from pyMELTScalc.Barom import *
from pyMELTScalc.Crystallise import *
try:
    from pyMELTScalc.Holland import *
except:
    pass
from pyMELTScalc.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange

def findLiq_multi(cores = None, Model = None, comp = None, T_initial_C = None, P_bar = None):
    '''
    Carry out multiple findLiq calculations in parallel.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default is 4.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: np.ndarray or pd.DataFrame
        Matrix containing all oxide values required for the calculations.

    T_initial_C: float or np.ndarray
        Initial 'guess' temperature for findLiq calculations (degrees C).

    P_bar: float or np.ndarray
        Single pressure or array of pressures equal to length of compositions in bulk. Specifies the pressure of calculation (bar).

    Returns:
    ----------
    T_Liq_C: np.ndarray
        Array of liquidus temperatures.

    H2O: np.ndarray
        Array of melt H2O contents at the liquidus.
    '''

    T_Liq = np.zeros(len(comp['SiO2_Liq'].values))
    T_in = np.zeros(len(comp['SiO2_Liq'].values))
    H2O_melt = np.zeros(len(comp['SiO2_Liq'].values))
    index = np.zeros(len(comp['SiO2_Liq'].values)) - 1

    if Model is None:
        Model == "MELTSv1.0.1"

    if cores is None:
        cores = 4

    if type(P_bar) == np.ndarray:
        A = len(P_bar)//cores
        B = len(P_bar) % cores
    else:
        A = len(comp['SiO2'])//cores
        B = len(comp['SiO2']) % cores

    Group = np.zeros(A) + cores
    Group = np.append(Group, B)

    qs = []
    q = Queue()
    for j in tqdm(range(len(Group))):
        ps = []
        for i in range(int(cores*j), int(cores*j + Group[j])):
            if type(comp) == np.ndarray:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': list(comp[:,i])})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': list(comp[:,i])})

            elif type(comp) == pd.core.frame.DataFrame:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict()})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict()})

            ps.append(p)
            p.start()

        for p in ps:
            try:
                ret = q.get(timeout = 90)
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

    for i in range(len(qs)):
        if len(qs[i])>0:
            T_Liq[i], H2O_melt[i], index[i], T_in[i] = qs[i]

    H2O = np.zeros(len(T_Liq))
    T_Liq_C = np.zeros(len(T_Liq))

    for i in range(len(index)):
        if len(T_Liq[index == i]) > 0:
            T_Liq_C[i] = T_Liq[index == i]
            H2O[i] = H2O_melt[index == i]

    return T_Liq_C, H2O


def findLiq(q, index,*, Model = None, P_bar = None, T_initial_C = None, comp = None):
    '''
    Searches for the liquidus of a melt at the pressure specified. Multiple instances of findLiq are typically initiated in parallel.

    Parameters:
    ----------
    q: Multiprocessing Queue instance
        Queue instance to place the output variables in

    index: int
        index of the calculation in the master code (e.g., position within a for loop) to aid indexing results after calculations are complete.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified by additing "v1.0.1", "v1.1.0", "v1.2.0", or "p" to "MELTS". Default "v.1.0.1".

    P_bar: float
        Specifies the pressure of the calculation (bar).

    T_initial_C: float
        Initial 'guess' temperature for findLiq calculations (degrees C).

    comp: dict
        Dictionary containing all oxide values required for the calculations.

    Returns:
    ----------
    T_Liq: float
        Estimated liquidus temperatures.

    H2O_Melt: float
        Melt water content at the liquidus

    index: int
        index of the calculation

    T_in: float
        Initial temperature of the calculation - used to ensure indexing is working correctly.
    '''

    T_Liq = 0
    T_in = T_initial_C
    H2O_Melt = 0

    if "MELTS" in Model:
        try:
            T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_bar, Model = Model, T_C_init = T_initial_C, comp = comp)
            q.put([T_Liq, H2O_Melt, index, T_in])
            return
        except:
            q.put([T_Liq, H2O_Melt, index, T_in])
            return

    if Model == "Holland":
        try:
            T_Liq = findLiq_holland(P_bar = P_bar, T_C_init = T_initial_C, comp = comp)
            q.put([T_Liq, H2O_Melt, index, T_in])
            return
        except:
            q.put([T_Liq, H2O_Melt, index, T_in])
            return

