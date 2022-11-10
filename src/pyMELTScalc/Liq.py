import numpy as np
import pandas as pd
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.MELTS import *
# try:
#     from pyMELTScalc.Holland import *
# except:
#     pass
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange

def findLiq_multi(cores = None, Model = None, comp = None, T_initial_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
    '''
    Carry out multiple findLiq calculations in parallel.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default will be determined using Multiprocessing.cpu_count().

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: pd.DataFrame
        Matrix containing all oxide values required for the calculations.

    T_initial_C: float or np.ndarray
        Initial 'guess' temperature for findLiq calculations (degrees C).

    P_bar: float or np.ndarray
        Single pressure or array of pressures equal to length of compositions in bulk. Specifies the pressure of calculation (bar).

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, and type(Fe3Fet_Liq) == np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their Fe3Fet_Liq value. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable or an oxygen fugacity buffer must be chosen.

    H2O_Liq: float or np.ndarray
        H2O content of the initial melt phase. If type(comp) == dict, and type(H2O_Liq) = np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their H2O_Liq value. If comp is a pd.DataFrame, a single H2O_Liq value may be passes (float) and will be used as the initial melt H2O content for all starting compositions. Alternatively, if an array of H2O_Liq values is passed, equal to the number of compositions specified in comp, a different initial melt H2O value will be passed for each sample. If None, H2O_Liq must be specified in the comp variable.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float or np.ndarray
        Offset from the buffer spcified in fO2_buffer (log units).

    Returns:
    ----------
    T_Liq_C: np.ndarray
        Array of liquidus temperatures.

    H2O: np.ndarray
        Array of melt H2O contents at the liquidus.
    '''
    if Model is None:
        Model = "MELTSv1.0.2"

    if Model == "Holland":
        import pyMAGEMINcalc as MM

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

    if type(comp) != dict or type(P_bar) == np.ndarray:
        if type(comp) != dict:
            T_Liq = np.zeros(len(comp['SiO2_Liq'].values))
            T_in = np.zeros(len(comp['SiO2_Liq'].values))
            H2O_melt = np.zeros(len(comp['SiO2_Liq'].values))
            index = np.zeros(len(comp['SiO2_Liq'].values)) - 1
        else:
            T_Liq = np.zeros(len(P_bar))
            T_in = np.zeros(len(P_bar))
            H2O_melt = np.zeros(len(P_bar))
            index = np.zeros(len(P_bar)) - 1
    else:
        T_Liq = 0
        T_in = 0
        H2O_melt = 0
        index = 0

    if T_initial_C is None:
        if type(comp) != dict or type(P_bar) == np.ndarray:
            if type(comp) != dict:
                T_initial_C = np.zeros(len(comp['SiO2_Liq'].values)) + 1300.0
            else:
                T_initial_C = np.zeros(len(P_bar)) + 1300.0
        else:
            T_initial_C = np.array([1300.0])


    if cores is None:
        cores = multiprocessing.cpu_count()

    if type(P_bar) == np.ndarray:
        A = len(P_bar)//cores
        B = len(P_bar) % cores
    elif type(P_bar) != np.ndarray and type(comp) == dict:
        A = 0
        B = 1
    else:
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq']) % cores

    if A > 0:
        Group = np.zeros(A) + cores
        if B > 0:
            Group = np.append(Group, B)
    else:
        Group = np.array([B])

    qs = []
    q = Queue()
    for j in tqdm(range(len(Group))):
        ps = []
        for i in range(int(cores*j), int(cores*j + Group[j])):
            if type(comp) == dict:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
            else:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

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

    if type(comp) != dict or type(P_bar) == np.ndarray:
        for i in range(len(qs)):
            if len(qs[i])>0:
                T_Liq[i], H2O_melt[i], index[i], T_in[i] = qs[i]

        H2O = np.zeros(len(T_Liq))
        T_Liq_C = np.zeros(len(T_Liq))

        for i in range(len(index)):
            if len(T_Liq[index == i]) > 0:
                T_Liq_C[i] = T_Liq[index == i]
                H2O[i] = H2O_melt[index == i]
    else:
        T_Liq_C, H2O, index, T_in = qs[0]

    return T_Liq_C, H2O


def findLiq(q, index,*, Model = None, P_bar = None, T_initial_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
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
            T_Liq, H2O_Melt = findLiq_MELTS(P_bar = P_bar, Model = Model, T_C_init = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
            q.put([T_Liq, H2O_Melt, index, T_in])
            return
        except:
            q.put([T_Liq, H2O_Melt, index, T_in])
            return

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        #try:
        T_Liq = MM.findLiq(P_bar = P_bar, T_C_init = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([T_Liq, H2O_Melt, index, T_in])
        return
        #except:
        #    q.put([T_Liq, H2O_Melt, index, T_in])
        #    return

