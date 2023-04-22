import numpy as np
import pandas as pd
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.Path import *
from pyMELTScalc.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange
import time

def phaseDiagram_calc(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, i_max = 25):
    """
    Calculate phase diagrams for igneous systems (rocks and magmas).

    Parameters:
    -----------
    cores : int, optional
        The number of CPU cores to use for parallel processing. Default is None, which
        sets the number of cores to the number of CPUs available on the machine.
    Model : str, optional
        The name of the thermodynamic model to use. Default is None, which sets the model
        to "MELTSv1.0.2".
    bulk : dict or pandas.Series
        The bulk composition of the system. If passed as a pandas.Series, it will first
        be converted to a dict. Default is None.
    T_C : array-like, optional
        The array of temperature values to use for the phase diagram, in degrees Celsius.
        Default is None.
    P_bar : array-like, optional
        The array of pressure values to use for the phase diagram, in bars. Default is None.
    T_min_C : float, optional
        The minimum temperature value to use for the phase diagram, in degrees Celsius.
        Default is None.
    T_max_C : float, optional
        The maximum temperature value to use for the phase diagram, in degrees Celsius.
        Default is None.
    T_num : int, optional
        The number of temperature values to use for the phase diagram. Default is None.
    P_min_bar : float, optional
        The minimum pressure value to use for the phase diagram, in bars. Default is None.
    P_max_bar : float, optional
        The maximum pressure value to use for the phase diagram, in bars. Default is None.
    P_num : int, optional
        The number of pressure values to use for the phase diagram. Default is None.
    Fe3Fet_Liq : float, optional
        The Fe3+/Fetot ratio for the liquid phase. Default is None.
    H2O_Liq : float, optional
        The water content of the liquid phase, in wt%. Default is None.
    fO2_buffer : str, optional
        The name of the oxygen buffer to use for the phase diagram. Default is None.
    fO2_offset : float, optional
        The offset to apply to the fO2 buffer value. Default is None.
    i_max : int, optional
        The maximum number of attempts to make at calculating the phase diagram. Default is 25.

    Returns:
    --------
    pandas.DataFrame
        A dataframe containing the phase diagram results.
    """

    comp = bulk.copy()

    if cores is None:
        cores = multiprocessing.cpu_count()

    if Model is None:
        Model = "MELTSv1.0.2"

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        import julia
        from julia import MAGEMinCalc

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

    if T_min_C is not None:
        T_C = np.linspace(T_min_C, T_max_C, T_num)

    if P_min_bar is not None:
        P_bar = np.linspace(P_min_bar, P_max_bar, P_num)

    T, P = np.meshgrid(T_C, P_bar)

    T_flat = np.round(T.flatten(),2)
    P_flat = np.round(P.flatten(),2)

    A = len(P_flat)//cores + 1

    c = 0
    j = 0
    while len(T_flat)>1:
        if j > i_max:
            break

        print('Attempt ' + str(j))

        qs = []
        q = Queue()
        ps = []
        j = j + 1
        for i in range(cores):
            T_path_C = T_flat[i*A:(i+1)*A]
            P_path_bar = P_flat[i*A:(i+1)*A]

            if len(T_path_C) > 100:
                T_path_C = T_path_C[:99]
                P_path_bar = P_path_bar[:99]

            p = Process(target = path, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_path_C': T_path_C, 'P_path_bar': P_path_bar, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

            ps.append(p)
            p.start()

        TIMEOUT = 180
        start = time.time()

        for p in ps:
            if time.time() - start < TIMEOUT - 10:
                try:
                    ret = q.get(timeout = TIMEOUT - (time.time()-start) + 10)
                except:
                    ret = []
            else:
                try:
                    ret = q.get(timeout = 10)
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

        # TIMEOUT = 150 #+ #0.5*len(T_flat)
        # start = time.time()
        # for p in ps:
        #     ret = []
        #     if time.time() - start <= TIMEOUT:
        #         while time.time() - start <= TIMEOUT:
        #             try:
        #                 ret = q.get(timeout = 2)
        #                 break
        #             except:
        #                 continue
        #
        #             # if p.is_alive():
        #             #     time.sleep(.1)
        #             # else:
        #             #     p.terminate()
        #             #     p.join()
        #             #     break
        #
        #         # if p.is_alive():
        #     else:
        #         try:
        #             ret = q.get(timeout = 5)
        #         except:
        #             ret =[]
        #
        #     p.terminate()
        #     p.join()

            # else:
            #     p.terminate()
            #     p.join()

            # try:
            #     ret = q.get(timeout = 5)
            # except:
            #     ret = []

            # qs.append(ret)
            # if p.is_alive():
            #     while time.time() - start <= TIMEOUT:
            #         if not p.is_alive():
            #             p.terminate()
            #             p.join()
            #             break
            #         time.sleep(.1)
            #     else:
            #         p.terminate()
            #         p.join()
            # else:
            #     p.terminate()
            #     p.join()
            #
            # try:
            #     ret = q.get(timeout = 2)
            # except:
            #     ret = []

            # qs.append(ret)

        for p in ps:
            p.kill()

        Results = {}
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Res, index = qs[i]
                Results['index = ' + str(index)] = Res

        if "MELTS" in Model:
            Results = stich(Results, multi = True, Model = Model)

        for i in Results.keys():
            if c == 0:
                Combined = Results[i]['All'].copy()
                c = 1
            else:
                Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)

        flat = np.round(np.array([T.flatten(), P.flatten()]).T,2)
        Res_flat = np.round(np.array([Combined['T_C'], Combined['P_bar']]).T,2)
        new_flat = flat[np.where((flat[:, None] == Res_flat).all(-1).any(-1) == False)]

        if len(new_flat.T[0]) > 1:
            T_flat = new_flat.T[0][:]
            P_flat = new_flat.T[1][:]
            if (j % 2) != 0:
                T_flat = np.flip(T_flat)
                P_flat = np.flip(P_flat)

            A = len(P_flat)//cores + 1
        else:
            break

    flat = np.round(np.array([T.flatten(), P.flatten()]).T,2)
    Res_flat = np.round(np.array([Combined['T_C'].values, Combined['P_bar'].values]).T,2)
    new_flat = flat[np.where((flat[:, None] == Res_flat).all(-1).any(-1) == False)]

    if np.shape(new_flat)[0] > 0.0:
        # df = pd.DataFrame(columns = ['T_C', 'P_bar'])
        # for i in range(len(new_flat)):
        #     df.loc[i] = new_flat[i]
        A = np.zeros((np.shape(new_flat)[0], np.shape(Combined.values)[1]))
        A[:,:2] = new_flat

        B = Combined.values

        C = np.concatenate((A, B))

        Combined = pd.DataFrame(columns = list(Combined.keys()), data = C)
        # Combined = pd.concat([Combined, df], axis = 0, ignore_index = True)

    Combined['T_C'] = np.round(Combined['T_C'], 2)
    Combined['P_bar'] = np.round(Combined['P_bar'], 2)
    Combined = Combined.sort_values(['T_C', 'P_bar'])
    Combined = Combined.reset_index(drop = True)

    return Combined


def phaseDiagram_eq(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, number_max = 50):
    """
    Calculates the phase diagram for a given bulk composition over a
    range of temperatures and pressures.

    Parameters:
    -----------
    cores : int or None, optional
        The number of cores to use for multiprocessing, default is set to the number of logical processors available.
    Model : str or None, optional
        The thermodynamic model to use, either 'MELTSv1.0.2' or 'Holland'. Default is None.
    bulk : dict, pandas.core.series.Series
        The bulk composition to use. Can be a dictionary with the element names and mole
        fractions, a pandas series.
    T_C : array-like or None, optional
        The temperature range (in degrees Celsius) to calculate the phase diagram over. If
        T_min_C is given and T_C is None, this array is generated automatically. Default is None.
    P_bar : array-like or None, optional
        The pressure range (in bars) to calculate the phase diagram over. If P_min_bar is given
        and P_bar is None, this array is generated automatically. Default is None.
    T_min_C : float or None, optional
        The minimum temperature (in degrees Celsius) to calculate the phase diagram over. Default is None.
    T_max_C : float or None, optional
        The maximum temperature (in degrees Celsius) to calculate the phase diagram over. Default is None.
    T_num : int or None, optional
        The number of temperature steps to calculate the phase diagram over. Default is None.
    P_min_bar : float or None, optional
        The minimum pressure (in bars) to calculate the phase diagram over. Default is None.
    P_max_bar : float or None, optional
        The maximum pressure (in bars) to calculate the phase diagram over. Default is None.
    P_num : int or None, optional
        The number of pressure steps to calculate the phase diagram over. Default is None.
    Fe3Fet_Liq : float or None, optional
        The Fe3+/Fet ratio of the liquid for the bulk composition. Default is None.
    H2O_Liq : float or None, optional
        The H2O content of the liquid for the bulk composition. Default is None.
    fO2_buffer : str or None, optional
        The oxygen buffer to use for MELTS. Default is None.
    fO2_offset : float or None, optional
        The oxygen buffer offset to use for MELTS. Default is None.
    number_max : int, optional
        The maximum number of calculations to perform in a single Holland model calculation.
        If there are more than this number of calculations to perform, multiprocessing is used.
        Default is 50.

    Returns:
    --------
    pandas.core.frame.DataFrame
        A dataframe containing the equilibrium phase assemblage(s) for the given bulk
        composition over the range of temperatures and pressures specified.
    """
    comp = bulk.copy()

    if cores is None:
        cores = multiprocessing.cpu_count()

    if Model is None:
        Model = "MELTSv1.0.2"

    if Model == "Holland":
        import pyMAGEMINcalc as MM

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

    if T_min_C is not None:
        T_C = np.linspace(T_min_C, T_max_C, T_num)

    if P_min_bar is not None:
        P_bar = np.linspace(P_min_bar, P_max_bar, P_num)

    T, P = np.meshgrid(T_C, P_bar)

    T_flat = np.round(T.flatten(),2)
    P_flat = np.round(P.flatten(),2)

    if "MELTS" in Model:
        A = len(P_flat)//cores
        B = len(P_flat) % cores
        if A > 0:
            Group = np.zeros(A) + cores
            if B > 0:
                Group = np.append(Group, B)
        else:
            Group = np.array([B])

        qs = []
        q = Queue()
        c = 0
        #Combined = pd.DataFrame(columns = ['T_C', 'P_bar'], data = np.zeros((1,2)))
        for j in tqdm(range(len(Group))):
            ps = []
            for i in range(int(cores*j), int(cores*j + Group[j])):
                p = Process(target = equilibrate, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_C': T_flat[i], 'P_bar': P_flat[i], 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

                ps.append(p)
                p.start()

            TIMEOUT = 60 #+ #0.5*len(T_flat)
            start = time.time()
            for p in ps:
                if TIMEOUT  - (time.time() - start) > 10:
                    try:
                        ret = q.get(timeout = TIMEOUT - (time.time() - start))
                    except:
                        ret = []

                else:
                    try:
                        ret = q.get(timeout = 10)
                    except:
                        ret = []

                qs.append(ret)

            for p in ps:
                p.kill()

        Results = {}
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Res, index = qs[i]
                Results['index = ' + str(index)] = Res

        if "MELTS" in Model:
            Results = stich(Results, multi = True, Model = Model)

        for i in Results.keys():
            if c == 0:
                if "MELTS" in Model:
                    Combined = Results[i]['All'].copy()
                else:
                    Combined = Results[i].copy()
                c = 1
            else:
                if "MELTS" in Model:
                    Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)
                else:
                    try:
                        Combined = pd.concat([Combined, Results[i]], axis = 0, ignore_index = True)
                    except:
                        continue

    if Model == "Holland":
        Combined = pd.DataFrame()
        if len(T_flat) < number_max:
            c = 0
            for i in tqdm(range(len(T_flat))):
                try:
                    Results = MM.equilibrate(Model = Model, P_bar = P_flat[i], T_C = T_flat[i], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
                    Combined = pd.concat([Combined, Results], axis = 0, ignore_index = True)
                except:
                    pass
        else:
            s = len(T_flat)//number_max
            A = s//cores
            B = s % cores

            if A > 0:
                Group = np.zeros(A) + cores
            if B > 0:
                if A > 0:
                    Group = np.append(Group, B)
                else:
                    Group = np.array([B])

            qs = []
            q = Queue()
            for j in tqdm(range(len(Group))):
                ps = []
                for i in range(int(cores*j), int(cores*j + Group[j])):
                    T2 = T_flat[i*number_max:i*number_max+number_max]
                    P2 = P_flat[i*number_max:i*number_max+number_max]
                    p = Process(target = equilibrate, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_C': T2, 'P_bar': P2, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

                    ps.append(p)
                    p.start()

                for p in ps:
                    try:
                        ret = q.get()
                    except:
                        ret = []

                    qs.append(ret)

                for p in ps:
                    p.kill()

            for i in range(len(qs)):
                if len(qs[i]) > 0:
                    Res, index = qs[i]
                    Combined = pd.concat([Combined, Res], axis = 0, ignore_index = True)


    if len(Combined['T_C']) < len(T_flat):
        flat = np.round(np.array([T.flatten(), P.flatten()]).T,2)
        Res_flat = np.round(np.array([Combined['T_C'].values, Combined['P_bar'].values]).T,2)
        new_flat = flat[np.where((flat[:, None] == Res_flat).all(-1).any(-1) == False)]

        # Combined['T_C'] = np.round(Combined['T_C'].values, 2)
        # Combined['P_bar'] = np.round(Combined['P_bar'].values, 2)

        for i in range(len(new_flat)):
            df = pd.DataFrame(columns = ['T_C', 'P_bar'])
            df.loc[0] = new_flat[i]

            Combined = pd.concat([Combined, df], axis = 0, ignore_index = True)

    return Combined


def equilibrate(q, index,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):

    if "MELTS" in Model:
        Results = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index])
        return

    if "Holland" in Model:
        import pyMAGEMINcalc as MM
        from tqdm.notebook import tqdm, trange
        Combined = pd.DataFrame()
        for i in range(len(T_C)):
            try:
                Results = MM.equilibrate(Model = Model, P_bar = P_bar[i], T_C = T_C[i], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
                Combined = pd.concat([Combined, Results], axis = 0, ignore_index = True)
            except:
                pass

        q.put([Combined, index])
        return
