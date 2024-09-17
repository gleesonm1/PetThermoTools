import numpy as np
import pandas as pd
from PetThermoTools.GenFuncs import *
from PetThermoTools.Plotting import *
from PetThermoTools.MELTS import *
# try:
#     from PetThermoTools.Holland import *
# except:
#     pass
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange

def equilibrate_multi(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, 
                      Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None,
                      timeout = None, merge_on = None):
    comp = bulk.copy()

    if Model is None:
        Model = "MELTSv1.0.2"

    if Model == "Holland":
        import pyMAGEMINcalc as MM

    if cores is None:
        cores = multiprocessing.cpu_count()

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


    if type(P_bar) != np.ndarray:
        if type(comp) == pd.core.frame.DataFrame:
            P_bar = np.zeros(len(comp['SiO2_Liq'])) + P_bar
        else:
            P_bar = np.array([P_bar])

    if type(T_C) != np.ndarray:
        if type(comp) == pd.core.frame.DataFrame:
            T_C = np.zeros(len(comp['SiO2_Liq'])) + T_C
        else:
            T_C = np.array([T_C])

    if type(comp) == pd.core.frame.DataFrame:
        A = len(P_bar)//cores
        B = len(P_bar) % cores

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
                p = Process(target = equilibrate, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_C': T_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

                ps.append(p)
                p.start()

            if timeout is None:
                TIMEOUT = 20
            else:
                TIMEOUT = timeout

            start = time.time()
            first = True
            for p in ps:
                if time.time() - start < TIMEOUT - 10:
                    try:
                        ret = q.get(timeout = TIMEOUT - (time.time()-start) + 10)
                    except:
                        ret = []
                else:
                    if first == True:
                        print('Timeout Reached - this likely indicates the calculation failed. \n You can try increasing the timeout limit using the "timeout" kwarg.')
                        first = False
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

        Output = {}
        Affinity = {}
        if len(qs) > 1:
            for i in range(len(qs)):
                if len(qs[i])>0:
                    Results, index = qs[i]

                    try:
                        Output[str(index)] = Results[0]
                        Affinity[str(index)] = Results[1]
                    except:
                        continue

            if "MELTS" in Model:
                Output = stich(Output, multi = True, Model = "MELTS")

            int_keys = list(map(int, Output.keys()))

            # Determine the size of the resulting DataFrame
            # max_index = max(int_keys)
            max_index = len(comp['SiO2_Liq'])

            # Initialize an empty list to maintain column order
            all_columns = []
            seen_columns = set()

            for df in Output.values():
                for col in df['All'].columns:
                    if col not in seen_columns:
                        all_columns.append(col)
                        seen_columns.add(col)

            Combined = pd.DataFrame(np.nan, index=range(max_index), columns=all_columns)

            # Populate the result DataFrame
            for key, df in Output.items():
                try:
                    Combined.loc[int(key)] = df['All'].iloc[0]
                except:
                    continue

            # Affinity combined
            int_keys = list(map(int, Affinity.keys()))

            # Determine the size of the resulting DataFrame
            # max_index = max(int_keys)
            max_index = len(comp['SiO2_Liq'])

            # Initialize an empty list to maintain column order
            all_columns = []
            seen_columns = set()

            for d in Affinity:
                for col in Affinity[d]:
                    if col not in seen_columns:
                        all_columns.append(col)
                        seen_columns.add(col)
                
                Affinity[d] = pd.Series(Affinity[d])

            Af_Combined = pd.DataFrame(np.nan, index=range(max_index), columns=all_columns)

            # Populate the result DataFrame
            for d in Affinity:
                try:
                    df = Affinity[d].to_frame().transpose()
                    Af_Combined.loc[int(d)] = df.iloc[0]
                except:
                    continue
            
        else:
            if len(qs[0]) > 0:
                Results, index = qs[0]
                Output = Results[0]
                Affinity = Results[1]
            
                if "MELTS" in Model:
                    Output = stich(Output, Model = Model)

    if type(comp) != pd.core.frame.DataFrame:
        A = len(P_bar)//cores
        B = len(P_bar) % cores

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
                p = Process(target = equilibrate, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_C': T_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

                ps.append(p)
                p.start()

            if timeout is None:
                TIMEOUT = 20
            else:
                TIMEOUT = timeout

            start = time.time()
            first = True
            for p in ps:
                if time.time() - start < TIMEOUT - 10:
                    try:
                        ret = q.get(timeout = TIMEOUT - (time.time()-start) + 10)
                    except:
                        ret = []
                else:
                    if first == True:
                        print('Timeout Reached - this likely indicates the calculation failed. \n You can try increasing the timeout limit using the "timeout" kwarg.')
                        first = False
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

        Output = {}
        Affinity = {}
        if len(qs) > 1:
            for i in range(len(qs)):
                if len(qs[i])>0:
                    Results, index = qs[i]

                    try:
                        Output[str(index)] = Results[0]
                        Affinity[str(index)] = Results[1]
                    except:
                        continue

            if "MELTS" in Model:
                Output = stich(Output, multi = True, Model = Model)

            # Output combined
            int_keys = list(map(int, Output.keys()))

            # Determine the size of the resulting DataFrame
            # max_index = max(int_keys)
            max_index = len(P_bar)

            # Initialize an empty list to maintain column order
            all_columns = []
            seen_columns = set()

            for df in Output.values():
                for col in df['All'].columns:
                    if col not in seen_columns:
                        all_columns.append(col)
                        seen_columns.add(col)

            Combined = pd.DataFrame(np.nan, index=range(max_index), columns=all_columns)

            # Populate the result DataFrame
            for key, df in Output.items():
                try:
                    Combined.loc[int(key)] = df['All'].iloc[0]
                except:
                    continue

            # Affinity combined
            int_keys = list(map(int, Affinity.keys()))

            # Determine the size of the resulting DataFrame
            # max_index = max(int_keys)
            max_index = len(P_bar)

            # Initialize an empty list to maintain column order
            all_columns = []
            seen_columns = set()

            for d in Affinity:
                for col in Affinity[d]:
                    if col not in seen_columns:
                        all_columns.append(col)
                        seen_columns.add(col)
                Affinity[d] = pd.Series(Affinity[d])

            Af_Combined = pd.DataFrame(np.nan, index=range(max_index), columns=all_columns)

            # Populate the result DataFrame
            for d in Affinity:
                try:
                    df = Affinity[d].to_frame().transpose()
                    Af_Combined.loc[int(d)] = df.iloc[0]
                except:
                    continue
        else:
            if len(qs[0]) > 0:
                Output, index = qs[0]
                Output = Results[0]
                Affinity = Results[1]
            
                if "MELTS" in Model:
                    Output = stich(Output, Model = Model)

    if merge_on is not None:
        if type(merge_on) == str:
            Combined.insert(0, merge_on, comp[merge_on])
            # Af_Combined.insert(0, merge_on, comp[merge_on])
        elif type(merge_on) == list:
            j = 0
            for i in merge_on:
                Combined.insert(j, i, comp[i])
                # Af_Combined.insert(j, i, comp[i])
                j = j + 1

    Af_Combined = Af_Combined.add_suffix('_affinity')
    Combined = pd.concat([Combined, Af_Combined], axis = 1)
    # Affinity = Af_Combined.copy()
    return Combined

def findCO2_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):

    comp = bulk.copy()

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
            T_Liq_C = np.zeros(len(comp['SiO2_Liq'].values))
            H2O_melt = np.zeros(len(comp['SiO2_Liq'].values))
            CO2_melt = np.zeros(len(comp['SiO2_Liq'].values))
            index = np.zeros(len(comp['SiO2_Liq'].values)) - 1
        else:
            T_Liq_C = np.zeros(len(P_bar))
            H2O_melt = np.zeros(len(P_bar))
            CO2_melt = np.zeros(len(P_bar))
            index = np.zeros(len(P_bar)) - 1
    else:
        T_Liq_C = 0
        H2O_melt = 0
        CO2_melt = 0
        index = 0

    if T_initial_C is None:
        if type(comp) != dict or type(P_bar) == np.ndarray:
            if type(comp) != dict:
                T_initial_C = np.zeros(len(comp['SiO2_Liq'].values)) + 1300.0
            else:
                T_initial_C = np.zeros(len(P_bar)) + 1300.0
        else:
            T_initial_C = np.array([1300.0])

    if type(T_initial_C) == int or type(T_initial_C) == float:
        T_initial_C = np.array([T_initial_C])

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
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                else:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
            else:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                else:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

            ps.append(p)
            p.start()

        TIMEOUT = 240

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

    if type(comp) != dict or type(P_bar) == np.ndarray:
        for i in range(len(qs)):
            if len(qs[i])>0:
                T_Liq_C[i], H2O_melt[i], CO2_melt[i], index[i] = qs[i]

        T_Liq = np.zeros(len(T_Liq_C))
        H2O = np.zeros(len(T_Liq_C))
        CO2 = np.zeros(len(T_Liq_C))

        for i in range(len(index)):
            if len(CO2[index == i]) > 0:
                T_Liq[i] = T_Liq_C[index == i]
                H2O[i] = H2O_melt[index == i]
                CO2[i] = CO2_melt[index == i]
    else:
        T_Liq, H2O, CO2, index = qs[0]

    return T_Liq, H2O, CO2

def findLiq_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None, Affinity = False):
    '''
    Perform multiple liquidus temperature (findLiq) calculations in parallel.

    Parameters:
    ----------
    cores: int, optional
        Number of processes to run in parallel. If not provided, the number of available CPU cores is used.

    Model: str, optional
        The thermodynamic model to use for calculations: "MELTS..." or "Holland". 
        If "MELTS..." is selected, the version needs to be specified (e.g., "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS").
        Default is "MELTSv1.0.2".

    bulk: pd.DataFrame
        A dataframe containing the oxide compositions required for the calculations.

    T_initial_C: float or np.ndarray, optional
        Initial guess for the temperature in Celsius for the findLiq calculations. If not provided, the default is 1300°C.

    P_bar: float or np.ndarray
        Pressure (in bars) for the calculations. Can be a single value or an array with length matching the number of compositions in the bulk.

    Fe3Fet_Liq: float or np.ndarray, optional
        Fe³⁺/total Fe ratio in the liquid. If the composition is provided as a dictionary, this can be an array to vary the Fe³⁺/Fe ratio across compositions. 
        If the composition is a DataFrame, a single value or array (matching the number of compositions) can be used. 
        If not provided, the redox state must be defined elsewhere (e.g., via the fO2 buffer or directly in the bulk variable).

    H2O_Liq: float or np.ndarray, optional
        Initial water content in the liquid phase. Can be a single value for all compositions or an array with values for each composition. 
        If not provided, H2O content must be included in the composition.

    CO2_Liq: float or np.ndarray, optional
        Initial CO₂ content in the liquid phase. Like H2O_Liq, this can be a single value or an array matching the number of compositions.

    fO2_buffer: str, optional
        Oxygen fugacity buffer to control oxidation states during crystallization or decompression. Options include "FMQ" (Fayalite-Magnetite-Quartz) or "NNO" (Nickel-Nickel Oxide).

    fO2_offset: float or np.ndarray, optional
        Offset from the specified fO2 buffer (in log units). Can be a single value or an array.

    Affinity: bool, optional
        If True, the function will return additional information related to chemical affinity (interaction between phases) during crystallization. Default is False.

    Returns:
    ----------
    Res: pd.DataFrame
        DataFrame containing liquidus temperature, liquidus phase, fluid saturation, and normalized melt chemistry.

    Af_Combined: pd.DataFrame, optional
        If Affinity is True, an additional DataFrame containing chemical affinity data for each sample is returned.
    '''
    comp = bulk.copy()

    if Model is None:
        Model = "MELTSv1.0.2"

    if Model == "Holland":
        import pyMAGEMINcalc as MM

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq)

    if type(comp) != dict or type(P_bar) == np.ndarray:
        if type(comp) != dict:
            T_Liq = np.zeros(len(comp['SiO2_Liq'].values))
            T_in = np.zeros(len(comp['SiO2_Liq'].values))
            index = np.zeros(len(comp['SiO2_Liq'].values)) - 1
        else:
            T_Liq = np.zeros(len(P_bar))
            T_in = np.zeros(len(P_bar))
            index = np.zeros(len(P_bar)) - 1
    else:
        T_Liq = 0
        T_in = 0
        index = 0

    if T_initial_C is None:
        if type(comp) != dict or type(P_bar) == np.ndarray:
            if type(comp) != dict:
                T_initial_C = np.zeros(len(comp['SiO2_Liq'].values)) + 1300.0
            else:
                T_initial_C = np.zeros(len(P_bar)) + 1300.0
        else:
            T_initial_C = np.array([1300.0])
    else:
        if type(T_initial_C) == int or type(T_initial_C) == float:
            if type(comp) != dict or type(P_bar) == np.ndarray:
                if type(comp) != dict:
                    T_initial_C = np.zeros(len(comp['SiO2_Liq'].values)) + T_initial_C
                else:
                    T_initial_C = np.zeros(len(P_bar)) + T_initial_C
            else:
                T_initial_C = np.array([T_initial_C])
            


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
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Affinity': Affinity})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Affinity': Affinity})
            else:
                if type(P_bar) == np.ndarray:
                    if len(comp['SiO2_Liq']) != len(P_bar):
                        raise Warning("The length of your composition and pressure variables are different. Please correct this before running the code.")

                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Affinity': Affinity})
                else:
                    p = Process(target = findLiq, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Affinity': Affinity})

            ps.append(p)
            p.start()

        TIMEOUT = 60

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

    Affin = {}
    if type(comp) != dict or type(P_bar) == np.ndarray:
        Results = pd.DataFrame(data = np.zeros((len(T_Liq), 17)), columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
        for i in range(len(qs)):
            if len(qs[i])>0:
                if Affinity is True:
                    Res, Af, index = qs[i]
                    Affin[str(index)] = Af
                else:
                    Res, index = qs[i]

                Results.loc[index,:] = Res
    else:
        if Affinity is True:
            Results, Affin, index = qs[0]
        else:
            Results, index = qs[0]

    Res = comp_fix(Model = Model, comp = Results)
    if type(Res) == dict:
        Res = pd.DataFrame.from_dict(Res, orient = "index").T

    Res = Res[['T_Liq', 'liquidus_phase', 'fluid_saturated', 
       'SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq','MnO_Liq', 'MgO_Liq', 'CaO_Liq',
       'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 
       'Fe3Fet_Liq', 'Fe2O3', 'FeO', ]]
    
    Res = Res.drop(columns = ['Fe2O3',  'FeO']) # 'Cr2O3',

    if Affinity is True:
        # Affinity combined
        int_keys = list(map(int, Affin.keys()))

        # Determine the size of the resulting DataFrame
        max_index = max(int_keys)

        # Initialize an empty list to maintain column order
        all_columns = []
        seen_columns = set()

        for d in Affin:
            for col in Affin[d]:
                if col not in seen_columns:
                    all_columns.append(col)
                    seen_columns.add(col)
            
            Affin[d] = pd.Series(Affin[d])

        Af_Combined = pd.DataFrame(np.nan, index=range(max_index + 1), columns=all_columns)

        # Populate the result DataFrame
        for d in Affin:
            try:
                df = Affin[d].to_frame().transpose()
                Af_Combined.loc[int(d)] = df.iloc[0]
            except:
                continue
        
        return Res, Af_Combined
    else:
        return Res

def findCO2(q, index, *, Model = None, P_bar = None, T_initial_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    T_Liq = 0
    H2O_Melt = 0
    CO2_Melt = 0

    if "MELTS" in Model:
        #try:
        T_Liq, H2O_Melt, CO2_Melt = findCO2_MELTS(P_bar = P_bar, Model = Model, T_C = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([T_Liq, H2O_Melt, CO2_Melt, index])
        return
        #except:
        #    q.put([T_Liq, H2O_Melt, CO2_Melt, index])
        #    return

def findLiq(q, index,*, Model = None, P_bar = None, T_initial_C = None, comp = None, fO2_buffer = None, fO2_offset = None, Affinity = False):
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

    T_in = T_initial_C

    
    if "MELTS" in Model:
        try:
            if Affinity is True:
                Results, Affin = findLiq_MELTS(P_bar = P_bar, Model = Model, T_C_init = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, Affinity = Affinity)
            else:
                Results = findLiq_MELTS(P_bar = P_bar, Model = Model, T_C_init = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, Affinity = Affinity)

            if Affinity is True:
                q.put([Results, Affin, index])
            else:
                q.put([Results, index])
            return
        except:
            Results = {} #pd.DataFrame(data = np.zeros((1, 17)), columns = ['T_Liq', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'])
            Affin = {}
            if Affinity is True:
                q.put([Results, Affin, index])
            else:
                q.put([Results, index])
            return

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        #try:
        Results = MM.findLiq(P_bar = P_bar, T_C_init = T_initial_C, comp = comp)#, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index, T_in])
        return
        #except:
        #    q.put([T_Liq, H2O_Melt, index, T_in])
        #    return

def equilibrate(q, i,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    Res = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
    q.put([Res, i])
    return