import numpy as np
import pandas as pd
from PetThermoTools.GenFuncs import *
from PetThermoTools.Plotting import *
from PetThermoTools.Path import *
from PetThermoTools.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange
import random
import time

def phaseDiagram_calc(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, i_max = 25, grid = True):
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

    if type(comp) == dict:
        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    if T_min_C is not None:
        T_C = np.linspace(T_min_C, T_max_C, T_num)

    if P_min_bar is not None:
        P_bar = np.linspace(P_min_bar, P_max_bar, P_num)

    if grid is True:
        P, T = np.meshgrid(P_bar, T_C)
    else:
        P, T = P_bar, T_C

    T_flat = np.round(T.flatten(),2)
    P_flat = np.round(P.flatten(),2)

    T_flat = np.flip(T_flat)
    P_flat = np.flip(P_flat)

    flat = np.array([T_flat, P_flat]).T

    A = len(P_flat)//cores + 1

    c = 0
    j = 0


    # subarray_length = len(T_flat) // cores

    # # Initialize an empty list to store subarrays
    # subarrays_T = []

    # Loop through the indices and create subarrays
    # for i in range(subarray_length):
    #     subarray_T = T_flat[i::subarray_length]
    #     subarrays_T.append(subarray_T)
        
    #     subarray_P = T_flat[i::subarray_length]
    #     subarrays_P.append(subarray_P)

    while len(T_flat)>1:
        if j > i_max:
            break

        # Initialize an empty list to store subarrays
        subarrays_T = []
        subarrays_P = []

        # Loop through the indices and create subarrays
        for i in range(cores):
            subarray_T = T_flat[i::cores]
            subarrays_T.append(subarray_T)
            
            subarray_P = P_flat[i::cores]
            subarrays_P.append(subarray_P)

        print('Attempt ' + str(j))

        qs = []
        q = Queue()
        ps = []
        j = j + 1

        for i in range(cores):
            T_path_C = np.array(subarrays_T[i])#T_flat[i*A:(i+1)*A]
            P_path_bar = np.array(subarrays_P[i])#P_flat[i*A:(i+1)*A]

            if len(T_path_C) > 150:
                T_path_C = T_path_C[:99]
                P_path_bar = P_path_bar[:99]

            if "MELTS" in Model:
                if j % 3 == 0:
                    ed = (5*np.random.random())
                    T_path_C = T_path_C[round(len(T_path_C)/ed):]
                    P_path_bar = P_path_bar[round(len(P_path_bar)/ed):]

            if j % 2 > 0:
                T_path_C = np.flip(T_path_C)
                P_path_bar = np.flip(P_path_bar)

            # if j > 5:
            #     com = list(zip(T, P))

            #     # Step 2: Randomize the order of the combined list
            #     random.shuffle(com)

            #     # Step 3: Separate the pairs back into two arrays
            #     T_randomized, P_randomized = zip(*com)

            #     T_path_C = np.array(T_randomized)
            #     P_path_bar = np.array(P_randomized)

            p = Process(target = path, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_path_C': T_path_C, 'P_path_bar': P_path_bar, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

            ps.append(p)
            p.start()

        TIMEOUT = 120
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
                    print('Timeout Reached. Calculation will continue in a new process.')
                    first = False
                try:
                    ret = q.get(timeout = 30)
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

        Results = stich(Results, multi = True, Model = Model)

        for i in Results.keys():
            if len(Results[i]['All']['T_C']) > 1:
                if c == 0:
                    Combined = Results[i]['All'].copy()
                    c = 1
                else:
                    Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)

            Combined['T_C'] = np.round(Combined['T_C'], 2)
            Combined['P_bar'] = np.round(Combined['P_bar'], 2)
            Combined = Combined.sort_values(['T_C', 'P_bar'])
            Combined = Combined.reset_index(drop = True)
            if "MELTS" in Model:
                Combined = Combined.dropna(subset = ['h'])

        Res_flat = np.array([Combined['T_C'], Combined['P_bar']]).T
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
    Combined = Combined.dropna(subset = ['T_C'])

    return Combined

def phaseDiagram_refine(Data = None, Model = None, bulk = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, i_max = 25):
    Combined = Data.copy()
    # find existing T,P data
    T_C = Combined['T_C'].unique()
    P_bar = Combined['P_bar'].unique()

    # calculate new T,P points
    T_C_new = np.round(np.linspace(np.nanmin(T_C), np.nanmax(T_C), 2*len(T_C)-1),2)
    P_bar_new = np.round(np.linspace(np.nanmin(P_bar), np.nanmax(P_bar), 2*len(P_bar)-1),2)

    # set up the new T,P grid
    T, P = np.meshgrid(T_C_new, P_bar_new)

    # extract the phase information from the current dataframe
    def combine_headers(row):
        return ','.join([col[5:] for col in Combined.loc[:, Combined.columns.str.contains('mass')].columns if row[col] > 0.0 or pd.isna(row[col])])

    # Apply the function to each row
    Combined['PhaseResults'] = Combined.apply(combine_headers, axis=1).tolist()

    ## find values that match T
    # Convert Combined['T_C'] to a set of rounded unique values for faster lookup
    unique_combined_T_C = set(np.round(Combined['T_C'].unique(), 2))
    unique_combined_P_bar = set(np.round(Combined['P_bar'].unique(), 2))

    # Flatten T and P and round their values
    flattened_T = np.round(T.flatten(), 2)
    flattened_P = np.round(P.flatten(), 2)

    # Check for matches and print them
    matching_values_T = [value for value in flattened_T if value in unique_combined_T_C]
    matching_values_P = [value for value in flattened_P if value in unique_combined_P_bar]

    matching_df = Data.copy()
    T_C = []
    P_bar = []
    A = None
    for i in range(len(flattened_T)):
        if flattened_T[i] in unique_combined_T_C and flattened_P[i] in unique_combined_P_bar:
            continue
        elif flattened_T[i] in unique_combined_T_C and flattened_P[i] not in unique_combined_P_bar:
            sorted_arr = np.sort(np.unique(flattened_P))
            idx = np.searchsorted(sorted_arr, flattened_P[i])

            below = sorted_arr[idx - 1]
            above = sorted_arr[idx + 1]

            if Combined['PhaseResults'].iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == below))[0][0]] == Combined['PhaseResults'].iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == above))[0][0]]:
                if A is None:
                    A = (matching_df.iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == below))[0][0],:].values + matching_df.iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == above))[0][0],:].values)/2
                    A[0] = flattened_T[i]
                    A[1] = flattened_P[i]
                else:

                    B = (matching_df.iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == below))[0][0],:].values + matching_df.iloc[np.where((Combined['T_C'].values == flattened_T[i]) & (Combined['P_bar'].values == above))[0][0],:].values)/2
                    B[0] = flattened_T[i]
                    B[1] = flattened_P[i]
                    A = np.vstack((A,B))

            else:
                T_C.append(flattened_T[i])
                P_bar.append(flattened_P[i])
        elif flattened_T[i] not in unique_combined_T_C and flattened_P[i] in unique_combined_P_bar:
            sorted_arr = np.sort(np.unique(flattened_T))
            idx = np.searchsorted(sorted_arr, flattened_T[i])

            below = sorted_arr[idx - 1]
            above = sorted_arr[idx + 1]

            if Combined['PhaseResults'].iloc[np.where((Combined['T_C'].values == below) & (Combined['P_bar'].values == flattened_P[i]))[0][0]] == Combined['PhaseResults'].iloc[np.where((Combined['T_C'].values == above) & (Combined['P_bar'].values == flattened_P[i]))[0][0]]:
                if A is None:
                    A = (matching_df.iloc[np.where((Combined['T_C'].values == below) & (Combined['P_bar'].values == flattened_P[i]))[0][0],:].values + matching_df.iloc[np.where((Combined['T_C'].values == above) & (Combined['P_bar'].values == flattened_P[i]))[0][0],:].values)/2
                    A[0] = flattened_T[i]
                    A[1] = flattened_P[i]
                else:

                    B = (matching_df.iloc[np.where((Combined['T_C'].values == below) & (Combined['P_bar'].values == flattened_P[i]))[0][0],:].values + matching_df.iloc[np.where((Combined['T_C'].values == above) & (Combined['P_bar'].values == flattened_P[i]))[0][0],:].values)/2
                    B[0] = flattened_T[i]
                    B[1] = flattened_P[i]
                    A = np.vstack((A,B))

            else:
                T_C.append(flattened_T[i])
                P_bar.append(flattened_P[i])
        else:
            sorted_arr = np.sort(np.unique(flattened_T))
            idx = np.searchsorted(sorted_arr, flattened_T[i])

            below_T = sorted_arr[idx - 1]
            above_T = sorted_arr[idx + 1]

            sorted_arr = np.sort(np.unique(flattened_P))
            idx = np.searchsorted(sorted_arr, flattened_P[i])

            below_P = sorted_arr[idx - 1]
            above_P = sorted_arr[idx + 1]

            index_1 = np.where((Combined['T_C'].values == below_T) & (Combined['P_bar'].values == below_P))[0][0]
            index_2 = np.where((Combined['T_C'].values == above_T) & (Combined['P_bar'].values == below_P))[0][0]
            index_3 = np.where((Combined['T_C'].values == below_T) & (Combined['P_bar'].values == above_P))[0][0]
            index_4 = np.where((Combined['T_C'].values == above_T) & (Combined['P_bar'].values == above_P))[0][0]

            if Combined['PhaseResults'].loc[index_1] == Combined['PhaseResults'].loc[index_2] == Combined['PhaseResults'].loc[index_3] == Combined['PhaseResults'].loc[index_4]:
                if A is None:
                    A = (matching_df.iloc[index_1,:] + matching_df.iloc[index_2,:] + matching_df.iloc[index_3,:] + matching_df.iloc[index_4,:])/4
                    A[0] = flattened_T[i]
                    A[1] = flattened_P[i]
                else:
                    B = (matching_df.iloc[index_1,:] + matching_df.iloc[index_2,:] + matching_df.iloc[index_3,:] + matching_df.iloc[index_4,:])/4
                    B[0] = flattened_T[i]
                    B[1] = flattened_P[i]
                    A = np.vstack((A,B))

            else:
                T_C.append(flattened_T[i])
                P_bar.append(flattened_P[i])

    B = matching_df.values

    C = np.concatenate((A, B))

    matching_df = pd.DataFrame(columns = list(matching_df.keys()), data = C)

    matching_df['T_C'] = np.round(matching_df['T_C'], 2)
    matching_df['P_bar'] = np.round(matching_df['P_bar'], 2)
    matching_df = matching_df.sort_values(['T_C', 'P_bar'])
    matching_df = matching_df.reset_index(drop = True)
    matching_df = matching_df.dropna(subset = ['T_C'])

    idx_add = np.where(matching_df['h'] == 0.0)[0]
    
    T_C = np.array(T_C)
    P_bar = np.array(P_bar)

    T_C = np.round(np.concatenate((T_C, matching_df.loc[idx_add,'T_C'].values)),2)
    P_bar = np.round(np.concatenate((P_bar, matching_df.loc[idx_add,'P_bar'].values)),2)

    matching_df = matching_df.loc[np.where(matching_df['h'] != 0.0)[0],:]
    matching_df = matching_df.reset_index(drop = True)

    New = phaseDiagram_calc(cores = 1, Model = Model, bulk = bulk, T_C = T_C, P_bar = P_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, i_max = i_max, grid = False)

    # A = New.values
    # B = matching_df.values

    # C = np.concatenate((A,B))

    out = pd.concat([New,matching_df])

    out['T_C'] = np.round(out['T_C'], 2)
    out['P_bar'] = np.round(out['P_bar'], 2)
    out = out.sort_values(['T_C', 'P_bar'])
    out = out.reset_index(drop = True)
    out = out.dropna(subset = ['T_C'])

    return out

def tidy_phaseDiagram(Data = None, Model = None, bulk = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, i_max = 25):
    Combined = Data.copy()
    T_C = Combined.loc[np.where(Combined['h'] != 0.0)[0], 'T_C'].values
    P_bar = Combined.loc[np.where(Combined['h'] != 0.0)[0], 'P_bar'].values

    Combined = Combined.loc[np.where(Combined['h'] != 0.0)[0],:]
    Combined = Combined.reset_index(drop = True)

    New = phaseDiagram_calc(cores = 1, Model = Model, bulk = bulk, T_C = T_C, P_bar = P_bar, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, i_max = i_max, grid = False)

    out = pd.concat([New,Combined])

    out['T_C'] = np.round(out['T_C'], 2)
    out['P_bar'] = np.round(out['P_bar'], 2)
    out = out.sort_values(['T_C', 'P_bar'])
    out = out.reset_index(drop = True)
    out = out.dropna(subset = ['T_C'])

    return out


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
        #from tqdm.notebook import tqdm, trange
        Combined = pd.DataFrame()
        for i in range(len(T_C)):
            try:
                Results = MM.equilibrate(Model = Model, P_bar = P_bar[i], T_C = T_C[i], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
                Combined = pd.concat([Combined, Results], axis = 0, ignore_index = True)
            except:
                pass

        q.put([Combined, index])
        return
