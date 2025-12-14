import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.Plotting import *
from petthermotools.Path import *
from petthermotools.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange
import random
import time


def make_grid_variables_from_phase_diagram(Results = None, x_col = "T_C", y_col = "P_bar"):
    """
    Grid all variables in a phaseDiagram_calc-style dataframe
    onto the same X–Y mesh.

    Parameters
    ----------
    df : pandas.DataFrame
        Output from ptt.phaseDiagram_calc
    x_col, y_col : str
        Grid axes (most likely T–P or P–T)

    Returns
    -------
    X, Y : np.ndarray
        Meshgrid arrays
    Z_dict : dict
        Dictionary of gridded variables for ALL other columns
    """
    if Results is None:
        print('Please provide results from the ptt.phaseDiagram_calc function.')
        
    df = Results.copy()
    # Unique grid axes
    x_vals = np.unique(df[x_col])
    y_vals = np.unique(df[y_col])

    X, Y = np.meshgrid(x_vals, y_vals)
    shape = X.shape

    # All columns except the grid axes
    z_cols = [c for c in df.columns if c not in (x_col, y_col)]

    Z_dict = {}
    for col in z_cols:
        try:
            Z_dict[col] = df[col].values.reshape(shape).T
        except ValueError:
            # Skip columns that cannot be reshaped (e.g. strings, metadata)
            continue

    return X, Y, Z_dict


def phaseDiagram_calc(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, 
                      T_min_C = None, T_max_C = None, T_num = None, 
                      P_min_bar = None, P_max_bar = None, P_num = None, 
                      Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None,
                      Fe3Fet_init = None, H2O_init = None, CO2_init = None,
                      fO2_buffer = None, fO2_offset = None, i_max = 15, grid = True, refine = None):
    """
    Calculates a Phase Diagram (P-T grid) for a given bulk composition using MELTS or MAGEMinCalc 
    via parallel processing.

    The function handles setting up the T-P grid, distributing calculations across multiple 
    cores, and iteratively re-attempting failed points up to a maximum number of attempts (`i_max`).

    Parameters:
    -----------
    cores : int, optional
        The number of CPU cores for parallel processing. Defaults to all available cores.
    Model : str, optional, default "MELTSv1.0.2"
        The name of the thermodynamic model to use (e.g., "MELTSv1.0.2", "pMELTS", "Green2025").
    bulk : dict or pandas.Series
        The bulk composition of the system (wt% oxides).
    T_C, P_bar : array-like, optional
        Explicit arrays of temperature (°C) and pressure (bar) values. If provided, they 
        override the min/max/num parameters.
    T_min_C, T_max_C, T_num : float, float, int, optional
        Minimum, maximum, and number of temperature points for the grid creation.
    P_min_bar, P_max_bar, P_num : float, float, int, optional
        Minimum, maximum, and number of pressure points for the grid creation.
    Fe3Fet_init, H2O_init, CO2_init : float, optional
        Initial $\text{Fe}^{3+}/\Sigma\text{Fe}$ ratio, $\text{H}_2\text{O}$ (wt%), and $\text{CO}_2$ (wt%) 
        content used to initialize the bulk system. These replace the deprecated `_Liq` arguments.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float, optional
        **DEPRECATED**. Use `Fe3Fet_init`, `H2O_init`, `CO2_init` instead.
    fO2_buffer : str, optional
        Oxygen buffer to constrain redox state (e.g., "FMQ", "NNO").
    fO2_offset : float, optional
        Offset to apply to the $\text{f}\text{O}_2$ buffer value (log units).
    i_max : int, optional, default 15
        The maximum number of times the parallel loop will re-attempt to calculate missing P-T points.
    grid : bool, optional, default True
        If True, a full 2D T-P meshgrid is created. If False, T\_C and P\_bar must be 1D arrays of the same length.
    refine : int, optional
        The number of times to run `phaseDiagram_refine` after the initial calculation to improve 
        phase boundary resolution.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the results for all successfully calculated P-T points, 
        including phase fractions, compositions, and thermodynamic properties.
    """

    ## make sure everything is a float
    T_C        = to_float(T_C)
    T_min_C   = to_float(T_min_C)
    T_max_C  = to_float(T_max_C)
    T_num   = to_float(T_num)

    P_bar      = to_float(P_bar)
    P_min_bar = to_float(P_min_bar)
    P_max_bar= to_float(P_max_bar)
    P_num  = to_float(P_num)

    Fe3Fet_init= to_float(Fe3Fet_init)
    Fe3Fet_Liq = to_float(Fe3Fet_Liq)
    H2O_init   = to_float(H2O_init)
    H2O_Liq    = to_float(H2O_Liq)
    CO2_init   = to_float(CO2_init)
    CO2_Liq    = to_float(CO2_Liq)
    fO2_offset = to_float(fO2_offset)

    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.0.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.0.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.0.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    comp = bulk.copy()

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"

    if cores is None:
        cores = multiprocessing.cpu_count()

    if Model is None:
        Model = "MELTSv1.0.2"

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq=CO2_init)

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
            if "MELTS" in Model:
                T_path_C = T_flat[i*A:(i+1)*A]
                P_path_bar = P_flat[i*A:(i+1)*A]
            else:
                T_path_C = np.array(subarrays_T[i])
                P_path_bar = np.array(subarrays_P[i])


            if "MELTS" in Model:
                if len(T_path_C) > 150:
                    T_path_C = T_path_C[:99]
                    P_path_bar = P_path_bar[:99]
            else:
                if len(T_path_C) > 500:
                    T_path_C = T_path_C[:499]
                    P_path_bar = P_path_bar[:499]

            if "MELTS" in Model:
                if j % 3 == 0:
                    ed = (5*np.random.random())
                    T_path_C = T_path_C[round(len(T_path_C)/ed):]
                    P_path_bar = P_path_bar[round(len(P_path_bar)/ed):]

            if j % 2 > 0:
                T_path_C = np.flip(T_path_C)
                P_path_bar = np.flip(P_path_bar)

            p = Process(target = path, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 
                                                               'T_path_C': T_path_C, 
                                                               'P_path_bar': P_path_bar, 
                                                               'fO2_buffer': fO2_buffer, 
                                                               'fO2_offset': fO2_offset,
                                                               'Suppress': ['rutile', 'tridymite']})

            ps.append(p)
            p.start()

        if "MELTS" in Model:
            TIMEOUT = 240
        else:
            TIMEOUT = 1200

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
                    Combined = Combined.dropna(subset = ['h_J'])

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
        A = np.zeros((np.shape(new_flat)[0], np.shape(Combined.values)[1]))
        A[:,:2] = new_flat

        B = Combined.values

        C = np.concatenate((A, B))

        Combined = pd.DataFrame(columns = list(Combined.keys()), data = C)

    Combined['T_C'] = np.round(Combined['T_C'], 2)
    Combined['P_bar'] = np.round(Combined['P_bar'], 2)
    Combined = Combined.sort_values(['T_C', 'P_bar'])
    Combined = Combined.reset_index(drop = True)
    Combined = Combined.dropna(subset = ['T_C'])

    if refine is not None:
        for i in range(refine):
            Combined = phaseDiagram_refine(Data = Combined, Model = Model, bulk = bulk, Fe3Fet_Liq=Fe3Fet_init,
                                           H2O_Liq=H2O_init, CO2_Liq=CO2_init, fO2_buffer=fO2_buffer, fO2_offset=fO2_offset, i_max = i_max)

    return Combined

def phaseDiagram_refine(Data = None, Model = None, bulk = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, 
                        fO2_buffer = None, fO2_offset = None, i_max = 15):
    """
    Refines an existing Phase Diagram calculation by identifying and calculating new P-T points 
    near inferred phase boundaries.

    The function first generates a finer grid of P and T points. It then uses the phase assemblage 
    information from the existing `Data` to interpolate values where phase assemblages are 
    homogeneous, and flags points where the phase assemblage changes (i.e., near phase boundaries) 
    for re-calculation using `phaseDiagram_calc`.

    Parameters:
    -----------
    Data : pandas.DataFrame
        The existing phase diagram calculation results (output from `phaseDiagram_calc`).
    Model : str, optional
        The thermodynamic model used for the initial calculation.
    bulk : dict or pandas.Series, optional
        The bulk composition of the system.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float, optional
        Initial system volatile contents and redox state used for the calculation.
    fO2_buffer : str, optional
        Oxygen buffer used for the calculation.
    fO2_offset : float, optional
        Offset to the $\text{f}\text{O}_2$ buffer value (log units).
    i_max : int, optional, default 15
        The maximum number of attempts for the sub-call to `phaseDiagram_calc` to resolve new points.

    Returns:
    --------
    pandas.DataFrame
        A consolidated DataFrame containing the original data, the newly interpolated data (for 
        homogeneous regions), and the newly calculated data (for phase boundary regions), 
        providing a refined P-T phase diagram.
    """
    
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
        return ','.join([col[5:] for col in Combined.loc[:, Combined.columns.str.contains('mass_g')].columns if row[col] > 0.0 and not pd.isna(row[col])])

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
                    A["T_C"] = flattened_T[i]
                    A["P_bar"] = flattened_P[i]
                else:
                    B = (matching_df.iloc[index_1,:] + matching_df.iloc[index_2,:] + matching_df.iloc[index_3,:] + matching_df.iloc[index_4,:])/4
                    B["T_C"] = flattened_T[i]
                    B["P_bar"] = flattened_P[i]
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

    
    T_C = np.array(T_C)
    P_bar = np.array(P_bar)

    if "MELTS" in Model:
        idx_add = np.where(matching_df['h_J'] == 0.0)[0]

        T_C = np.round(np.concatenate((T_C, matching_df.loc[idx_add,'T_C'].values)),2)
        P_bar = np.round(np.concatenate((P_bar, matching_df.loc[idx_add,'P_bar'].values)),2)

        matching_df = matching_df.loc[np.where(matching_df['h_J'] != 0.0)[0],:]
        matching_df = matching_df.reset_index(drop = True)

    New = phaseDiagram_calc(cores = multiprocessing.cpu_count(), 
                            Model = Model, bulk = bulk, T_C = T_C, P_bar = P_bar, Fe3Fet_init = Fe3Fet_Liq, H2O_init = H2O_Liq, CO2_init = CO2_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, i_max = i_max, grid = False)

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


# def phaseDiagram_eq(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, number_max = 50):
#     """
#     Calculates the phase diagram for a given bulk composition over a
#     range of temperatures and pressures.

#     Parameters:
#     -----------
#     cores : int or None, optional
#         The number of cores to use for multiprocessing, default is set to the number of logical processors available.
#     Model : str or None, optional
#         The thermodynamic model to use, either 'MELTSv1.0.2' or 'Holland'. Default is None.
#     bulk : dict, pandas.core.series.Series
#         The bulk composition to use. Can be a dictionary with the element names and mole
#         fractions, a pandas series.
#     T_C : array-like or None, optional
#         The temperature range (in degrees Celsius) to calculate the phase diagram over. If
#         T_min_C is given and T_C is None, this array is generated automatically. Default is None.
#     P_bar : array-like or None, optional
#         The pressure range (in bars) to calculate the phase diagram over. If P_min_bar is given
#         and P_bar is None, this array is generated automatically. Default is None.
#     T_min_C : float or None, optional
#         The minimum temperature (in degrees Celsius) to calculate the phase diagram over. Default is None.
#     T_max_C : float or None, optional
#         The maximum temperature (in degrees Celsius) to calculate the phase diagram over. Default is None.
#     T_num : int or None, optional
#         The number of temperature steps to calculate the phase diagram over. Default is None.
#     P_min_bar : float or None, optional
#         The minimum pressure (in bars) to calculate the phase diagram over. Default is None.
#     P_max_bar : float or None, optional
#         The maximum pressure (in bars) to calculate the phase diagram over. Default is None.
#     P_num : int or None, optional
#         The number of pressure steps to calculate the phase diagram over. Default is None.
#     Fe3Fet_Liq : float or None, optional
#         The Fe3+/Fet ratio of the liquid for the bulk composition. Default is None.
#     H2O_Liq : float or None, optional
#         The H2O content of the liquid for the bulk composition. Default is None.
#     fO2_buffer : str or None, optional
#         The oxygen buffer to use for MELTS. Default is None.
#     fO2_offset : float or None, optional
#         The oxygen buffer offset to use for MELTS. Default is None.
#     number_max : int, optional
#         The maximum number of calculations to perform in a single Holland model calculation.
#         If there are more than this number of calculations to perform, multiprocessing is used.
#         Default is 50.

#     Returns:
#     --------
#     pandas.core.frame.DataFrame
#         A dataframe containing the equilibrium phase assemblage(s) for the given bulk
#         composition over the range of temperatures and pressures specified.
#     """
#     comp = bulk.copy()

#     if cores is None:
#         cores = multiprocessing.cpu_count()

#     if Model is None:
#         Model = "MELTSv1.0.2"

#     if Model == "Holland":
#         import pyMAGEMINcalc as MM

#     # if comp is entered as a pandas series, it must first be converted to a dict
#     if type(comp) == pd.core.series.Series:
#         comp = comp.to_dict()

#     # ensure the bulk composition has the correct headers etc.
#     comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

#     if T_min_C is not None:
#         T_C = np.linspace(T_min_C, T_max_C, T_num)

#     if P_min_bar is not None:
#         P_bar = np.linspace(P_min_bar, P_max_bar, P_num)

#     T, P = np.meshgrid(T_C, P_bar)

#     T_flat = np.round(T.flatten(),2)
#     P_flat = np.round(P.flatten(),2)

#     if "MELTS" in Model:
#         A = len(P_flat)//cores
#         B = len(P_flat) % cores
#         if A > 0:
#             Group = np.zeros(A) + cores
#             if B > 0:
#                 Group = np.append(Group, B)
#         else:
#             Group = np.array([B])

#         qs = []
#         q = Queue()
#         c = 0
#         #Combined = pd.DataFrame(columns = ['T_C', 'P_bar'], data = np.zeros((1,2)))
#         for j in tqdm(range(len(Group))):
#             ps = []
#             for i in range(int(cores*j), int(cores*j + Group[j])):
#                 p = Process(target = equilibrate, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_C': T_flat[i], 'P_bar': P_flat[i], 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

#                 ps.append(p)
#                 p.start()

#             TIMEOUT = 60 #+ #0.5*len(T_flat)
#             start = time.time()
#             for p in ps:
#                 if TIMEOUT  - (time.time() - start) > 10:
#                     try:
#                         ret = q.get(timeout = TIMEOUT - (time.time() - start))
#                     except:
#                         ret = []

#                 else:
#                     try:
#                         ret = q.get(timeout = 10)
#                     except:
#                         ret = []

#                 qs.append(ret)

#             for p in ps:
#                 p.kill()

#         Results = {}
#         for i in range(len(qs)):
#             if len(qs[i]) > 0:
#                 Res, index = qs[i]
#                 Results['index = ' + str(index)] = Res

#         if "MELTS" in Model:
#             Results = stich(Results, multi = True, Model = Model)

#         for i in Results.keys():
#             if c == 0:
#                 if "MELTS" in Model:
#                     Combined = Results[i]['All'].copy()
#                 else:
#                     Combined = Results[i].copy()
#                 c = 1
#             else:
#                 if "MELTS" in Model:
#                     Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)
#                 else:
#                     try:
#                         Combined = pd.concat([Combined, Results[i]], axis = 0, ignore_index = True)
#                     except:
#                         continue

#     if Model == "Holland":
#         Combined = pd.DataFrame()
#         if len(T_flat) < number_max:
#             c = 0
#             for i in tqdm(range(len(T_flat))):
#                 try:
#                     Results = MM.equilibrate(Model = Model, P_bar = P_flat[i], T_C = T_flat[i], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
#                     Combined = pd.concat([Combined, Results], axis = 0, ignore_index = True)
#                 except:
#                     pass
#         else:
#             s = len(T_flat)//number_max
#             A = s//cores
#             B = s % cores

#             if A > 0:
#                 Group = np.zeros(A) + cores
#             if B > 0:
#                 if A > 0:
#                     Group = np.append(Group, B)
#                 else:
#                     Group = np.array([B])

#             qs = []
#             q = Queue()
#             for j in tqdm(range(len(Group))):
#                 ps = []
#                 for i in range(int(cores*j), int(cores*j + Group[j])):
#                     T2 = T_flat[i*number_max:i*number_max+number_max]
#                     P2 = P_flat[i*number_max:i*number_max+number_max]
#                     p = Process(target = equilibrate, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_C': T2, 'P_bar': P2, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

#                     ps.append(p)
#                     p.start()

#                 for p in ps:
#                     try:
#                         ret = q.get()
#                     except:
#                         ret = []

#                     qs.append(ret)

#                 for p in ps:
#                     p.kill()

#             for i in range(len(qs)):
#                 if len(qs[i]) > 0:
#                     Res, index = qs[i]
#                     Combined = pd.concat([Combined, Res], axis = 0, ignore_index = True)

#     if len(Combined['T_C']) < len(T_flat):
#         flat = np.round(np.array([T.flatten(), P.flatten()]).T,2)
#         Res_flat = np.round(np.array([Combined['T_C'].values, Combined['P_bar'].values]).T,2)
#         new_flat = flat[np.where((flat[:, None] == Res_flat).all(-1).any(-1) == False)]

#         # Combined['T_C'] = np.round(Combined['T_C'].values, 2)
#         # Combined['P_bar'] = np.round(Combined['P_bar'].values, 2)

#         for i in range(len(new_flat)):
#             df = pd.DataFrame(columns = ['T_C', 'P_bar'])
#             df.loc[0] = new_flat[i]

#             Combined = pd.concat([Combined, df], axis = 0, ignore_index = True)

#     return Combined

# def equilibrate(q, index,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):

#     if "MELTS" in Model:
#         Results = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
#         q.put([Results, index])
#         return

#     if "Holland" in Model:
#         import pyMAGEMINcalc as MM
#         #from tqdm.notebook import tqdm, trange
#         Combined = pd.DataFrame()
#         for i in range(len(T_C)):
#             try:
#                 Results = MM.equilibrate(Model = Model, P_bar = P_bar[i], T_C = T_C[i], comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
#                 Combined = pd.concat([Combined, Results], axis = 0, ignore_index = True)
#             except:
#                 pass

#         q.put([Combined, index])
#         return
