import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.Plotting import *
from petthermotools.MELTS import *
# try:
#     from petthermotools.Holland import *
# except:
#     pass
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange
from pathlib import Path

def equilibrate_multi(cores = multiprocessing.cpu_count(), Model = "MELTSv1.0.2", bulk = None, T_C = None, P_bar = None, 
                      Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None,
                      timeout = None, copy_columns = None, Suppress = None):
    '''
    Runs single-step phase equilibrium calculations (isothermal and isobaric) for a batch of 
    compositions and/or P-T conditions in parallel using MELTS or MAGEMinCalc.
    
    This function handles input type conversion, broadcasting of scalar P-T inputs to 
    match DataFrame compositions, parallel process management, result stitching, and unit standardization.

    Parameters
    ----------
    cores : int, default multiprocessing.cpu_count()
        Number of parallel processes to use.
    Model : str, default "MELTSv1.0.2"
        The thermodynamic model to use (e.g., "MELTSv1.0.2", "pMELTS", or MAGEMin variant).
    bulk : dict or pd.DataFrame
        Input bulk composition(s). If a DataFrame, each row is a separate calculation.
        If a dict, it will be replicated for scalar P/T arrays.
    T_C : float or np.ndarray
        Temperature(s) in Celsius for the equilibrium step. Must match the number of compositions  
        if `bulk` is a DataFrame.
    P_bar : float or np.ndarray
        Pressure(s) in bars for the equilibrium step. Must match the number of compositions
        if 'bulk' is a DataFrame.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or np.ndarray, optional
        Overrides for initial liquid redox state and volatile contents.
    fO2_buffer : str, optional
        Oxygen fugacity buffer ("FMQ" or "NNO").
    fO2_offset : float, optional
        Offset in log units from the specified fO2 buffer.
    timeout : int, optional
        Timeout (in seconds) for each individual calculation process. Default is 20s.
    copy_columns : str or list of str, optional
        Column name(s) from the original `bulk` DataFrame to include in the final output 
        DataFrame (e.g., 'Sample_ID').
    Suppress : list of str, optional
        List of phases (e.g., ['rutile']) to exclude from the calculation (MELTS only).

    Returns
    -------
    Combined : pd.DataFrame
        A comprehensive DataFrame where each row is an equilibrium calculation result.
        It contains stitched results (Conditions, phase compositions, and phase properties)
        and, for MELTS, the affinity values (suffix `_affinity`).
    '''

    T_C = to_float(T_C)

    P_bar = to_float(P_bar)

    H2O_Liq = to_float(H2O_Liq)
    CO2_Liq = to_float(CO2_Liq)
    Fe3Fet_Liq = to_float(Fe3Fet_Liq)
    fO2_offset = to_float(fO2_offset)

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"

    if "MELTS" in Model:
        try:
            from meltsdynamic import MELTSdynamic
        except:
            Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

    comp = bulk.copy()

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq=CO2_Liq)

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

    if "MELTS" in Model:
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
                    p = Process(target = equilibrate, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_C': T_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Suppress': Suppress})

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
                        if first is True:
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
                    p = Process(target = equilibrate, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_C': T_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Suppress': Suppress})

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
                        if first is True:
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
                # try:
                Results, index = qs[0]
                Output = Results[0]
                Affinity = pd.DataFrame([Results[1]])
                Af_Combined = Affinity.copy()
            
                if "MELTS" in Model:
                    Output = stich(Output, Model = Model)

                Combined = Output['All'].copy()
                # except:
                #     return Output

        if copy_columns is not None:
            if type(copy_columns) == str:
                Combined.insert(0, copy_columns, comp[copy_columns])
            elif type(copy_columns) == list:
                j = 0
                for i in copy_columns:
                    Combined.insert(j, i, comp[i])
                    j = j + 1

        Af_Combined = Af_Combined.add_suffix('_affinity')
        Combined = pd.concat([Combined, Af_Combined], axis = 1)

        Output['Affinity'] = Af_Combined
        return Combined
    else:
        from juliacall import Main as jl, convert as jlconvert
        env_dir = Path.home() / ".petthermotools_julia_env"
        jl_env_path = env_dir.as_posix()

        jl.seval(f"""
            import Pkg
            Pkg.activate("{jl_env_path}")
            """)

        jl.seval("using MAGEMinCalc")

        comp['O'] = comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq'])

        if Model == "Weller2024":
            bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq']].astype(float).values
        else:
            bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq', 'H2O_Liq']].astype(float).values

        print(np.shape(bulk))
        bulk_jl = jl.seval("collect")(bulk)

        if type(T_C) == np.ndarray:
            T_C = jl.seval("collect")(T_C)
        if type(P_bar) == np.ndarray:
            P_kbar = jl.seval("collect")(P_bar/1000.0)
        else:
            P_kbar = P_bar/1000.0
        if type(fO2_offset) == np.ndarray:
            fO2_offset = jl.seval("collect")(fO2_offset)

        Output = jl.MAGEMinCalc.equilibrate(bulk = bulk_jl, P_kbar = P_kbar, T_C = T_C, fo2_buffer = fO2_buffer, fo2_offset = fO2_offset, Model = Model)
        Output = dict(Output)
        Combined = stich(Output, Model = Model)

        if copy_columns is not None:
            if type(copy_columns) == str:
                Combined['All'].insert(0, copy_columns, comp[copy_columns])
                # Af_Combined.insert(0, copy_columns, comp[copy_columns])
            elif type(copy_columns) == list:
                j = 0
                for i in copy_columns:
                    Combined['All'].insert(j, i, comp[i])
                    # Af_Combined.insert(j, i, comp[i])
                    j = j + 1

        return Combined['All']

def multi_equilibrate(q, index, *, Model = None, comp = None,
            T_C = None, P_bar = None, fO2_buffer = None, fO2_offset = None,
            Suppress = None, Suppress_except = None, trail = True):
    """
    NOT CURRENTLY IN USE
    
    Worker function to run a subset of equilibration models (MELTS or MAGEMin) in parallel.

    This function is intended to be run in a separate process. It takes a set of indices representing model runs,
    executes them using the appropriate model interface, and returns the results via a multiprocessing queue.

    Parameters
    ----------
    q : multiprocessing.Queue
        Output queue for sending back results.
    index : list of int
        Indices of the simulations to be run by this worker.
    Model : str
        The thermodynamic model ("MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", "pMELTS", or MAGEMin variant: "Green2025" or "Weller2024").
    comp : dict or pd.DataFrame
        Starting compositions. Either a single dictionary or a DataFrame with one row per simulation.
    Frac_solid : bool
        If True, removes solids at each step.
    Frac_fluid : bool
        If True, removes fluids at each step.
    T_C, T_path_C, T_start_C, T_end_C, dt_C : float or np.ndarray
        Temperature constraints or paths for each simulation.
    P_bar, P_path_bar, P_start_bar, P_end_bar, dp_bar : float or np.ndarray
        Pressure constraints or paths for each simulation.
    isenthalpic, isentropic, isochoric : bool
        Apply respective thermodynamic constraints.
    find_liquidus : bool
        If True, finds the liquidus temperature before starting.
    fO2_buffer : str
        Oxygen fugacity buffer ("FMQ" or "NNO").
    fO2_offset : float or array
        Offset from specified fO2 buffer in log units.
    fluid_sat : bool
        If True, terminates runs at fluid saturation.
    Crystallinity_limit : float
        Ends run when crystallinity exceeds this value.
    Suppress : list of str
        Phases to exclude from results.
    Suppress_except : bool
        If True, treat `Suppress` as a whitelist.
    trail : bool
        If True, include trailing properties from model output.

    Returns
    -------
    None
        The function returns results using `q.put()`:
        q.put([idx, results])
        where `idx` is the list of completed indices and `results` is a dictionary of output per run.
    """

    results = {}
    idx = []

    if len(index) > 15:
        index = index[:15]
    
    if "MELTS" in Model:
        from meltsdynamic import MELTSdynamic

        if Model is None or Model == "MELTSv1.0.2":
            melts = MELTSdynamic(1)
        elif Model == "pMELTS":
            melts = MELTSdynamic(2)
        elif Model == "MELTSv1.1.0":
            melts = MELTSdynamic(3)
        elif Model == "MELTSv1.2.0":
            melts = MELTSdynamic(4)
    else:
        from juliacall import Main as jl
        env_dir = Path.home() / ".petthermotools_julia_env"
        jl_env_path = env_dir.as_posix()

        jl.seval(f"""
            import Pkg
            Pkg.activate("{jl_env_path}")
            """)

        jl.seval("using MAGEMinCalc")

    for i in index:
        try:
            if "MELTS" in Model:
                if type(comp) == dict:
                    Results, tr = equilibrate_MELTS(Model = Model, comp = comp, 
                                            T_C = T_C[i], P_bar = P_bar[i], 
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                            Suppress = Suppress, Suppress_except=Suppress_except, trail = trail, melts = melts)
                else:
                    Results, tr = equilibrate_MELTS(Model = Model, comp = comp.loc[i].to_dict(), 
                                            T_C = T_C[i], P_bar = P_bar[i], 
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                            Suppress = Suppress, Suppress_except=Suppress_except, trail = trail, melts = melts)
            
                idx.append(i)

                results[f"Run {i}"] = Results

                if tr is False:
                    break
        except:
            idx.append(i)
            break
    
    if "MELTS" not in Model:
        jl.seval("using MAGEMinCalc")

        comp['O'] = comp['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq'])

        if Model == "Weller2024":
            bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq']].astype(float).values
        else:
            bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq', 'H2O_Liq']].astype(float).values

        print(np.shape(bulk))
        bulk_jl = jl.seval("collect")(bulk)

        if type(T_C) == np.ndarray:
            T_C = jl.seval("collect")(T_C)
        if type(P_bar) == np.ndarray:
            P_kbar = jl.seval("collect")(P_bar/1000.0)
        else:
            P_kbar = P_bar/1000.0
        if type(fO2_offset) == np.ndarray:
            fO2_offset = jl.seval("collect")(fO2_offset)

        Output = jl.MAGEMinCalc.equilibrate(bulk = bulk_jl, P_kbar = P_kbar, T_C = T_C, fo2_buffer = fO2_buffer, fo2_offset = fO2_offset, Model = Model)
        results = dict(Output)
        # results = stich(Output, Model = Model)

        idx = index

    q.put([idx, results])
    return


def findCO2_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
    '''
    Calculates the $\text{CO}_2$ saturation limit (solubility) in the melt at 
    given P-T conditions using a $\text{CO}_2$-incorporating thermodynamic model (MELTSv1.1.0 or MELTSv1.2.0).

    This function runs the calculation for multiple compositions or pressure points in parallel.

    Note: This function is currently only available for **MELTS** thermodynamic models 
    that incorporate $\text{CO}_2$ into the fluid (MELTSv1.1.0 or MELTSv1.2.0).

    Parameters
    ----------
    cores : int, optional
        Number of parallel processes to use. Defaults to all available CPU cores.
    Model : str
        The thermodynamic model to use (must be a MELTS version supporting $\text{CO}_2$).
    bulk : dict or pd.DataFrame
        Input bulk composition(s).
    T_initial_C : float or np.ndarray, optional
        Starting temperature(s) in Celsius for the calculation. If None, defaults to $1300\,^{\circ}\text{C}$.
    P_bar : float or np.ndarray
        Pressure(s) in bars at which to determine the $\text{CO}_2$ saturation.
    Fe3Fet_Liq : float or np.ndarray, optional
        Initial $\text{Fe}^{3+}/\text{Fe}_{\text{total}}$ ratio override.
    H2O_Liq : float or np.ndarray, optional
        Initial $\text{H}_2\text{O}$ content in the system (wt%).
    fO2_buffer : str, optional
        Oxygen fugacity buffer ("FMQ" or "NNO").
    fO2_offset : float or np.ndarray, optional
        Offset in log units from the specified $\text{f}\text{O}_2$ buffer.

    Returns
    -------
    T_Liq : np.ndarray
        Liquid temperature ($^{\circ}\text{C}$) at which the $\text{CO}_2$ saturation was determined.
    H2O : np.ndarray
        $\text{H}_2\text{O}$ content (wt%) in the liquid at $\text{CO}_2$ saturation.
    CO2 : np.ndarray
        $\text{CO}_2$ content (wt%) in the liquid at saturation.
    '''

    # try:
    #     from meltsdynamic import MELTSdynamic
    # except:
    #     Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"

    comp = bulk.copy()

    if Model is None:
        Model = "MELTSv1.0.2"

    if "MELTS" not in Model:
        return "the find CO2 function is only available for thermodynamic models incorporating CO2 into the fluid"
    # if Model == "Holland":
    #     import pyMAGEMINcalc as MM

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

def findLiq_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None, 
                  Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None, Affinity = False):
    '''
    Performs multiple liquidus (`findLiq`) calculations for a batch of 
    compositions and/or P-fO2-H2O conditions in parallel using the MELTS dynamic library.

    It manages the creation of parallel processes, distributes the workload, collects 
    results from the queue, and compiles them into a clean DataFrame.

    Note: Currently, liquidus calculation is only fully implemented for MELTS models.

    Parameters:
    ----------
    cores : int, optional
        Number of processes to run in parallel. If not provided, uses all available CPU cores.
    Model : str, optional
        The thermodynamic model to use (e.g., "MELTSv1.0.2", "pMELTS"). Default is "MELTSv1.0.2".
    bulk : pd.DataFrame or dict
        A DataFrame (preferred for multiple runs) or dictionary containing oxide compositions (wt%).
    T_initial_C : float or np.ndarray, optional
        Initial guess temperature in Celsius for the liquidus search. Default is $1300\,^{\circ}\text{C}$.
    P_bar : float or np.ndarray
        Pressure(s) in bars for the calculations. If an array, the length must match the number of compositions 
        if `bulk` is a DataFrame.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or np.ndarray, optional
        Overrides for initial liquid redox state and volatile contents (wt%).
    fO2_buffer : str, optional
        Oxygen fugacity buffer ("FMQ" or "NNO").
    fO2_offset : float or np.ndarray, optional
        Offset from the specified $\text{f}\text{O}_2$ buffer (in log units).
    Affinity : bool, optional, default False
        If True, returns additional DataFrame containing chemical **affinity** (driving force for crystallization) 
        data for all phases.

    Returns:
    ----------
    Res : pd.DataFrame
        DataFrame containing the calculated liquidus temperature ($T_{\text{Liq}}$), the liquidus phase, 
        fluid saturation status, and the normalized melt composition at the liquidus.
    Af_Combined : pd.DataFrame, optional
        Returned only if `Affinity=True`. Contains chemical affinity data for the liquidus phase.
    '''
    # try:
    #     from meltsdynamic import MELTSdynamic
    # except:
    #     Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

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

    if Model is None:
        Model = "MELTSv1.0.2"

    # if Model == "Holland":
    #     import pyMAGEMINcalc as MM

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
            
    if "MELTS" in Model:
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

        Res = comp_fix(Model = Model, comp = Results, keep_columns=True)
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
    else:
        # T_Liq = MM.findLiq_multi(P_bar = P_bar, T_initial_C = T_initial_C, comp = comp)
        return "find liquidus calculations are currently not available through the MAGEMin models. This is an issue I'm working to fix as soon as possible."

def findCO2(q, index, *, Model = None, P_bar = None, T_initial_C = None, comp = None, fO2_buffer = None, fO2_offset = None):
    '''
    Worker function to determine $\text{CO}_2$ saturation (solubility) 
    in a melt at specified P-T conditions using a $\text{CO}_2$-incorporating MELTS model.
    
    This function is executed in a separate process by `findCO2_multi`.

    Parameters:
    ----------
    q : multiprocessing.Queue
        Output queue for sending back results.
    index : int
        Index of the calculation (row number) in the master dataset for result indexing.
    Model : str, optional
        MELTS thermodynamic model version (must support $\text{CO}_2$, e.g., "MELTSv1.1.0").
    P_bar : float
        Pressure in bars.
    T_initial_C : float
        Starting temperature in Celsius.
    comp : dict
        Dictionary containing the oxide composition (wt%).
    fO2_buffer : str, optional
        Oxygen fugacity buffer.
    fO2_offset : float, optional
        Offset from the specified $\text{f}\text{O}_2$ buffer (log units).

    Returns:
    ----------
    None
        The function returns results via `q.put()`:
        `q.put([T_Liq, H2O_Melt, CO2_Melt, index])`
    '''
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
    Worker function that executes the liquidus temperature search (`findLiq_MELTS`) for a single 
    composition at a specified pressure.

    This function is executed in a separate process by `findLiq_multi`.

    Parameters:
    ----------
    q : multiprocessing.Queue
        Output queue for placing the results.
    index : int
        Index of the calculation (row number) in the master code for result indexing.
    Model : str, optional
        Thermodynamic model (e.g., "MELTSv1.0.2").
    P_bar : float
        Pressure of the calculation (bar).
    T_initial_C : float
        Initial 'guess' temperature in Celsius for the liquidus search.
    comp : dict
        Dictionary containing the oxide composition (wt%) for the calculation.
    fO2_buffer : str, optional
        Oxygen fugacity buffer.
    fO2_offset : float, optional
        Offset from the specified $\text{f}\text{O}_2$ buffer (log units).
    Affinity : bool, optional, default False
        If True, calculates and returns the chemical affinity of the liquidus phase.

    Returns:
    ----------
    None
        The function returns results via `q.put()`:
        - If `Affinity=True`: `q.put([Results, Affin, index])`
        - If `Affinity=False`: `q.put([Results, index])`
        where `Results` is a dictionary/series of liquidus properties and `Affin` is the affinity data.
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

def equilibrate(q, i,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None, Suppress = None):
    '''
    Worker function to execute a single-step phase equilibrium calculation (`equilibrate_MELTS`) 
    at specified P-T conditions.

    This function is executed in a separate process by `equilibrate_multi`.

    Parameters:
    ----------
    q : multiprocessing.Queue
        Output queue for sending back results.
    i : int
        Index of the calculation (row number) in the master code for result indexing.
    Model : str, optional
        MELTS thermodynamic model version.
    P_bar : float, optional
        Pressure in bars.
    T_C : float, optional
        Temperature in Celsius.
    comp : dict, optional
        Dictionary containing the oxide composition (wt%).
    fO2_buffer : str, optional
        Oxygen fugacity buffer.
    fO2_offset : float, optional
        Offset from the specified $\text{f}\text{O}_2$ buffer (log units).
    Suppress : list of str, optional
        List of phases to suppress (exclude) from the calculation.

    Returns:
    ----------
    None
        The function returns results via `q.put()`:
        `q.put([Res, i])`
        where `Res` is the raw output dictionary from the MELTS calculation.
    '''
    Res = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, Suppress = Suppress)
    q.put([Res, i])
    return