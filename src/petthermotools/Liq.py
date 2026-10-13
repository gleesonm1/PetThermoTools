import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.GenFuncs import _ensure_julia_ready
from petthermotools.GenFuncs import _ensure_julia_main_ready
from petthermotools.GenFuncs import _quiet_engine_output
from petthermotools.GenFuncs import _permanently_quiet_engine_output
from petthermotools.core_config import MAX_WORKERS, resolve_cores
from petthermotools.mp_utils import run_tasks, TaskStatus
import contextlib
from petthermotools.Plotting import *
from petthermotools.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.auto import tqdm, trange
from pathlib import Path
import time
import warnings
import platform
from typing import NamedTuple

def supCalc(Model = "MELTSv1.0.2", bulk = None, phase = None, T_C = None, P_bar = None,
             Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None,
             melts = None, quiet_engine = True):
    '''
    A wrapper function to run a extrac the thermodynamic properties of a user specified
    phase in alphaMELTS for Python for a given composition and specified P-T-fO2 conditions.

    The function first standardizes the composition using `comp_fix` and then passes 
    the parameters to the underlying MELTS calculation function (`supCalc_MELTS`).

    Parameters:
    ----------
    Model : str, default "MELTSv1.0.2"
        Thermodynamic model to use (e.g., "MELTSv1.0.2", "pMELTS").
    bulk : dict or pd.Series
        Bulk composition of the starting material (oxides in wt%).
    phase : str, optional
        Phase name (e.g., 'olivine1') to monitor or target in the calculation.
    T_C : float
        Temperature in Celsius.
    P_bar : float
        Pressure in bars.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float, optional
        Overrides for initial redox state and volatile content in the liquid.
    fO2_buffer : str, optional
        Oxygen fugacity buffer (e.g., "FMQ", "NNO").
    fO2_offset : float, optional
        Offset in log units from the specified fO2 buffer.
    melts : MELTSdynamic object, optional
        An existing MELTSdynamic instance, reduces the overhead associated with initiating new MELTS objects.
    quiet_engine : bool, default True
        If True, redirects the MELTS engine's stdout/stderr to the null device for the duration of
        the calculation, which can give a speed-up by avoiding large volumes of text being written
        to the console/notebook. Set to False to inspect the raw engine output for debugging purposes.

    Returns:
    ----------
    Results : dict
        The raw results dictionary from the underlying MELTS simulation.
    '''

    comp = bulk.copy()

    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    comp = comp_fix(Model = Model, comp = comp, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, Fe3Fet_Liq = Fe3Fet_Liq)

    with (_quiet_engine_output() if quiet_engine else contextlib.nullcontext()):
        Results = supCalc_MELTS(Model = Model, comp = comp, phase = phase, T_C = T_C, P_bar = P_bar,
                 fO2_buffer = fO2_buffer, fO2_offset = fO2_offset,
                 melts = melts)

    return Results

def equilibrate_multi(cores = MAX_WORKERS, Model = "MELTSv1.0.2", bulk = None, T_C = None, P_bar = None,
                      Fe3Fet_init = None, H2O_init = None, CO2_init = None, 
                      Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None,
                      timeout = 20, copy_columns = None, Suppress = None, Suppress_except = None, quiet_engine = True):
    '''
    Runs single-step phase equilibrium calculations (isothermal and isobaric) for a batch of 
    compositions and/or P-T conditions in parallel using MELTS or MAGEMinCalc.
    
    This function handles input type conversion, broadcasting of scalar P-T inputs to 
    match DataFrame compositions, parallel process management, result stitching, and unit standardization.

    Parameters
    ----------
    cores : int, default petthermotools.core_config.MAX_WORKERS
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
        List of phases (e.g., ['rutile']) to exclude from the calculation. Works for both
        MELTS and MAGEMin (wired through MAGEMin_C's `remove_phases` for a MAGEMin Model).
    Suppress_except : bool, optional
        MELTS only -- treat `Suppress` as a whitelist instead of a blacklist. Has no MAGEMin
        equivalent (MAGEMin's suppress is always a blacklist); raises ValueError if set for a
        MAGEMin Model rather than being silently ignored.
    quiet_engine : bool, optional
        If True (default), the MELTS engine's own progress text is suppressed for speed. Set to
        False to see its raw per-step output for debugging.

    Returns
    -------
    Combined : pd.DataFrame
        A comprehensive DataFrame where each row is an equilibrium calculation result.
        It contains stitched results (Conditions, phase compositions, and phase properties)
        and, for MELTS, the affinity values (suffix `_affinity`).
    '''
    cores = resolve_cores(cores)

    T_C = to_float(T_C)
    T_C = check_array(T_C)

    P_bar = to_float(P_bar)
    P_bar = check_array(P_bar)

    H2O_Liq = to_float(H2O_Liq)
    CO2_Liq = to_float(CO2_Liq)
    Fe3Fet_Liq = to_float(Fe3Fet_Liq)
    fO2_offset = to_float(fO2_offset)

    H2O_Liq = check_array(H2O_Liq)
    CO2_Liq = check_array(CO2_Liq)
    Fe3Fet_Liq = check_array(Fe3Fet_Liq)
    fO2_offset = check_array(fO2_offset)

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

    ### NGIBBS BRANCH ###
    _ng_base = Model[:-6] if Model.endswith('NoProp') else Model
    if Model.startswith('n'):
        from petthermotools.nGibbs_bridge import nGibbsAPI, _nGibbs_models
    if Model.startswith('n') and _ng_base in _nGibbs_models:
        output = nGibbsAPI(Model = Model, comp = bulk, T_C = T_C, P_bar = P_bar,
                      Fe3Fet_init = Fe3Fet_init, H2O_init = H2O_init, CO2_init = CO2_init, 
                      Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset,
                      Suppress = None, Suppress_except = None)
        if copy_columns is not None: # Copy Columns as requested
            for df_name, df in output.items():
                for i, col in enumerate(copy_columns):
                    if col in bulk.columns:
                        df.insert(i, col, bulk[col].values)
        return output

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                if fO2_buffer != "IW":
                    raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"


    comp = bulk.copy()

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq=CO2_init)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            warnings.warn("If not performing subsolidus or melting calculations adding small amounts of H2O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H2O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    
    if type(fO2_offset) == np.ndarray:
        if "MELTS" in Model:
            length = len(fO2_offset)
        else:
            raise ValueError("Currently the code only allows a single value of fO2_offset to be supplied to MAGEMin. Please supply a single value and/or use the Fe3Fet_init kwarg")
    elif type(comp) == pd.core.frame.DataFrame:
        length = len(comp['SiO2_Liq'])
    elif type(P_bar) == np.ndarray:
        length = len(P_bar)
    elif type(T_C) == np.ndarray:
        length = len(T_C)
    else:
        length = 1

    if length > 1:
        if type(P_bar) != np.ndarray:
            P_bar = np.zeros(length) + P_bar

        if type(T_C) != np.ndarray:
            T_C = np.zeros(length) + T_C

        if fO2_offset is not None:
            if type(fO2_offset) != np.ndarray:
                if "MELTS" in Model:
                    fO2_offset = np.zeros(length) + fO2_offset
        else:
            if "MELTS" in Model:
                fO2_offset = [None]*length
            else:
                # raise Warning("Currently the code only allows a single value of fO2_offset to be supplied to MAGEMin. Please supply a single value and/or use the Fe3Fet_init kwarg")
                fO2_offset = None

    if length == 1 and "MELTS" not in Model:    
        if type(P_bar) != np.ndarray:
            P_bar = np.zeros(length) + P_bar

        if type(T_C) != np.ndarray:
            T_C = np.zeros(length) + T_C

    
    if length > 1:
        if type(comp) == pd.core.frame.DataFrame:
            if "MELTS" in Model:
                if len(comp['SiO2_Liq']) != len(T_C) or len(comp['SiO2_Liq']) != len(P_bar) or len(comp['SiO2_Liq']) != len(fO2_offset):
                    raise ValueError("Inconsistent length between loaded variables")
            else:
                if len(comp['SiO2_Liq']) != len(T_C) or len(comp['SiO2_Liq']) != len(P_bar):
                    raise ValueError("Inconsistent length between loaded variables")
        else:
            if "MELTS" in Model:
                if len(P_bar) != len(T_C) or len(P_bar) != len(fO2_offset):
                    raise ValueError("Inconsistent length between loaded variables")
            else:
                if len(P_bar) != len(T_C):
                    raise ValueError("Inconsistent length between loaded variables")


    if "MELTS" in Model:
        if length == 1:
            q = Queue()
            ps = []
            p = Process(target = equilibrate, args = (q, 0), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_C': T_C,
                                                                       'comp': comp, 'fO2_buffer': fO2_buffer,
                                                                       'fO2_offset': fO2_offset, 'Suppress': Suppress,
                                                                       'Suppress_except': Suppress_except, 'quiet_engine': quiet_engine})
            ps.append(p)
            p.start()

            if timeout is None:
                TIMEOUT = 20
            else:
                TIMEOUT = timeout
            
            start = time.time()
            for p in ps:
                try:
                    ret = q.get(timeout = TIMEOUT)
                except:
                    ret = []

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

            if len(ret)>1:
                Results, index = ret
                print(Results)
                Output = Results[0]
                Affinity = pd.DataFrame([Results[1]])
                Af_Combined = Affinity.copy()
            
                if "MELTS" in Model:
                    Output = stich(Output, Model = Model)

                Output['Affinity'] = Af_Combined

            return Output
        else:
            # if type(comp) == dict:

            timeout_main = timeout

            index_in = np.arange(int(length))
            combined_results = {}
            combined_affinity = {}
            index_out = np.array([], dtype=int)

            time_track = 0
            while len(index_out) < len(index_in):
                index = np.setdiff1d(index_in, index_out)
                groups = np.array_split(index, cores)
                non_empty_groups = [g for g in groups if g.size > 0]
                groups = non_empty_groups

                timeout = timeout_main*len(groups[0])

                processes = []
                Start = time.time()
                for i in range(len(groups)):
                    q = Queue()
                    p = Process(target = multi_equilibrate, args = (q, groups[i]),
                                    kwargs = {'Model': Model, 'comp': comp,
                                            'T_C': T_C,
                                            'P_bar': P_bar,
                                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset,
                                            'Suppress': Suppress, 'Suppress_except': Suppress_except,
                                            'quiet_engine': quiet_engine})
                    
                    processes.append((p, q, groups[i]))
                    p.start()

                for p, q, group in processes:
                    try:
                        if time.time() - Start < timeout + 5:
                            res = q.get(timeout=timeout)
                        else:
                            res = q.get(timeout=10)
                    except:
                        res = []
                        idx_chunks = np.array([group[0]], dtype = int)

                    p.join(timeout = 2)
                    p.terminate()

                    if len(res) > 0.0:
                        idx_chunks, results, affinity = res
                        combined_results.update(results)
                        combined_affinity.update(affinity)

                    index_out = np.hstack((index_out, idx_chunks))

                print(f"Completed {100*len(index_out)/len(index_in)} %")
                if time_track > 0:
                    if saved == len(index_out)/len(index_in):
                        print('Likely at least 1 calculation did not complete')
                        break

                time_track = time_track + 1
                saved = len(index_out)/len(index_in)
            
            Final_Combined = {}

            # identify unique dataframes across all calculations
            all_dict_keys = set()
            for run_data in combined_results.values():
                all_dict_keys.update(run_data.keys())

            #Consolidate data into combined DataFrams
            for key in all_dict_keys:
                for it in combined_results:
                    if key in combined_results[it].keys():
                        Final_Combined[key] = pd.DataFrame(np.nan, index = np.arange(length), columns = combined_results[it][key].columns)

                        # Final_Combined[key][["liquidus_phase","fluid_saturated"]] = (
                        #     Final_Combined[key][["liquidus_phase","fluid_saturated"]].astype("string")
                        # )

                for run_name, run_data in combined_results.items():
                    if key in run_data:
                        run_idx = int(run_name.split(' ')[1])
                        data = run_data[key]

                        Final_Combined[key].loc[run_idx,data.columns] = data.values[0]

            Output = stich(Final_Combined, Model = Model)

            return Output
    else:
        # Loads MAGEMinCalc onto the main embedded Julia session (and sets
        # PYTHON_JULIACALL_THREADS before juliacall's first import, so
        # multi_point_minimization's internal threading actually runs) if
        # not already done -- previously this call was commented out, so
        # this branch only worked if the user had separately called
        # activate_petthermotools_env() first; see _ensure_julia_main_ready's
        # own docstring for why this is the lightweight, no-Distributed path.
        _ensure_julia_main_ready()
        from juliacall import Main as jl, convert as jlconvert

        if Suppress_except is not None:
            raise ValueError("Suppress_except has no MAGEMin equivalent -- MAGEMinCalc's "
                              "suppress (wired through MAGEMin_C's remove_phases) is always "
                              "a blacklist of phases to remove, not a whitelist. Use Suppress "
                              "alone for a MAGEMin Model.")

        comp['O'] = comp['Fe3Fet_Liq']*(((159.69/2)/71.844)*comp['FeOt_Liq'] - comp['FeOt_Liq'])

        if type(comp) == dict:
            if Model == "Weller2024":
                bulk = np.array([[comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['O'], comp['Cr2O3_Liq']]], dtype=float)
            else:
                bulk = np.array([[comp['SiO2_Liq'], comp['Al2O3_Liq'], comp['CaO_Liq'], comp['MgO_Liq'], comp['FeOt_Liq'], comp['K2O_Liq'], comp['Na2O_Liq'], comp['TiO2_Liq'], comp['O'], comp['Cr2O3_Liq'], comp['H2O_Liq']]], dtype = float)
            
            # If a length is specified and is greater than 1, copy the row 'length' times
            if length > 1:
                bulk = np.tile(bulk, (length, 1))
        else:
            if Model == "Weller2024":
                bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq']].astype(float).values
            else:
                bulk = comp[['SiO2_Liq', 'Al2O3_Liq', 'CaO_Liq', 'MgO_Liq', 'FeOt_Liq', 'K2O_Liq', 'Na2O_Liq', 'TiO2_Liq', 'O', 'Cr2O3_Liq', 'H2O_Liq']].astype(float).values

        # print(np.shape(bulk))
        bulk_jl = jl.seval("collect")(bulk)

        if type(T_C) == np.ndarray:
            T_C = jl.seval("collect")(T_C)
        if type(P_bar) == np.ndarray:
            P_kbar = jl.seval("collect")(P_bar/1000.0)
        else:
            P_kbar = P_bar/1000.0
        if type(fO2_offset) == np.ndarray:
            fO2_offset = jl.seval("collect")(fO2_offset)

        # A plain Python list arrives in juliacall as a lazy PyList{Any},
        # which does NOT satisfy the Julia Union{Nothing, Vector{String}}
        # type -- confirmed directly (TypeError at the call below without
        # this). jl.seval("Vector{String}") is the actual Julia type/
        # constructor; calling it forces element-wise str -> String
        # conversion, not just a Vector{Any} wrapper.
        suppress_jl = jl.seval("Vector{String}")(Suppress) if Suppress is not None else None

        Output = jl.MAGEMinCalc.equilibrate(bulk = bulk_jl, P_kbar = P_kbar, T_C = T_C, fo2_buffer = fO2_buffer, fo2_offset = fO2_offset, Model = Model, suppress = suppress_jl)
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

        # def combine_headers(row):
        #     return ','.join([col[7:] for col in Combined['All'].loc[:, Combined['All'].columns.str.contains('mass_g_')].columns if row[col] > 0.0 and not pd.isna(row[col])])

        # Combined['All'].insert(0, 'Phase List', Combined['All'].apply(combine_headers, axis=1).tolist())

        return Combined

def multi_equilibrate(q, index, *, Model = None, comp = None,
            T_C = None, P_bar = None, fO2_buffer = None, fO2_offset = None,
            Suppress = None, Suppress_except = None, trail = True, quiet_engine = True):
    """
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

    # print(index)
    results = {}
    affinity = {}
    idx = []

    if len(index) > 150:
        index = index[:150]
    
    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's docstring for why
        # the old `with _quiet_engine_output():` here is intentionally gone, not just unused --
        # it let the engine's own buffered banner/status text leak through after the block
        # restored the real fd 1/2, confirmed directly (see MAINTENANCE_NOTES.md). MUST run
        # before MELTSdynamic() is constructed below -- that construction itself is what
        # prints the "alphaMELTS for Python version..." banner, confirmed directly (a version
        # of this fix that called this AFTER construction still leaked the banner).

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
    for i in index:
        melts = melts.addNodeAfter()
        try:
            if "MELTS" in Model:
                if type(comp) == dict:
                    Results, Affinity, tr = equilibrate_MELTS(Model = Model, comp = comp,
                                            T_C = T_C[i], P_bar = P_bar[i],
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i],
                                            Suppress = Suppress, Suppress_except=Suppress_except, trail = trail, melts = melts)
                else:
                    Results, Affinity, tr = equilibrate_MELTS(Model = Model, comp = comp.loc[i].to_dict(),
                                            T_C = T_C[i], P_bar = P_bar[i],
                                            fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i],
                                            Suppress = Suppress, Suppress_except=Suppress_except, trail = trail, melts = melts)

                idx.append(i)

                results[f"Run {i}"] = Results
                affinity[f"Run {i}"] = Affinity

                if tr is False:
                    break
        except:
            idx.append(i)
            break

    q.put([idx, results, affinity])
    return


# TODO(phase 1/2): findCO2_multi and its worker findCO2 (below) rely on the older per-worker, non-batched
# structure, are slow, and have no caller in any maintained notebook, test or doc. Revisit once the
# supervisor rewrite (Phase 2) and backend registry (Phase 1) land, rather than updating them in place now.
def findCO2_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, quiet_engine = True):
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
        Starting temperature(s) in Celsius for the calculation. If None, defaults to $1300\\,^{\\circ}\text{C}$.
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
    quiet_engine : bool, optional
        If True (default), the MELTS engine's own progress text is suppressed for speed. Set to
        False to see its raw per-step output for debugging.

    Returns
    -------
    T_Liq_C : np.ndarray
        Liquid temperature ($^{\\circ}\text{C}$) at which the $\text{CO}_2$ saturation was determined.
    H2O : np.ndarray
        $\text{H}_2\text{O}$ content (wt%) in the liquid at $\text{CO}_2$ saturation.
    CO2 : np.ndarray
        $\text{CO}_2$ content (wt%) in the liquid at saturation.
    '''


    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                if fO2_buffer != "IW":
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

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq)

    if type(comp) != dict or type(P_bar) == np.ndarray:
        if type(comp) != dict:
            T_Liq = np.zeros(len(comp['SiO2_Liq'].values))
            H2O_melt = np.zeros(len(comp['SiO2_Liq'].values))
            CO2_melt = np.zeros(len(comp['SiO2_Liq'].values))
            index = np.zeros(len(comp['SiO2_Liq'].values)) - 1
        else:
            T_Liq = np.zeros(len(P_bar))
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

    cores = resolve_cores(cores)

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
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'quiet_engine': quiet_engine})
                else:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'quiet_engine': quiet_engine})
            else:
                if type(P_bar) == np.ndarray:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar[i], 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'quiet_engine': quiet_engine})
                else:
                    p = Process(target = findCO2, args = (q, i), kwargs = {'Model': Model, 'P_bar': P_bar, 'T_initial_C': T_initial_C[i], 'comp': comp.loc[i].to_dict(), 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'quiet_engine': quiet_engine})

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
                T_Liq[i], H2O_melt[i], CO2_melt[i], index[i] = qs[i]

        T_Liq_C = np.zeros(len(T_Liq))
        H2O = np.zeros(len(T_Liq))
        CO2 = np.zeros(len(T_Liq))

        for i in range(len(index)):
            if len(CO2[index == i]) > 0:
                T_Liq_C[i] = T_Liq[index == i]
                H2O[i] = H2O_melt[index == i]
                CO2[i] = CO2_melt[index == i]
    else:
        T_Liq_C, H2O, CO2, index = qs[0]

    return T_Liq_C, H2O, CO2

# TODO(phase 1/2): worker for findCO2_multi; see the note above that function.
def findCO2(q, index, *, Model = None, P_bar = None, T_initial_C = None, comp = None, fO2_buffer = None, fO2_offset = None, quiet_engine = True):
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
        `q.put([T_Liq_C, H2O_Melt, CO2_Melt, index])`
    '''
    T_Liq_C = 0
    H2O_Melt = 0
    CO2_Melt = 0

    if "MELTS" in Model:
        #try:
        if quiet_engine:
            _permanently_quiet_engine_output()   # one-way; see that function's own docstring
        T_Liq_C, H2O_Melt, CO2_Melt = findCO2_MELTS(P_bar = P_bar, Model = Model, T_C = T_initial_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([T_Liq_C, H2O_Melt, CO2_Melt, index])
        return
        #except:
        #    q.put([T_Liq_C, H2O_Melt, CO2_Melt, index])
        #    return

class _RowTask(NamedTuple):
    """
    One row's own data for findLiq_multi, carried through
    mp_utils.run_tasks's `indices` (pickled ONCE per task, when it goes
    onto the task queue) instead of through `task_kwargs` (a dict shared
    by -- and re-pickled/re-transmitted to -- EVERY worker spawn).

    This replaced passing the FULL batch's P_bar/T_initial_C/comp arrays
    through task_kwargs, which worked fine at small N but was directly
    confirmed (via instrumented p.start() timing, isolated from any other
    load) to make every single spawn progressively, severely slower as N
    grows -- ~0.01s/spawn at N=40, ~0.05-0.2s at N=300, ~1.1-1.5s at
    N=1000 -- because `max_tasks_per_worker`'s default of 1 means EVERY
    task pays that cost once. At N=1000 that alone adds roughly
    1000 * 1.2s / cores of pure spawn overhead, on top of the actual MELTS
    computation, and compounds with a separate mp_utils-level race
    (see mp_utils.run_tasks's drain-loop comment) whose own cost scales
    with how expensive a wasted spawn is -- together these fully explained
    a batch that should take minutes instead projecting multiple hours.
    Giving each task only its own row fixes this at the root: a `_RowTask`
    is a few small scalars and one composition dict, not an O(n)-sized
    structure, so spawn cost no longer depends on how many rows are in the
    whole batch.
    """
    row_id: int
    P_bar: float
    T_initial_C: float
    comp: dict

    def __hash__(self):
        # mp_utils.run_tasks uses `indices` entries as dict keys / set
        # members (attempts, outstanding, results, tracked_by_pid lookups),
        # which requires hashability -- a plain NamedTuple's default hash
        # would hash every field, including `comp` (a dict, unhashable:
        # confirmed live, `TypeError: unhashable type: 'dict'` straight out
        # of mp_utils's `{i: 1 for i in indices}`). `row_id` alone is
        # already a unique identity for a task within one run_tasks() call
        # (never duplicated within a single round), so hashing/comparing on
        # it alone is correct, not just a workaround.
        return hash(self.row_id)

    def __eq__(self, other):
        return isinstance(other, _RowTask) and self.row_id == other.row_id


_cached_melts = None         # per-worker-process cache of the last SUCCESSFULLY-used
_cached_melts_model = None   # MELTSdynamic engine, see _findLiq_task's docstring.


def _findLiq_task(q, task, *, Model = None,
                   fO2_buffer = None, fO2_offset = None, Affinity = False, quiet_engine = True):
    """
    task_fn for findLiq_multi's mp_utils.run_tasks() call (see that
    function's docstring for the full picture, and `_RowTask`'s docstring
    for why `task` carries only this one row's own data instead of indexing
    into full-batch arrays passed via task_kwargs).

    Deliberately does NOT reuse findLiq()'s body: findLiq() itself, and
    findLiq_MELTS() one level deeper, both swallow a non-converging search
    with a bare `except:` and return an empty `{}` -- completely normal
    behaviour for "the liquidus search didn't converge", not an exceptional
    one, so it does NOT raise on the common failure path (confirmed by
    reading findLiq_MELTS: every non-success branch does `return Results`
    where `Results = {}` from the top, and `Results['T_Liq_C'] = ...` is
    only ever set on the genuine success path). Reporting that as "OK"
    here would defeat mp_utils.run_tasks's CRASHED detection and this
    function's own "retry at a higher starting T" bookkeeping, so this
    wrapper calls findLiq_MELTS directly and explicitly raises when it
    comes back without a liquidus temperature.

    Engine reuse: when this worker runs a SECOND (or later) task --
    relevant only once max_tasks_per_worker > 1 -- this reuses the SAME
    MELTSdynamic engine via `melts.addNodeAfter()`, the same pattern
    multi_path/saturation_pressure already use, instead of (what this used
    to do) constructing a brand-new MELTSdynamic() from scratch every
    single call. The difference matters: a brand-new MELTSdynamic()
    reinitialises the whole C library/engine (loading thermodynamic data
    tables etc.), not just cheap Python-level state, so discarding and
    rebuilding it on every task was paying that cost repeatedly even when
    the whole point of reusing the worker PROCESS was to avoid repeated
    initialisation. Confirmed safe to do here specifically because
    _worker_loop retires (exits) this process immediately on ANY exception
    -- so by the time a second task ever runs in this same process, the
    FIRST one must have succeeded; this cache is only ever read (and only
    ever written) on the success path, so reuse here is never reached
    after an incomplete/failed search, which is the one case confirmed
    UNSAFE earlier this session (a real crash: `TypeError: '>' not
    supported between instances of 'MELTSstatus' and 'int'`, from reusing
    a fresh MELTSdynamic right after a give-up). Reuse after success is
    the one case separately confirmed SAFE (three successive successful
    calculations in one process, all converging correctly).
    """
    global _cached_melts, _cached_melts_model

    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's docstring for why
        # the old `with _quiet_engine_output():` here is intentionally gone, not just unused.

    melts = None
    if _cached_melts is not None and _cached_melts_model == Model:
        melts = _cached_melts.addNodeAfter()

    if Affinity is True:
        Results, Affin, melts = findLiq_MELTS(P_bar = task.P_bar, Model = Model, T_C_init = task.T_initial_C,
                                               comp = task.comp, melts = melts, fO2_buffer = fO2_buffer,
                                               fO2_offset = fO2_offset, Affinity = True, bulk_return = True)
    else:
        Results, melts = findLiq_MELTS(P_bar = task.P_bar, Model = Model, T_C_init = task.T_initial_C,
                                        comp = task.comp, melts = melts, fO2_buffer = fO2_buffer,
                                        fO2_offset = fO2_offset, Affinity = False, bulk_return = True)
        Affin = None

    if not Results or 'T_Liq_C' not in Results:
        raise RuntimeError(f"findLiq_MELTS did not converge at T_init={task.T_initial_C:.1f} C (empty result)")

    # Only cache a melts object that just produced a genuine success -- see
    # this function's own docstring for why that is also the only time
    # this cache is ever read.
    _cached_melts = melts
    _cached_melts_model = Model

    if Affinity is True:
        q.put([Results, Affin, task.row_id])
    else:
        q.put([Results, task.row_id])


def findLiq_multi(cores = None, Model = None, bulk = None, T_initial_C = None, P_bar = None,
                     Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None,
                     Affinity = False, quiet_engine = True,
                     timeout = None, retry_T_offset_C = 25.0, max_domain_retries = 1,
                     max_tasks_per_worker = 1, retry_timeout_multiplier = 1.5,
                     retry_max_tasks_per_worker = 1, Suppress = None):
    '''
    Performs multiple liquidus calculations for a batch of compositions
    and/or P-fO2-H2O conditions, using either MELTS (findLiq_MELTS, via the
    persistent-worker supervisor in mp_utils.py) or MAGEMin
    (MAGEMinCalc.findLiq_julia, via the main embedded Julia session) --
    chosen from `Model` exactly as the rest of the package does ("MELTS" in
    Model -> MELTS; otherwise MAGEMin). The two engines share this
    function's input normalisation and 17-column output shape, but their
    search algorithms, parallelism, and several kwargs below are specific
    to one engine -- see the MAGEMin-specific caveats under `cores`,
    `timeout`, the `retry_*`/`max_tasks_per_worker` family, `quiet_engine`
    and `Affinity`.

    This MELTS implementation replaced an earlier one (one fresh process per group per
    round, round-barrier scheduling, no per-task timeout recovery) after a
    full correctness/performance comparison on the real LiquidusTest_all
    dataset (2048 rows, MAINTENANCE_NOTES.md 2026-10-10): this version is
    2.6x faster (635.6s vs 1678.8s), converges on strictly more rows (2044
    vs 2032 -- zero regressions), and was independently verified more
    accurate on every row where the two disagreed (checked by equilibrating
    at each reported T_Liq_C and confirming the liquid is actually near-
    liquidus -- near-zero solid fraction -- whereas the old implementation's
    disagreeing answers showed exactly zero crystallisation even several
    degrees below their own reported "liquidus", i.e. nowhere near a real
    phase boundary). Behaviourally:

    1. A hung calculation (the C engine stuck inside calcEquilibriumState)
       is killed after `timeout` seconds, but a worker that finishes early
       immediately picks up the next not-yet-attempted composition instead
       of idling until every worker in a fixed-size round has finished or
       been killed -- a persistent pool of workers pulls from a shared
       queue, rather than "spawn `cores` fresh processes, wait for the
       slowest of them, repeat".
    2. A calculation that fails -- either killed on timeout, or
       findLiq_MELTS simply not converging, which (see _findLiq_task's
       docstring) is the common, "successful-looking" failure mode, not an
       exception -- is NOT retried immediately/individually. Every
       still-pending index runs as one batch per round (one `run_tasks`
       call per round, which blocks until that whole batch is done); only
       once a round finishes are its failures collected and handed to the
       NEXT round together. Each domain-retry round changes THREE things
       at once, not just the starting temperature: `T_initial_C +
       retry_T_offset_C` (compounding each round), the per-task `timeout`
       multiplied by `retry_timeout_multiplier` (compounding each round --
       a calculation that timed out may simply have needed more time, not
       a different guess), and `max_tasks_per_worker` forced down to
       `retry_max_tasks_per_worker` (full process isolation -- a
       calculation that already failed once is more likely to be
       genuinely slow/pathological, not just squeezed by contention, so
       it gets a process with no neighbours competing for the same cores).
       Up to `max_domain_retries` rounds. This is a liquidus-search-
       specific heuristic, unrelated to mp_utils.run_tasks's own (disabled
       here, max_retries=0) same-conditions retry.
    3. A worker CAN run up to `max_tasks_per_worker` SUCCESSFUL calculations
       before it is retired and replaced, reusing the same MELTS engine via
       `melts.addNodeAfter()` rather than constructing a fresh one -- in
       principle avoiding both the process-spawn cost AND the engine-init
       cost of starting a new process per row. Both are confirmed safe for
       successive successful calculations (bit-identical results to a
       fresh process/engine) and a calculation that fails always retires
       its worker immediately regardless of this setting, so a bad engine
       state can never carry into the next calculation. HOWEVER: the
       default is 1 (no reuse), because at real batch scale reuse measured
       SLOWER and LESS reliable, not faster -- see `max_tasks_per_worker`'s
       own entry below for the measurements, and MAINTENANCE_NOTES.md: raw
       engine construction itself is cheap (~1.3ms), so the real driver is
       CPU contention under sustained multi-core saturation, which reuse
       makes worse, not process/engine-init overhead, which reuse avoids.
       The knob is kept (not removed) for anyone on lightly-loaded/high-
       core-count hardware where the isolated-test benefit might actually
       materialise, but it is not the default trade-off for this workload.
    4. With quiet_engine=True (the default), a worker's stdout/stderr are
       redirected to the null device ONCE and never restored, rather than
       per-call -- confirmed necessary: the MELTS engine buffers its own
       stderr output internally and only flushes it at process exit (or
       when a SIGTERM-triggered exit runs), which is AFTER a narrow
       per-call redirect would already have put the real fd 2 back. That
       delayed flush is also what let "Caught SIGTERM: a termination
       request was sent to the program" through on every timeout-kill --
       both are now fully suppressed (confirmed against 5/5 real SIGTERM
       kills mid-calculation, zero leaked bytes). A tqdm progress bar
       tracks overall completion (one tick per row that converges, or that
       permanently fails after exhausting its domain retries) independent
       of this -- the bar uses the parent process's own stdout, which is
       never touched by the above.
    5. Each task carries only its own row's data (a `_RowTask`: row id,
       P_bar, T_initial_C, comp), not the whole batch, passed through
       mp_utils.run_tasks's `indices` (pickled once per task) rather than
       its shared `task_kwargs` (which would otherwise be re-pickled and
       re-sent to every single worker spawn). Confirmed necessary: with the
       whole batch in task_kwargs, measured spawn cost was ~0.01s at N=40
       but ~1.1-1.5s at N=1000 -- a batch that should take minutes was
       projecting multiple hours, purely from this. See `_RowTask`'s own
       docstring for the full measurement and the mp_utils.run_tasks
       drain-loop comment for a related race this also made much cheaper
       (a worker that already succeeded could have its exit misread as a
       crash before its success message was drained; harmless to the
       returned data, since that goes through on_progress separately, but
       it wasted a respawn every time -- cheap at small N, expensive at
       large N before this fix).

    Parameters
    ----------
    cores : int, optional
        MELTS only. Number of worker processes to run concurrently. If not
        provided, uses all available CPU cores. Meaningless for MAGEMin
        (findLiq_julia is one Julia call, internally thread-parallel across
        every row from a single shared context -- no per-row processes).
    Model : str, optional
        The thermodynamic model to use: a MELTS variant (e.g.,
        "MELTSv1.0.2", "pMELTS") or a MAGEMin variant ("Green2025" or
        "Weller2024"), matching the rest of the package's "MELTS" in Model
        convention. Default is "MELTSv1.0.2".
    bulk : pd.DataFrame or dict
        A DataFrame (preferred for multiple runs) or dictionary containing
        oxide compositions (wt%).
    T_initial_C : float or np.ndarray, optional
        Initial guess temperature in Celsius for the liquidus search.
        Estimated from the composition (Putirka 2008 eq22, +60 C) if not
        given.
    P_bar : float or np.ndarray
        Pressure(s) in bars for the calculations. If an array, the length
        must match the number of compositions if `bulk` is a DataFrame.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or np.ndarray, optional
        Overrides for initial liquid redox state and volatile contents (wt%).
    fO2_buffer : str, optional
        Oxygen fugacity buffer ("FMQ", "NNO" or "IW").
    fO2_offset : float or np.ndarray, optional
        Offset from the specified fO2 buffer (in log units).
    Suppress : list of str, optional
        MAGEMin only (findLiq_MELTS has no suppress support, unlike
        equilibrate_MELTS -- raises ValueError if set for a MELTS Model
        rather than being silently ignored). List of phases (e.g.,
        ["ol"]) to exclude from the liquidus search, wired through
        MAGEMin_C's `remove_phases`.
    Affinity : bool, optional, default False
        MELTS only. If True, returns an additional DataFrame containing
        chemical affinity (driving force for crystallisation) data for all
        phases. Raises ValueError if True for a MAGEMin Model -- there is
        no MAGEMin equivalent, so this is never silently ignored.
    quiet_engine : bool, default True
        MELTS only. If True, the MELTS engine's own progress text is
        permanently redirected to the null device for each worker process
        (see behavioural note 4 above). Set to False to see its raw
        per-step output for debugging. Meaningless for MAGEMin.
    timeout : float or None, default None
        MELTS only (bounds the mp_utils worker supervisor; a MAGEMin call
        has no comparable open-ended failure mode to bound -- it's one
        deterministic, bounded bracket+bisection search, see
        findLiq_julia's own docstring in Crystallisation.jl).
        Per-calculation budget in seconds for the FIRST attempt, before a
        stuck worker is killed. None (the default) resolves from the
        platform: 20s on macOS/Linux, 45s on Windows (Windows's process
        spawn and MELTS-engine start-up are both measurably slower -- see
        the Round 5-7 Windows benchmarking in MAINTENANCE_NOTES.md; Linux
        is defaulted the same as macOS since both are fork/spawn-cheap
        POSIX platforms and no Linux measurement was taken specifically --
        flag if that's wrong for your hardware). Each domain-retry round's
        timeout is the previous round's multiplied by
        `retry_timeout_multiplier`.
    retry_T_offset_C : float, default 25.0
        MELTS only. Added to a failed calculation's starting temperature
        for each subsequent domain-retry round (compounding: round 2 is
        +25, round 3 is +50, ...).
    max_domain_retries : int, default 1
        MELTS only. How many extra rounds (T_initial_C + retry_T_offset_C, timeout
        * retry_timeout_multiplier, max_tasks_per_worker forced to
        retry_max_tasks_per_worker -- all three change together, not just
        temperature) to attempt for calculations still failing. Every
        round is a single batch: one `run_tasks` call over everything
        still pending, which blocks until that whole round is done; only
        the round's failures (collected together, not retried one at a
        time as they occur) are carried into the next round.
    max_tasks_per_worker : int, default 1
        MELTS only. How many successful calculations one worker process runs in the
        FIRST round before it is retired and replaced with a fresh one --
        in principle, the knob that trades process-spawn/engine-init
        overhead against timeout risk under CPU contention. Raising it
        reduces how often a new process has to be started, and in total
        isolation (cores=1, zero contention, 10 identical fast
        compositions) reuse measured a real 7.2x speedup (max=1 11.4s vs
        max=10 1.6s). That benefit does NOT survive real multi-core
        batches, though -- measured directly on a real 150-row batch
        (retries disabled, so this is round-0-only, with no retry-round
        cost mixed in):

          cores=8:  max=1 -> 38.4s, 150/150 converged
                    max=2 -> 56.3s, 140/150 converged
                    max=3 -> 56.1s, 136/150 converged
          cores=4:  max=1 -> 61.0s, 150/150 converged
                    max=2 -> 66.7s, 146/150 converged
                    max=3 -> 59.1s, 142/150 converged
                    max=5 -> 68.7s, 139/150 converged
                    max=10-> 70.9s, 138/150 converged

        At both core counts max=1 is both the fastest AND the most
        reliable; every level of reuse above 1 is slower AND loses
        conversions (to contention-driven timeouts -- confirmed spurious,
        not a correctness issue: most re-converge fine on their own, and
        WHICH rows fail changes when timeout is raised, the signature of
        contention rather than a genuinely bad composition). Reducing
        cores does not create a sweet spot for reuse either: it just makes
        every setting slower across the board. So although the parameter
        is kept (for anyone on hardware with enough headroom that
        contention stops dominating), 1 -- i.e. no reuse -- is the measured
        best default for real batches. Only applies to the first round --
        every domain-retry round uses `retry_max_tasks_per_worker` instead,
        regardless of this value.
    retry_timeout_multiplier : float, default 1.5
        MELTS only. Each domain-retry round's `timeout` is the previous round's
        multiplied by this (compounding -- round 2 is 1.5x the first
        round's timeout, round 3 is 2.25x, ...).
    retry_max_tasks_per_worker : int, default 1
        MELTS only. `max_tasks_per_worker` forced for every domain-retry round,
        regardless of the first round's setting. A calculation that's
        already failed once is more likely to be genuinely slow or
        pathological rather than just squeezed by contention, so its
        retry attempt(s) get a process with no neighbours competing for
        the same cores -- the safety margin that matters most for exactly
        the compositions most likely to need it.

    Returns
    -------
    Res : pd.DataFrame
        DataFrame containing the calculated liquidus temperature
        (T_Liq_C), the liquidus phase, fluid saturation status, and the
        normalised melt composition at the liquidus. Rows that never
        converged even after all retries are NaN.
    Af_Combined : pd.DataFrame, optional
        Returned only if `Affinity=True`. Contains chemical affinity data
        for the liquidus phase.
    '''
    if timeout is None:
        timeout = 45.0 if platform.system() == "Windows" else 20.0

    comp = bulk.copy()

    if fO2_buffer is not None:
        if fO2_buffer not in ("NNO", "FMQ", "IW"):
            raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if Model is None:
        Model = "MELTSv1.0.2"

    if Affinity is True and "MELTS" not in Model:
        raise ValueError("Affinity=True is MELTS-specific (chemical affinity has no MAGEMin "
                          "equivalent here); it is meaningless for a MAGEMin Model and is not "
                          "silently ignored.")

    if Suppress is not None and "MELTS" in Model:
        raise ValueError("Suppress is MAGEMin-only here -- findLiq_MELTS has no suppress "
                          "support (unlike equilibrate_MELTS), so this is not silently "
                          "ignored for a MELTS Model.")

    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq)

    # --- normalise T_initial_C / P_bar / n for either a single dict
    # composition or a DataFrame batch ---
    if type(comp) != dict or type(P_bar) == np.ndarray:
        if type(comp) != dict:
            n = len(comp['SiO2_Liq'])
        else:
            n = len(P_bar)
    else:
        n = 1

    if T_initial_C is None:
        if type(comp) != dict or type(P_bar) == np.ndarray:
            T_initial_C = estimate_t(comp = comp, P_bar = P_bar).values + 60.0
        else:
            T_initial_C = np.array([estimate_t(comp = comp, P_bar = P_bar) + 60.0])
    elif type(T_initial_C) == int or type(T_initial_C) == float:
        T_initial_C = np.zeros(n) + T_initial_C
    T_initial_C = np.asarray(T_initial_C, dtype = float)

    if type(comp) == dict:
        comps = [comp] * n
    else:
        comps = [comp.loc[i].to_dict() for i in range(n)]

    if type(P_bar) == np.ndarray:
        if len(P_bar) != n:
            raise Warning("The length of your composition and pressure variables are different. Please correct this before running the code.")
        P_list = list(P_bar)
    else:
        P_list = [P_bar] * n

    if "MELTS" in Model:
        cores = resolve_cores(cores)

        pending_tasks = [_RowTask(i, P_list[i], T_initial_C[i], comps[i]) for i in range(n)]
        collected: dict = {}   # index -> (Results, Affin_or_None)
        retry_log: list = []   # (index, round, status, error, T_used, timeout_used) -- for diagnostics/benchmarking

        # One tick per row that is DONE, whether that's a converged result (ticked
        # live, as soon as its _PROGRESS message arrives -- see mp_utils.run_tasks's
        # on_progress) or a permanent failure (ticked only once max_domain_retries
        # is exhausted for it -- a row still queued for another retry round isn't
        # "done" yet, so it must not tick early). `ticked` guards against a rare
        # but real race at the mp_utils level: a task that converges right at the
        # timeout boundary can have its result message sent just before the
        # supervisor's own bookkeeping kills it as a timeout anyway, so the
        # domain-retry logic re-queues and re-runs an index that already
        # succeeded -- harmless to the returned results (collected[idx] just gets
        # overwritten with an equivalent value) but would otherwise double-tick
        # the bar for that row.
        pbar = tqdm(total = n, desc = "findLiq_multi")
        ticked: set = set()

        def _capture(item):
            if Affinity is True:
                Results, Affin, idx = item
                collected[idx] = (Results, Affin)
            else:
                Results, idx = item
                collected[idx] = (Results, None)
            if idx not in ticked:
                ticked.add(idx)
                pbar.update(1)

        current_timeout = timeout
        try:
            for round_num in range(max_domain_retries + 1):
                if not pending_tasks:
                    break

                # First round: the caller's own max_tasks_per_worker, at `timeout`.
                # Every domain-retry round: full process isolation
                # (retry_max_tasks_per_worker, default 1) at a timeout that's
                # already been multiplied up for this round -- not just a
                # different starting temperature. All of pending_tasks runs as
                # ONE batch (one run_tasks call, blocking until the whole round is
                # done); failures are collected and only handed to the NEXT round
                # together, never retried individually/immediately.
                round_max_tasks = max_tasks_per_worker if round_num == 0 else retry_max_tasks_per_worker

                out = run_tasks(
                    _findLiq_task,
                    pending_tasks,  # list of _RowTask -- each carries only its own row's
                                    # data, NOT the full batch (see _RowTask's docstring)
                    task_kwargs = dict(Model = Model, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset,
                                        Affinity = Affinity, quiet_engine = quiet_engine),
                    cores = cores,
                    timeout = current_timeout,
                    max_retries = 0,           # domain-level retry (below) replaces this
                    max_tasks_per_worker = round_max_tasks,  # see this function's own docstring
                                                # for the throughput-vs-contention trade-off and the
                                                # measurements behind the defaults.
                    on_progress = _capture,
                )

                next_round = []
                for task, res in out.items():
                    if res.status is TaskStatus.OK:
                        continue
                    retry_log.append((task.row_id, round_num, res.status.value, res.error,
                                       task.T_initial_C, current_timeout))
                    if round_num < max_domain_retries:
                        next_round.append(_RowTask(task.row_id, task.P_bar,
                                                    task.T_initial_C + retry_T_offset_C, task.comp))
                    elif task.row_id not in ticked:
                        ticked.add(task.row_id)
                        pbar.update(1)   # exhausted its retries -- permanently failed, now "done"

                pending_tasks = next_round
                current_timeout = current_timeout * retry_timeout_multiplier  # for the NEXT round, if any
        finally:
            pbar.close()

        n_ok = len(collected)
        n_failed = n - n_ok
        print(f"findLiq_multi: {n_ok}/{n} converged, {n_failed} did not "
              f"(after {max_domain_retries} retry round(s) at +{retry_T_offset_C} C and "
              f"x{retry_timeout_multiplier} timeout each).")

        Results_df = pd.DataFrame(
            data = np.full((n, 17), np.nan),
            columns = ['T_Liq_C', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3',
                       'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'],
        )
        Results_df[['liquidus_phase', 'fluid_saturated']] = Results_df[['liquidus_phase', 'fluid_saturated']].astype("string")

        Affin = {}
        for idx, (Res, Af) in collected.items():
            Results_df.loc[idx, :] = Res
            if Affinity is True and Af is not None:
                Affin[str(idx)] = Af
    else:
        # --- MAGEMin branch: MAGEMinCalc.findLiq_julia, via the main
        # embedded Julia session (no mp_utils/Distributed -- see that
        # function's own docstring in Crystallisation.jl: multi_point_
        # minimization is internally thread-parallel across all n points
        # from a single Initialize_MAGEMin context, so a batch of
        # thousands is one call, not one process per row). `cores`,
        # `timeout`, `retry_*`, `max_tasks_per_worker`, `quiet_engine` and
        # `Affinity` are MELTS-only knobs for the process-isolation-after-
        # an-opaque-failure problem that doesn't exist here (a single
        # synchronous, deterministic Julia call) -- silently ignored, not
        # remapped onto anything.
        _ensure_julia_main_ready()
        from juliacall import Main as jl

        fO2_buffer_jl = fO2_buffer
        if fO2_buffer == "FMQ":
            fO2_buffer_jl = "qfm"
        elif fO2_buffer == "NNO":
            fO2_buffer_jl = "nno"

        # Same ferric-iron -> "O" placeholder convention as
        # equilibrate_multi's existing MAGEMin branch (Liq.py, ~line 409):
        # MAGEMin has no separate Fe2O3/FeO oxides, just a single "FeO"
        # slot (total Fe as FeO) plus this "O" entry that encodes the
        # ferric fraction.
        fe_c = (159.69 / 2) / 71.844
        if Model == "Weller2024":
            bulk_rows = [[d['SiO2_Liq'], d['Al2O3_Liq'], d['CaO_Liq'], d['MgO_Liq'], d['FeOt_Liq'],
                          d['K2O_Liq'], d['Na2O_Liq'], d['TiO2_Liq'],
                          d['Fe3Fet_Liq'] * (fe_c * d['FeOt_Liq'] - d['FeOt_Liq']), d['Cr2O3_Liq']]
                         for d in comps]
        else:
            bulk_rows = [[d['SiO2_Liq'], d['Al2O3_Liq'], d['CaO_Liq'], d['MgO_Liq'], d['FeOt_Liq'],
                          d['K2O_Liq'], d['Na2O_Liq'], d['TiO2_Liq'],
                          d['Fe3Fet_Liq'] * (fe_c * d['FeOt_Liq'] - d['FeOt_Liq']), d['Cr2O3_Liq'],
                          d['H2O_Liq']]
                         for d in comps]

        bulk_jl = jl.seval("collect")(np.array(bulk_rows, dtype = float))
        P_kbar_jl = jl.seval("collect")(np.asarray(P_list, dtype = float) / 1000.0)
        T_start_jl = jl.seval("collect")(T_initial_C)

        if fO2_offset is None:
            fO2_offset_jl = 0.0
        elif isinstance(fO2_offset, np.ndarray):
            fO2_offset_jl = jl.seval("collect")(fO2_offset.astype(float))
        else:
            fO2_offset_jl = float(fO2_offset)

        # See equilibrate_multi's identical conversion for why this can't
        # just be the plain Python list -- juliacall passes that as a lazy
        # PyList{Any}, which fails findLiq_julia's Vector{String} type.
        suppress_jl = jl.seval("Vector{String}")(Suppress) if Suppress is not None else None

        jl_result = jl.MAGEMinCalc.findLiq_julia(bulk = bulk_jl, P_kbar = P_kbar_jl, T_start_C = T_start_jl,
                                                    Model = Model, fo2_buffer = fO2_buffer_jl, fo2_offset = fO2_offset_jl,
                                                    suppress = suppress_jl)
        jl_result = dict(jl_result)
        T_Liq_C_out = list(jl_result["T_Liq_C"])
        converged_out = list(jl_result["converged"])
        liquidus_phase_out = list(jl_result["liquidus_phase"])
        fluid_saturated_out = list(jl_result["fluid_saturated"])
        liquid_oxides_out = list(jl_result["liquid_oxides"])
        liquid_comp_out = [list(row) for row in jl_result["liquid_comp"]]
        ox_idx = {name: j for j, name in enumerate(liquid_oxides_out)}

        n_ok = sum(1 for c in converged_out if c)
        print(f"findLiq_multi: {n_ok}/{n} converged (MAGEMin, {Model}).")

        Results_df = pd.DataFrame(
            data = np.full((n, 17), np.nan),
            columns = ['T_Liq_C', 'liquidus_phase', 'fluid_saturated', 'SiO2', 'TiO2', 'Al2O3',
                       'Fe2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2'],
        )
        Results_df[['liquidus_phase', 'fluid_saturated']] = Results_df[['liquidus_phase', 'fluid_saturated']].astype("string")

        for i in range(n):
            if not converged_out[i]:
                continue
            row = liquid_comp_out[i]

            def _ox(name):
                return row[ox_idx[name]] * 100.0 if name in ox_idx else np.nan

            # Invert the same ferric-iron <-> "O" formula used to build the
            # input bulk, applied to the LIQUID's own converged composition
            # (not the input bulk -- a buffer can shift the melt's
            # oxidation state away from what was specified), to recover
            # Fe2O3/FeO in the convention comp_fix/the shared tail below
            # already expects (matching what MELTS's own engine returns).
            FeOt_pct = _ox('FeO')
            O_pct = _ox('O')
            if FeOt_pct > 0.0 and not np.isnan(O_pct):
                Fe3Fet_i = min(max(O_pct / (FeOt_pct * (fe_c - 1)), 0.0), 1.0)
            else:
                Fe3Fet_i = 0.0
            Fe2O3_pct = Fe3Fet_i * FeOt_pct / 0.8998
            FeO_pct = FeOt_pct - 0.8998 * Fe2O3_pct

            Results_df.loc[i, 'T_Liq_C'] = T_Liq_C_out[i]
            Results_df.loc[i, 'liquidus_phase'] = liquidus_phase_out[i]
            Results_df.loc[i, 'fluid_saturated'] = "Yes" if fluid_saturated_out[i] else "No"
            Results_df.loc[i, 'SiO2'] = _ox('SiO2')
            Results_df.loc[i, 'TiO2'] = _ox('TiO2')
            Results_df.loc[i, 'Al2O3'] = _ox('Al2O3')
            Results_df.loc[i, 'Fe2O3'] = Fe2O3_pct
            Results_df.loc[i, 'Cr2O3'] = _ox('Cr2O3')
            Results_df.loc[i, 'FeO'] = FeO_pct
            Results_df.loc[i, 'MgO'] = _ox('MgO')
            Results_df.loc[i, 'CaO'] = _ox('CaO')
            Results_df.loc[i, 'Na2O'] = _ox('Na2O')
            Results_df.loc[i, 'K2O'] = _ox('K2O')
            Results_df.loc[i, 'H2O'] = _ox('H2O')
            # MnO, P2O5, CO2 stay NaN -- MAGEMin has no such oxides.

        Affin = {}

    # comp_fix only renames MnO/P2O5/CO2 -> *_Liq for MELTS models (they're
    # not in its "Columns_bad" list for a MAGEMin Model, confirmed by
    # reading GenFuncs.py's comp_fix directly) -- do it ourselves so the
    # column selection below works for both engines without branching.
    # No-op for MELTS (comp_fix's own rename just finds nothing left to
    # rename); for MAGEMin these 3 columns stay NaN either way (MAGEMin has
    # no MnO/P2O5/CO2 oxides), just correctly named now.
    Results_df = Results_df.rename(columns={'MnO': 'MnO_Liq', 'P2O5': 'P2O5_Liq', 'CO2': 'CO2_Liq'})

    Res = comp_fix(Model = Model, comp = Results_df, keep_columns = True)
    if type(Res) == dict:
        Res = pd.DataFrame.from_dict(Res, orient = "index").T

    Res = Res[['T_Liq_C', 'liquidus_phase', 'fluid_saturated',
    'SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq','MnO_Liq', 'MgO_Liq', 'CaO_Liq',
    'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq',
    'Fe3Fet_Liq', 'Fe2O3', 'FeO', ]]

    Res = Res.drop(columns = ['Fe2O3', 'FeO'])

    if Affinity is True:
        int_keys = list(map(int, Affin.keys()))
        if int_keys:
            max_index = max(int_keys)
            all_columns = []
            seen_columns = set()
            for d in Affin:
                for col in Affin[d]:
                    if col not in seen_columns:
                        all_columns.append(col)
                        seen_columns.add(col)
                Affin[d] = pd.Series(Affin[d])

            Af_Combined = pd.DataFrame(np.nan, index = range(max_index + 1), columns = all_columns)
            for d in Affin:
                try:
                    df = Affin[d].to_frame().transpose()
                    Af_Combined.loc[int(d)] = df.iloc[0]
                except:
                    continue
        else:
            Af_Combined = pd.DataFrame()

        return Res, Af_Combined
    else:
        return Res


def equilibrate(q, i,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None, Suppress = None, Suppress_except=None, quiet_engine = True):
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
    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's own docstring
    Res = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, Suppress = Suppress, Suppress_except=Suppress_except)
    q.put([Res, i])
    return