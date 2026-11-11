import contextlib
import warnings
import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.GenFuncs import _ensure_julia_ready
from petthermotools.GenFuncs import _ensure_julia_workers
from petthermotools.GenFuncs import _quiet_engine_output
from petthermotools.GenFuncs import _permanently_quiet_engine_output
from petthermotools.GenFuncs import _reset_comp_index
from petthermotools.Plotting import *
from petthermotools.MELTS import *
from petthermotools.MELTS import _MODEL_CODES
from petthermotools.core_config import MAX_WORKERS, resolve_cores
from petthermotools.mp_utils import run_tasks, TaskStatus as _TaskStatus
# Underscore-prefixed deliberately -- see the identical comment in Liq.py's
# own import of this: TaskStatus must not leak into the wildcard-exported
# public API (enum.Enum's own __call__ signature differs across Python
# versions, which made a bare import here break CI on 3.10/3.11 while
# passing on 3.12/3.13). MAINTENANCE_NOTES.md 2026-10-21.
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.auto import tqdm, trange
from pathlib import Path
from queue import Empty
from typing import NamedTuple
import platform
import os

# nGibbs registry and input builder live in nGibbs_bridge.py (single
# instantiation point, shared with __init__.py without circular imports).
from petthermotools.nGibbs_bridge import _nGibbs_models, nGibbsAPI


try:
    shell = get_ipython().__class__.__name__  # noqa: F821
    if shell == 'ZMQInteractiveShell':
        from wurlitzer import sys_pipes as pipes  # Jupyter notebook or qtconsole
    elif shell == 'TerminalInteractiveShell':
        pipes = None  # Terminal running IPython
    else:
        pipes = None  # Other type (?)
except Exception:
    pipes = None  # Standard Python interpreter (or wurlitzer not available)


def multi_path_v1(cores = None, Model = None, bulk = None, comp = None, Frac_solid = None, Frac_fluid = None,
               T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None,
               P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None,
               Fe3Fet_init = None, Fe3Fet_Liq = None, H2O_init = None, H2O_Liq = None, CO2_init = None, CO2_Liq = None,
               isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = False,
               fO2_buffer = None, fO2_offset = None,
               Print_suppress = None, fluid_sat = False, Crystallinity_limit = None, Combined = None,
               label = None, timeout = None, print_label = True, Suppress = ['rutile', 'tridymite'], Suppress_except=False,
               multi_processing = True, quiet_engine = True):
    """
    Perform single or multiple MELTS/MAGEMin crystallization or degassing calculations in series or parallel.

    Supports isothermal, isenthalpic, isentropic, and isochoric processes over fixed or variable pressure and 
    temperature paths. Works with both MELTS and MAGEMin models. Automatically handles composition expansion when
    using arrays for Fe3+/FeT, H2O, or CO2.

    Parameters
    ----------
    cores : int, optional
        Number of CPU cores to use for multiprocessing. Defaults to total available.
    Model : str
        Thermodynamic model. MELTS variants: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", "pMELTS";
        or MAGEMin: "Green2025", "Weller2024".
    bulk or comp : dict or pd.DataFrame
        Starting composition(s) for the model run(s). If `dict`, single composition for all runs.
    Frac_solid : bool, optional
        If True, remove solids at each step (fractional crystallization).
    Frac_fluid : bool, optional
        If True, remove fluids at each step.
    T_C, P_bar : float or np.ndarray, optional
        Fixed temperature/pressure values for isothermal/isobaric runs.
    T_path_C, P_path_bar : np.ndarray or 2D array, optional
        User-defined temperature/pressure paths for each model run.
    T_start_C, T_end_C, dt_C : float or array, optional
        Starting, ending, and increment temperature values for path calculations.
    P_start_bar, P_end_bar, dp_bar : float or array, optional
        Starting, ending, and increment pressure values for path calculations.
    Fe3Fet_init : float or array, optional
        Initial Fe3+/FeT ratio(s). Required if redox state is not specified in `comp` or `fO2_buffer`.
    H2O_init, CO2_init : float or array, optional
        Initial melt volatile contents. Required if not specified in `comp`.
    isenthalpic, isentropic, isochoric : bool, optional
        If True, applies respective thermodynamic constraint.
    find_liquidus : bool, optional
        If True, searches for the liquidus before starting.
    fO2_buffer : {"FMQ", "NNO"}, optional
        Specifies redox buffer for constraining oxygen fugacity.
    fO2_offset : float or array, optional
        Offset (in log units) from the specified fO2 buffer.
    Print_suppress : bool, optional
        If True, suppresses status messages.
    fluid_sat : bool, optional
        If True, starts the system at fluid saturation.
    Crystallinity_limit : float, optional
        Ends run when crystallinity (excluding fluids) exceeds this threshold.
    Combined : unused
        Placeholder for backwards compatibility.
    label : str or array, optional
        If provided, labels output dictionaries using `label`.
    timeout : float, optional
        Max allowed runtime for each subprocess (in seconds).
    print_label : bool, default=True
        If True, prints label with status updates.
    Suppress : list of str
        List of phases to exclude from results.
    Suppress_except : bool
        If True, `Suppress` acts as a whitelist instead of a blacklist.
    multi_processing : bool
        If False, run sequentially in main process.
    quiet_engine : bool, default=True
        If True, the MELTS engine's own progress text (written to stderr, ~400 KB per run) is
        redirected to the null device for speed -- in a notebook/IDE kernel this alone is worth
        roughly 40% on a batch of runs, since every worker's output otherwise funnels through one
        shared pipe. Set to False to see the engine's raw per-step output for debugging; this does
        not affect PTT's own status messages (see `Print_suppress`) or a worker's own traceback.

    Returns
    -------
    Results : dict
        Dictionary with each run's label as key. Values are sub-dictionaries of phase/property DataFrames.
        Includes `Input` key summarizing model configuration per run.
    """
    bulk = _reset_comp_index(bulk)
    comp = _reset_comp_index(comp)

    ## make sure everything is a float
    T_C        = to_float(T_C)
    T_path_C   = to_float(T_path_C)
    T_start_C  = to_float(T_start_C)
    T_end_C    = to_float(T_end_C)
    dt_C       = to_float(dt_C)

    P_bar      = to_float(P_bar)
    P_path_bar = to_float(P_path_bar)
    P_start_bar= to_float(P_start_bar)
    P_end_bar  = to_float(P_end_bar)
    dp_bar     = to_float(dp_bar)

    Fe3Fet_init= to_float(Fe3Fet_init)
    Fe3Fet_Liq = to_float(Fe3Fet_Liq)
    H2O_init   = to_float(H2O_init)
    H2O_Liq    = to_float(H2O_Liq)
    CO2_init   = to_float(CO2_init)
    CO2_Liq    = to_float(CO2_Liq)
    fO2_offset = to_float(fO2_offset)

    T_C        = check_array(T_C)
    T_start_C  = check_array(T_start_C)
    T_end_C    = check_array(T_end_C)
    dt_C       = check_array(dt_C)

    P_bar      = check_array(P_bar)
    P_start_bar= check_array(P_start_bar)
    P_end_bar  = check_array(P_end_bar)
    dp_bar     = check_array(dp_bar)

    Fe3Fet_init= check_array(Fe3Fet_init)
    Fe3Fet_Liq = check_array(Fe3Fet_Liq)
    H2O_init   = check_array(H2O_init)
    H2O_Liq    = check_array(H2O_Liq)
    CO2_init   = check_array(CO2_init)
    CO2_Liq    = check_array(CO2_Liq)
    fO2_offset = check_array(fO2_offset)

    # Resolve bulk→comp alias and set default Model before any routing.
    if bulk is not None:
        comp = bulk.copy()

    if Model is None:
        Model = "MELTSv1.0.2"

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

    if "MELTS" in Model and multi_processing:
        if platform.system() == "Windows" and "juliacall" in sys.modules:
            raise Exception("On Windows MELTS cannot be run after MAGEMin \nhas been initiated. Please reformat your \nnotebook so that all MELTS calculations take place \nbefore running ptt.activate_petthermotools_env()")

    # ── nGibbs neural-network emulator branch ───────────────────────────────
    # Routed BEFORE the PTT length-consistency check so that nGibbs grid mode
    # (len(P_bar) != len(T_C)) is not incorrectly rejected.
    # Model strings start with 'n' (e.g. "nMELTSv1.0.2", "nMELTSv1.0.2NoProp").
    _ng_base = Model[:-6] if Model.endswith('NoProp') else Model
    if Model.startswith('n') and _ng_base in _nGibbs_models:
        return nGibbsAPI(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
               T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C,
               P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar,
               Fe3Fet_init = Fe3Fet_init, Fe3Fet_Liq = Fe3Fet_Liq, H2O_init = H2O_init, H2O_Liq = H2O_Liq, CO2_init = CO2_init, CO2_Liq = CO2_Liq,
               isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus,
               fO2_buffer = fO2_buffer, fO2_offset = fO2_offset,
               Print_suppress = Print_suppress, fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Combined = Combined,
               label = label, Suppress = Suppress, Suppress_except=Suppress_except)

    ## add check to ensure that there arn't multiple lists provided with different lengths
    l = []
    for item in [Fe3Fet_init,Fe3Fet_Liq,H2O_init,H2O_Liq,CO2_init,CO2_Liq,fO2_offset,P_bar,P_start_bar, P_end_bar, dp_bar, T_C, T_start_C, T_end_C, dt_C]:
        if type(item) == np.ndarray:
            l.append(len(item))

    if type(comp) == pd.core.frame.DataFrame:
        l.append(len(comp))

    if len(l) > 1:
        if len(set(l))>1:
            raise ValueError("Error: non-identical length in provided conditions. \nPlease check lists or arrays are the same length and/or specify single values.")

    if len(l) == 0:
        ln = 1
    else:
        ln = l[0]

    # check label is allowed
    if label is not None:
        if label not in ["CO2", "CO2_init", "pressure", "P", "P_bar", "P_start_bar", "fO2", "fO2_offset", "H2O", "H2O_init", "Fe3Fet", "Fe3Fet_init"]:
            raise ValueError(f"Do not recognise the label chosen. Availabily labels are {['CO2', 'CO2_init', 'pressure', 'P', 'P_bar', 'P_start_bar', 'fO2', 'fO2_offset', 'H2O', 'H2O_init', 'Fe3Fet', 'Fe3Fet_init']}")

    # make sure a timeout is specified
    if timeout is None:
        timeout = 300
        print('A default timeout of 5 minutes has been specified. If you are not getting any results try increasing this using the timeout kwarg.')

    timeout_main = timeout

    if Frac_solid is False:
        Frac_solid = None

    if Frac_fluid is False:
        Frac_fluid = None

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                if fO2_buffer != "IW":
                    raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ', 'IW', or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"
        _ensure_julia_ready()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq = CO2_init)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            warnings.warn("Adding small amounts of H2O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H2O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.", stacklevel=2)

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    cores = resolve_cores(cores)

    if "MELTS" not in Model:
        cores = memory_limit(cores = cores)

    if "MELTS" not in Model:
        if ln > 1:
            P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset = control_path_parameters(ln, P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset)

        _ensure_julia_workers()
        from juliacall import Main as jl 
        julia_Suppress = None                
        if Suppress == ['rutile', 'tridymite']:
            Suppress = None
            julia_Suppress = None

        if Suppress is not None:
            jl.Suppress = Suppress
            jl.seval("Suppress = Vector{String}(Suppress)")
            julia_Suppress = jl.Suppress

        if isinstance(comp, dict):
            # If it's already a dict, convert it to a Julia Dict
            julia_comp = jl.seval("Dict")(comp)
        elif hasattr(comp, 'to_dict'):  # Detects a Pandas DataFrame
            # Convert the entire DataFrame into a list of nested dictionaries (one per row)
            # orient='records' turns it into: [{'SiO2': 50, 'MgO': 10}, {'SiO2': 52, 'MgO': 8}, ...]
            py_dict_list = comp.to_dict(orient='records')
            
            # Convert the Python list of dicts into a Julia Vector of Dicts
            # This makes it incredibly easy for Julia to read row 'i'
            julia_comp = [jl.seval("Dict")(row) for row in py_dict_list]
        else:
            julia_comp = comp

        args = {
            "Model": Model,
            "comp": julia_comp,
            "Frac_solid": Frac_solid,
            "Frac_fluid": Frac_fluid,

            "T_C": T_C,
            "T_path_C": T_path_C,
            "T_start_C": T_start_C,
            "T_end_C": T_end_C,
            "dt_C": dt_C,
            "find_liquidus": find_liquidus,

            "P_bar": P_bar,
            "P_path_bar": P_path_bar,
            "P_start_bar": P_start_bar,
            "P_end_bar": P_end_bar,
            "dp_bar": dp_bar,

            "fO2_buffer": fO2_buffer,
            "fO2_offset": fO2_offset,

            "isenthalpic": isenthalpic,
            "isentropic": isentropic,
            "isochoric": isochoric,

            "fluid_sat": fluid_sat,
            "Crystallinity_limit": Crystallinity_limit,
            "length": ln,

            "timeout": timeout,
            "suppress": julia_Suppress}
        
        args = jl.seval("Dict")(args)
        Results = jl.MAGEMinCalc.multi_path_julia(args)
        if ln == 1:
            Results = stich(Res=Results, Model=Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
            return Results
        else:
            Results = stich(Res=Results, multi=True, Model=Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
            for r in Results:
                i = int(r.split(' ')[1].strip())
                if type(comp) == dict:
                    Results[r]['Input'] = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C[i], 'dt_C': dt_C[i],
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar[i], 'dp_bar': dp_bar[i],
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                        'Suppress': Suppress, 'Suppress_except': Suppress_except}
                else:
                    Results[r]['Input'] = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C[i], 'dt_C': dt_C[i],
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar[i], 'dp_bar': dp_bar[i],
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                        'Suppress': Suppress, 'Suppress_except': Suppress_except}
            
            if label is not None:
                new_out = label_results(Results,label)
                return new_out
            else:
                return Results
            
    else:
        # specify the number of calculations to be performed in each sequence
        One = 0
        if type(comp) == pd.core.frame.DataFrame: # simplest scenario - one calculation per bulk composition imported
            A = len(comp['SiO2_Liq'])//cores
            B = len(comp['SiO2_Liq'])%cores
        else:
            if P_bar is not None and type(P_bar) == np.ndarray: #one calculation per P loaded
                A = len(P_bar)//cores
                B = len(P_bar)%cores
            elif T_C is not None and type(T_C) == np.ndarray: # one calculation per T loaded.
                A = len(T_C)//cores
                B = len(T_C)%cores
            elif P_start_bar is not None and type(P_start_bar) == np.ndarray: # one calculation per starting P
                A = len(P_start_bar)//cores
                B = len(P_start_bar)%cores
            elif T_start_C is not None and type(T_start_C) == np.ndarray: # one calculation per starting T
                A = len(T_start_C)//cores
                B = len(T_start_C)%cores
            elif P_path_bar is not None and len(np.shape(P_path_bar)) > 1: # one calculation per P path
                A = np.shape(P_path_bar)[0]//cores
                B = np.shape(P_path_bar)[0]%cores
            elif T_path_C is not None and len(np.shape(T_path_C)) > 1: # one calculation per T path
                A = np.shape(T_path_C)[0]//cores
                B = np.shape(T_path_C)[0]%cores
            elif fO2_offset is not None and type(fO2_offset) == np.ndarray: # one calculation per offset
                A = len(fO2_offset)//cores
                B = len(fO2_offset)%cores

            else: # just one calculation
                One = 1
                A = 1
                B = 0

        Group = np.zeros(A) + cores
        if B > 0:
            Group = np.append(Group, B)

        q = Queue()
        
        # perform calculation if only 1 calculation is specified
        if One == 1:
            Results = {}
            if Print_suppress is None:
                print("Running " + Model + " calculation...", end = "", flush = True)
                s = time.time()

            if multi_processing:
                p = Process(target = path, args = (q, 1),
                            kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_C': T_C, 'T_path_C': T_path_C, 'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_bar': P_bar, 'P_path_bar': P_path_bar, 'P_start_bar': P_start_bar, 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                                    'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                    'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                    'Suppress': Suppress, 'Suppress_except': Suppress_except, 'quiet_engine': quiet_engine})

                p.start()

                start = time.time()
                if "MELTS" in Model:
                    # drain queue while waiting for calculation to finish
                    while True:
                        try:
                            item = q.get(timeout=0.1)

                            run_index, step_i, data = item
                            Results[step_i] = data

                        except Empty:
                            time.sleep(0.001)
                            pass

                        # Timeout condition
                        if time.time() - start > timeout_main:
                            print("Timeout — terminating process")
                            p.terminate()
                            break

                    # If process finished and queue empty → break
                        if not p.is_alive():
                            try:
                                # Drain remaining items
                                while True:
                                    item = q.get_nowait()
                                    run_index, step_i, data = item
                                    Results[step_i] = data
                            except Empty:
                                time.sleep(0.001)
                                pass
                            break

                    p.join()

                    if Results:
                        new_results = compile_results(Results)

                        # stich and return results
                        Results = stich(new_results, Model = Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
                        if Print_suppress is None:
                            print(" Complete (time taken = " + str(round(time.time() - s,2)) + " seconds)", end = "", flush = True)
                        return Results
                    else:
                        print("No calculations returned")
                        return Results


                if Print_suppress is None:
                    print(" Complete (time taken = " + str(round(time.time() - s,2)) + " seconds)", end = "", flush = True)

                Results = stich(Results, Model = Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)

                return Results
            else:
                with (_quiet_engine_output() if quiet_engine else contextlib.nullcontext()):
                    Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
                                        T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C,
                                        dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar,
                                        P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic,
                                        isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus,
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat,
                                        Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except)
                
                if Results:
                    # identify all unique keys to set up output dataframes
                    unique_keys = set()
                    properties = {}

                    for step_data in Results.values():
                        for key, value in step_data.items():
                            if "_keys" not in key:
                                unique_keys.add(key)
                            else:
                                base_key = key[:-5]
                                if base_key not in properties:
                                    properties[base_key] = value
                                # properties[base_key].update(value)

                    unique_keys = sorted(unique_keys)

                    sorted_steps = sorted(Results.keys())
                    n_steps = len(sorted_steps)

                    index_map = {step: i for i, step in enumerate(sorted_steps)}

                    # create new output to populate with results
                    new_results = {}

                    for key in unique_keys:
                        if key == "Conditions":
                            columns = ['temperature', 'pressure', 'h', 's', 'v', 'mass', 'dvdp', 'logfO2']
                            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

                        elif "_prop" not in key:
                            columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO',
                                    'MnO', 'MgO', 'CaO', 'Na2O', 'K2O',
                                    'P2O5', 'H2O', 'CO2']
                            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

                        elif "_key" not in key:
                            columns = properties[key]
                            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

                    # fill output
                    for step, step_data in Results.items():

                        row_i = index_map[step]

                        for key, value in step_data.items():

                            if key.endswith("_keys"):
                                continue

                            if key not in new_results:
                                continue
                            
                            if len(value) == len(new_results[key].iloc[row_i]):
                                new_results[key].iloc[row_i] = value
                            else:
                                print(f"Size mismatch, potentially errored results on step {row_i}")
                    Results = stich(new_results, Model = Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
                return Results

        else: # perform multiple crystallisation calculations
            # first make sure everything is the right length
            L = np.sum(Group)
            P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset = control_path_parameters(L, P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset)

            index_in = np.arange(int(L))
            Results = {}
            index_out = np.array([], dtype=int)

            while len(index_out) < len(index_in):
                index = np.setdiff1d(index_in, index_out)
                groups = np.array_split(index, cores)
                non_empty_groups = [g for g in groups if g.size > 0]
                groups = non_empty_groups

                if len(groups[-1])< 3:
                    timeout = len(groups[-1])*timeout_main
                else:
                    timeout = 3*timeout_main
                
                processes = []
                start = time.time()
                for i in range(len(groups)):
                    q = Queue()
                    p = Process(target = path_multi, args = (q, groups[i]),
                                    kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                            'T_C': T_C, 'T_path_C': T_path_C, 'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                            'P_bar': P_bar, 'P_path_bar': P_path_bar, 'P_start_bar': P_start_bar, 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                                            'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                            'Suppress': Suppress, 'Suppress_except': Suppress_except, 'quiet_engine': quiet_engine})
                    
                    p.start()
                    time.sleep(0.1)
                    
                    processes.append((p, q, groups[i]))

                    
                # drain queue while waiting for calculation to finish
                if "MELTS" in Model:
                    for p, q, group in processes:
                        while True:
                            try:
                                # 1. Prioritize pulling data. This blocks for 0.1s if empty.
                                item = q.get(timeout=0.1)
                                
                                run_index, step_i, data = item
                                Results.setdefault(run_index, {})[step_i] = data

                            except Empty:
                                # 2. We only check exit conditions when the queue is dry.
                                
                                # Condition A: Process finished naturally
                                if not p.is_alive():
                                    break  
                                
                                # Condition B: Global timeout exceeded
                                if time.time() - start > timeout:
                                    print("Timeout — terminating process")
                                    p.terminate()
                                    p.join()  # Best practice: clean up zombie processes
                                    
                                    # Final sweep: Catch anything that slipped in right at termination
                                    while not q.empty():
                                        try:
                                            item = q.get_nowait()
                                            run_index, step_i, data = item
                                            Results.setdefault(run_index, {})[step_i] = data
                                        except Empty:
                                            break
                                    break
                    # for p, q, group in processes:

                    #     while True:
                    #         try:
                    #             item = q.get(timeout=0.1)

                    #             run_index, step_i, data = item
                    #             Results.setdefault(run_index, {})[step_i] = data

                    #         except Empty:
                    #             time.sleep(0.001)
                    #             pass


                    index_out = np.array(list(Results.keys()))
                

                print(f"Completed {100*len(index_out)/len(index_in)} %")

            if "MELTS" in Model:
                full_results = {}
                for run in Results:
                    if Results[run]:
                        del Results[run][-1]
                        if Results[run]:
                            results = Results[run]
                            new_results = compile_results(results)

                            full_results[f"Run {run}"] = new_results

                Results = stich(Res=full_results, multi=True, Model=Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
            else:
                Results = stich(Res=Results, multi=True, Model=Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
        
            for r in Results:
                i = int(r.split(' ')[1].strip())
                if type(comp) == dict:
                    Results[r]['Input'] = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C[i], 'dt_C': dt_C[i],
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar[i], 'dp_bar': dp_bar[i],
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                        'Suppress': Suppress, 'Suppress_except': Suppress_except}
                else:
                    Results[r]['Input'] = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C[i], 'dt_C': dt_C[i],
                                        'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i], 'P_end_bar': P_end_bar[i], 'dp_bar': dp_bar[i],
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                        'Suppress': Suppress, 'Suppress_except': Suppress_except}
            
            if label is not None:
                new_out = label_results(Results,label)
                return new_out
            else:
                return Results


class _PathRowTask(NamedTuple):
    """
    One row's own data for multi_path_v2, carried through
    mp_utils.run_tasks's `indices` -- same rationale as every other _v2
    row task this session built (findLiq_multi_v2's _RowTask,
    equilibrate_multi_v2's _EquilibrateRowTask, saturation_pressure_v2's
    _SaturationRowTask): a few scalars/arrays and one composition dict per
    row, not the whole batch passed through task_kwargs (which would be
    re-pickled and re-sent to every single worker spawn).
    """
    row_id: int
    comp: dict
    T_C: float
    T_path_C: object
    T_start_C: float
    T_end_C: float
    dt_C: float
    P_bar: float
    P_path_bar: object
    P_start_bar: float
    P_end_bar: float
    dp_bar: float
    fO2_offset: float

    def __hash__(self):
        return hash(self.row_id)

    def __eq__(self, other):
        return isinstance(other, _PathRowTask) and self.row_id == other.row_id


def _path_task_v2(q, task, *, Model=None, fO2_buffer=None, Suppress=None, Suppress_except=None,
                   Frac_solid=None, Frac_fluid=None, isenthalpic=None, isentropic=None,
                   isochoric=None, find_liquidus=None, fluid_sat=None, Crystallinity_limit=None, MgO_limit=None,
                   quiet_engine=True):
    """
    task_fn for multi_path_v2's mp_utils.run_tasks() call. Always
    constructs a brand-new MELTSdynamic() -- deliberately never reads or
    writes a cross-task cache, so this is correct regardless of
    max_tasks_per_worker (if raised above 1 by a caller experimenting with
    it, every task STILL gets a fresh engine here; nothing in this
    function reuses one across rows). That keeps multi_path_v2 safe by
    construction rather than by relying on the default not being
    overridden.
    """
    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; must run before MELTSdynamic()
        # is constructed below, same reasoning as every other task_fn this session.

    from meltsdynamic import MELTSdynamic
    melts = MELTSdynamic(_MODEL_CODES.get(Model, 1))

    stop_info = {'complete': False, 'reason': None}
    Results, tr = path_MELTS(
        Model=Model, comp=task.comp, Frac_solid=Frac_solid, Frac_fluid=Frac_fluid,
        T_initial_C=np.random.random() * 200 + 1200,
        T_C=task.T_C, T_path_C=task.T_path_C, T_start_C=task.T_start_C,
        T_end_C=task.T_end_C, dt_C=task.dt_C,
        P_bar=task.P_bar, P_path_bar=task.P_path_bar, P_start_bar=task.P_start_bar,
        P_end_bar=task.P_end_bar, dp_bar=task.dp_bar,
        isenthalpic=isenthalpic, isentropic=isentropic, isochoric=isochoric,
        find_liquidus=find_liquidus, fO2_buffer=fO2_buffer, fO2_offset=task.fO2_offset,
        fluid_sat=fluid_sat, Crystallinity_limit=Crystallinity_limit,
        Suppress=Suppress, Suppress_except=Suppress_except, trail=True, melts=melts,
        q=q, index=task.row_id, MgO_limit=MgO_limit, stop_info=stop_info
    )

    # `tr` ("trail") is NOT a per-row success/failure flag, even though it
    # looks like one -- confirmed by reading how the ORIGINAL v1 worker
    # (path_multi, this file) treats it: `if tr is False: break` only stops
    # trying FURTHER rows in the same engine-reused batch group, a safety
    # measure against continuing to chain off a possibly-troubled engine.
    # It never discards the row that returned False -- that row's own
    # already-streamed q.put() step data is kept. multi_path_v2 never
    # chains engines across rows in the first place (fresh engine per row,
    # by design), so that concern doesn't apply here at all, and treating
    # `tr is False` as a hard failure was simply wrong: path_MELTS sets it
    # False on transient engine status messages (e.g. "Quadratic
    # iterations") that don't actually invalidate the Results already
    # accumulated, so a fully-populated, legitimate run was being reported
    # as CRASHED (confirmed live: "0/4 succeeded" with all 4 rows actually
    # complete, find_liquidus=True). Only a genuinely empty Results means
    # nothing usable came back.
    if not Results:
        raise RuntimeError(f"path_MELTS reported an empty result for row {task.row_id}")

    # A path only counts as completed if it reached its final step or stopped at a requested
    # limit (Crystallinity_limit, MgO_limit). Otherwise the engine failed part-way: raising here
    # reports the row as CRASHED, and the steps already streamed via q.put are kept, so
    # multi_path counts it as "returned partial results after a crash".
    if not stop_info['complete']:
        last = Results[max(Results)]['Conditions']
        raise RuntimeError(f"path for row {task.row_id} stopped early after {len(Results)} steps "
                           f"(last step: T = {last[0]:.1f} C, P = {last[1]:.1f} bar)")

    # Distinct from path_MELTS's own per-step q.put([index, step_i, data]) calls
    # (step_i is always an int) -- lets on_progress tick the bar once per ROW
    # done, not once per internal step, without guessing from step numbers.
    q.put([task.row_id, "__DONE__", None])


def _add_trace_results(Results, trace_params):
    """
    Runs the trace element engine on finished path results (one run or several) and returns
    (Results, trace_conc, trace_d). Shared by every model branch of multi_path.
    """
    from petthermotools.trace_engine.trace_engine import calc_trace_conc
    trace_conc, trace_d = calc_trace_conc(results=Results,
                                          initial_conc=trace_params['initial_conc'],
                                          partition_models=trace_params.get('partition_models'),
                                          D_coeffs=trace_params.get('D_coeffs', "D_Gibson_2010"))
    return Results, trace_conc, trace_d


def multi_path(cores=None, Model=None, bulk=None, comp=None, Frac_solid=None, Frac_fluid=None,
                   T_C=None, T_path_C=None, T_start_C=None, T_end_C=None, dt_C=None,
                   P_bar=None, P_path_bar=None, P_start_bar=None, P_end_bar=None, dp_bar=None,
                   Fe3Fet_init=None, Fe3Fet_Liq=None, H2O_init=None, H2O_Liq=None,
                   CO2_init=None, CO2_Liq=None,
                   isenthalpic=None, isentropic=None, isochoric=None, find_liquidus=False,
                   fO2_buffer=None, fO2_offset=None, fluid_sat=False, Crystallinity_limit=None,
                   label=None, timeout=None, Suppress=['rutile', 'tridymite'], Suppress_except=False,
                   quiet_engine=True, max_tasks_per_worker=1, max_retries=0, stall_timeout=None,
                   multi_processing=True, Print_suppress=None, MgO_limit=None, trace_params = False):
    """
    Runs one or more crystallisation, decompression or heating paths (isobaric, polybaric,
    isochoric, isentropic, ...). This is the function behind all of the path wrappers
    (isobaric_crystallisation, polybaric_crystallisation_path, isothermal_decompression, etc.).

    A batch of calculations is run when `bulk` is a DataFrame (one run per row) and/or when any
    of the conditions below is given as an array (one run per element); all arrays and the
    DataFrame must then have the same length. MELTS paths run in parallel worker processes, each
    path in its own fresh MELTS engine; MAGEMin and nGibbs models use their own back ends.

    Parameters
    ----------
    cores : int or {"performance", "all"}, optional
        Number of worker processes. Default (None) is the number of performance cores.
    Model : str, optional
        Thermodynamic model. Default "MELTSv1.0.2". Options:
            MELTS: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", "pMELTS"
            MAGEMin: "Green2025", "Weller2024"
            nGibbs: "nMELTSv1.0.2", "nMELTSv1.2.0" (add "NoProp" to skip phase properties)
    bulk, comp : dict, pandas.Series or pandas.DataFrame
        Starting composition(s) in wt% oxides with "_Liq" suffixes (e.g. "SiO2_Liq").
        A dict or Series gives one composition; a DataFrame gives one run per row (its index is
        reset to 0, 1, 2, ... if needed). `bulk` is used if both are given.
    Frac_solid : bool, optional
        If True, solid phases are removed after each step (fractional crystallisation).
        Default False (equilibrium crystallisation).
    Frac_fluid : bool, optional
        If True, the fluid phase is removed after each step. Default False.
    T_path_C : array, optional
        Temperature at every step (°C); overrides T_start_C/T_end_C/dt_C. A 2D array gives one
        path per row.
    T_start_C, T_end_C, dt_C : float or array, optional
        Starting temperature, final temperature and temperature step (°C) of a linear path.
    T_C : float or array, optional
        Fixed temperature (°C), used as T_start_C when T_start_C is not given (e.g. isothermal
        decompression).
    find_liquidus : bool, optional
        If True, the path starts at the liquidus temperature found for the starting
        composition and pressure. Either find_liquidus=True or a starting temperature
        (T_start_C, T_C or T_path_C) must be given. Default False.
    P_path_bar : array, optional
        Pressure at every step (bar); overrides P_start_bar/P_end_bar/dp_bar. A 2D array gives
        one path per row.
    P_start_bar, P_end_bar, dp_bar : float or array, optional
        Starting pressure, final pressure and pressure step (bar) of a linear path.
    P_bar : float or array, optional
        Constant pressure (bar) for an isobaric path, or the starting pressure of an isochoric
        path.
    Fe3Fet_init : float or array, optional
        Initial Fe3+/total Fe ratio of the system. If None, the value in the composition is used.
    H2O_init, CO2_init : float or array, optional
        Initial H2O and CO2 contents (wt%). If None, the values in the composition are used.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or array, optional
        Older names for Fe3Fet_init, H2O_init and CO2_init; still accepted.
    fO2_buffer : {"FMQ", "NNO", "IW"}, optional
        Oxygen fugacity buffer to follow along the path. Without a buffer the system is closed
        to oxygen, so the initial Fe3+/total Fe must not be zero.
    fO2_offset : float or array, optional
        Offset from fO2_buffer in log units.
    isochoric : bool, optional
        If True, the path is at constant volume (P_bar is the starting pressure).
    isentropic : bool, optional
        If True, the path is at constant entropy (e.g. adiabatic decompression).
    isenthalpic : bool, optional
        If True, the path is at constant enthalpy.
    fluid_sat : bool, optional
        If True, the calculation starts at fluid saturation. Default False.
    Crystallinity_limit : float, optional
        Stops a run once the volume fraction of crystals (0-1, relative to crystals + liquid)
        exceeds this value.
    MgO_limit : float, optional
        Stops a run once the MgO content of the liquid (wt%) falls below this value (or no
        liquid is left); the first step below the limit is kept. MELTS models only, ignored with
        a warning for MAGEMin and nGibbs.
    label : str, optional
        For a batch, names the runs after the input that varies instead of "Run 0", "Run 1", ...
        Options: 'P_bar' (or 'P', 'pressure'), 'P_start_bar', 'H2O_init' (or 'H2O'),
        'CO2_init' (or 'CO2'), 'Fe3Fet_init' (or 'Fe3Fet'), 'fO2_offset' (or 'fO2').
    Suppress : str or list of str, optional
        Phases that are not allowed to form. Default ['rutile', 'tridymite']; "All" suppresses
        every solid phase.
    Suppress_except : bool or str or list of str, optional
        If given, every solid phase except these is suppressed. Default False.
    timeout : float, optional
        Time limit per run (seconds) before a stuck MELTS run is stopped. Default 60 s
        (180 s on Windows). For MAGEMin the default is 5 minutes.
    max_retries : int, optional
        Number of times a MELTS run that failed or timed out is retried with the same inputs.
        Default 0.
    stall_timeout : float, optional
        Stops (and retries, if max_retries > 0) a MELTS run that has made no progress for this
        many seconds, even if `timeout` has not been reached. Default None (off).
    max_tasks_per_worker : int, optional
        Number of runs a MELTS worker process handles before it is replaced. Default 1.
    quiet_engine : bool, optional
        If True (default), the MELTS engine's own output is hidden, which is faster. Set to False
        to see it when debugging.
    Print_suppress : optional
        If not None, the timeout and summary messages are not printed.
    multi_processing : bool, optional
        Kept for compatibility; MELTS paths always run in worker processes.
    trace_params : dict, optional
        If given, the trace element engine is run on the results (any model). Keys:
        'initial_conc' (required, e.g. {"Ni": 250.0}), 'partition_models' (optional, e.g.
        {"cpx": "Sun2012", "ol": "Sun2014"}) and 'D_coeffs' (optional, a name such as
        "D_Gibson_2010" or a dict of phase -> {element: D}; default "D_Gibson_2010").

    Returns
    -------
    Results : dict
        For a single run, a dict of DataFrames ('All', 'Conditions', 'mass_g', one per phase,
        ...). For a batch, a dict of those keyed "Run 0", "Run 1", ... (or by `label`), each with
        an 'Input' entry recording its inputs.
    trace_conc, trace_d : dict
        Only when `trace_params` is given, in which case (Results, trace_conc, trace_d) is
        returned: trace element concentrations and partition coefficients per phase, keyed the
        same way as Results.

    Raises
    ------
    ValueError
        If neither a starting temperature (T_start_C, T_C or T_path_C) nor find_liquidus=True
        is given, if fO2_buffer or label is not recognised, or if array inputs differ in length.
    Warning
        For a single MELTS composition with zero Fe3+ and no fO2 buffer (MELTS usually fails).

    Warns
    -----
    UserWarning
        For a single MELTS composition with zero H2O (the calculation still runs).

    Examples
    --------
    >>> Results = multi_path(Model="MELTSv1.0.2", bulk=comp, P_bar=2000, find_liquidus=True,
    ...                      T_end_C=900, dt_C=5, Fe3Fet_init=0.15, H2O_init=1.0)
    """
    if Model is None:
        Model = "MELTSv1.0.2"

    bulk = _reset_comp_index(bulk)
    comp = _reset_comp_index(comp)

    # Checked here, before any model branch: without it a MELTS run starts anyway and each row
    # fails inside its worker with an unhelpful TypeError (path_MELTS's own check only fires for
    # find_liquidus=None, but the default here is False)
    if not find_liquidus and T_start_C is None and T_C is None and T_path_C is None:
        raise ValueError("A starting temperature must be specified (T_start_C, T_C or T_path_C), "
                         "or set find_liquidus=True to start at the liquidus.")

    if MgO_limit is not None and "MELTS" not in Model:
        warnings.warn(f"MgO_limit is only implemented for the MELTS models and will be ignored for Model = '{Model}'.")

    if "MELTS" not in Model:
        Results =  multi_path_v1(cores=cores, Model=Model, bulk=bulk, comp=comp, Frac_solid=Frac_solid,
                           Frac_fluid=Frac_fluid, T_C=T_C, T_path_C=T_path_C, T_start_C=T_start_C,
                           T_end_C=T_end_C, dt_C=dt_C, P_bar=P_bar, P_path_bar=P_path_bar,
                           P_start_bar=P_start_bar, P_end_bar=P_end_bar, dp_bar=dp_bar,
                           Fe3Fet_init=Fe3Fet_init, Fe3Fet_Liq=Fe3Fet_Liq, H2O_init=H2O_init,
                           H2O_Liq=H2O_Liq, CO2_init=CO2_init, CO2_Liq=CO2_Liq,
                           isenthalpic=isenthalpic, isentropic=isentropic, isochoric=isochoric,
                           find_liquidus=find_liquidus, fO2_buffer=fO2_buffer, fO2_offset=fO2_offset,
                           fluid_sat=fluid_sat, Crystallinity_limit=Crystallinity_limit,
                           label=label, timeout=timeout, Suppress=Suppress,
                           Suppress_except=Suppress_except, multi_processing=multi_processing,
                           Print_suppress=Print_suppress, quiet_engine=quiet_engine)

        if trace_params:
            return _add_trace_results(Results, trace_params)
        return Results

    if Model.startswith('n'):
        _ng_base = Model[:-6] if Model.endswith('NoProp') else Model
        if _ng_base in _nGibbs_models:
            Results = multi_path_v1(cores=cores, Model=Model, bulk=bulk, comp=comp, Frac_solid=Frac_solid,
                               Frac_fluid=Frac_fluid, T_C=T_C, T_path_C=T_path_C, T_start_C=T_start_C,
                               T_end_C=T_end_C, dt_C=dt_C, P_bar=P_bar, P_path_bar=P_path_bar,
                               P_start_bar=P_start_bar, P_end_bar=P_end_bar, dp_bar=dp_bar,
                               Fe3Fet_init=Fe3Fet_init, Fe3Fet_Liq=Fe3Fet_Liq, H2O_init=H2O_init,
                               H2O_Liq=H2O_Liq, CO2_init=CO2_init, CO2_Liq=CO2_Liq,
                               isenthalpic=isenthalpic, isentropic=isentropic, isochoric=isochoric,
                               find_liquidus=find_liquidus, fO2_buffer=fO2_buffer, fO2_offset=fO2_offset,
                               fluid_sat=fluid_sat, Crystallinity_limit=Crystallinity_limit,
                               label=label, timeout=timeout, Suppress=Suppress,
                               Suppress_except=Suppress_except, Print_suppress=Print_suppress,
                               quiet_engine=quiet_engine)
            if trace_params:
                return _add_trace_results(Results, trace_params)
            return Results

    ## make sure everything is a float -- identical to multi_path's own preamble
    T_C = check_array(to_float(T_C))
    T_path_C = to_float(T_path_C)
    T_start_C = check_array(to_float(T_start_C))
    T_end_C = check_array(to_float(T_end_C))
    dt_C = check_array(to_float(dt_C))

    P_bar = check_array(to_float(P_bar))
    P_path_bar = to_float(P_path_bar)
    P_start_bar = check_array(to_float(P_start_bar))
    P_end_bar = check_array(to_float(P_end_bar))
    dp_bar = check_array(to_float(dp_bar))

    Fe3Fet_init = check_array(to_float(Fe3Fet_init))
    Fe3Fet_Liq = check_array(to_float(Fe3Fet_Liq))
    H2O_init = check_array(to_float(H2O_init))
    H2O_Liq = check_array(to_float(H2O_Liq))
    CO2_init = check_array(to_float(CO2_init))
    CO2_Liq = check_array(to_float(CO2_Liq))
    fO2_offset = check_array(to_float(fO2_offset))

    if bulk is not None:
        comp = bulk.copy()

    if H2O_Liq is not None and H2O_init is None:
        H2O_init = H2O_Liq
    if CO2_Liq is not None and CO2_init is None:
        CO2_init = CO2_Liq
    if Fe3Fet_Liq is not None and Fe3Fet_init is None:
        Fe3Fet_init = Fe3Fet_Liq

    if Frac_solid is False:
        Frac_solid = None
    if Frac_fluid is False:
        Frac_fluid = None

    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    if fO2_buffer is not None and fO2_buffer not in ("NNO", "FMQ", "IW"):
        raise ValueError("fO2 buffer specified is not an allowed input. This argument can only be "
                          "'FMQ', 'IW', or 'NNO' -- if you want to offset from these buffers use the "
                          "'fO2_offset' argument.")

    if label is not None:
        allowed_labels = ['CO2', 'CO2_init', 'pressure', 'P', 'P_bar', 'P_start_bar', 'fO2',
                           'fO2_offset', 'H2O', 'H2O_init', 'Fe3Fet', 'Fe3Fet_init']
        if label not in allowed_labels:
            raise ValueError(f"Do not recognise the label chosen. Available labels are {allowed_labels}")

    if timeout is None:
        timeout = 180.0 if platform.system() == "Windows" else 60.0
        if Print_suppress is None:
            print(f"multi_path_v2: timeout not specified, using the platform default of {timeout:.0f}s "
                  f"per path (pass timeout=<seconds> to override).")
    else:
        if Print_suppress is None:
            print(f"multi_path_v2: using a {timeout:.0f}s timeout per path.")

    comp = comp_fix(Model=Model, comp=comp, Fe3Fet_Liq=Fe3Fet_init, H2O_Liq=H2O_init, CO2_Liq=CO2_init)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0:
            warnings.warn("Adding small amounts of H2O may improve the ability of MELTS to accurately "
                          "reproduce the saturation of oxide minerals. Additionally, sufficient H2O is "
                          "required in the model for MELTS to predict the crystallisation of apatite, "
                          "rather than whitlockite.", stacklevel=2)
        if comp['Fe3Fet_Liq'] == 0.0 and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in "
                           "the starting composition and an fO2 buffer is not set.")

    cores = resolve_cores(cores)

    ## figure out how many runs are being requested -- identical logic to multi_path's own
    l = []
    for item in [Fe3Fet_init, Fe3Fet_Liq, H2O_init, H2O_Liq, CO2_init, CO2_Liq, fO2_offset,
                 P_bar, P_start_bar, P_end_bar, dp_bar, T_C, T_start_C, T_end_C, dt_C]:
        if type(item) == np.ndarray:
            l.append(len(item))
    if type(comp) == pd.core.frame.DataFrame:
        l.append(len(comp))
    if len(l) > 1 and len(set(l)) > 1:
        raise ValueError("Error: non-identical length in provided conditions. Please check lists or "
                          "arrays are the same length and/or specify single values.")
    n = l[0] if len(l) > 0 else 1

    P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset = \
        control_path_parameters(n, P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C,
                                 dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset)

    if type(comp) == dict:
        comps = [comp] * n
    else:
        comps = [comp.loc[i].to_dict() for i in range(n)]

    pending_tasks = [
        _PathRowTask(i, comps[i], T_C[i], T_path_C[i], T_start_C[i], T_end_C[i], dt_C[i],
                     P_bar[i], P_path_bar[i], P_start_bar[i], P_end_bar[i], dp_bar[i], fO2_offset[i])
        for i in range(n)
    ]

    collected: dict = {}   # row_id -> {step_i: step_data, ...}
    pbar = tqdm(total=n, desc="multi_path_v2")
    ticked: set = set()

    def _tick(row_id):
        if row_id not in ticked:
            ticked.add(row_id)
            pbar.update(1)

    def _capture(item):
        row_id, step_i, data = item
        if step_i == "__DONE__":
            _tick(row_id)
            return
        collected.setdefault(row_id, {})[step_i] = data

    try:
        out = run_tasks(
            _path_task_v2,
            pending_tasks,
            task_kwargs=dict(Model=Model, fO2_buffer=fO2_buffer, Suppress=Suppress,
                              Suppress_except=Suppress_except, Frac_solid=Frac_solid,
                              Frac_fluid=Frac_fluid, isenthalpic=isenthalpic, isentropic=isentropic,
                              isochoric=isochoric, find_liquidus=find_liquidus, fluid_sat=fluid_sat,
                              Crystallinity_limit=Crystallinity_limit, quiet_engine=quiet_engine, MgO_limit=MgO_limit),
            cores=cores,
            timeout=timeout,
            stall_timeout=stall_timeout,
            max_retries=max_retries,
            max_tasks_per_worker=max_tasks_per_worker,
            on_progress=_capture,
        )
    finally:
        for row_id in range(n):
            _tick(row_id)   # any row never reaching __DONE__ (failed/timed out) still ticks once
        pbar.close()

    # Four-way breakdown instead of a plain succeeded/didn't: `out[row_id].status`
    # alone can't tell a row that crashed (or timed out) with real, already-
    # streamed partial step data apart from one that returned nothing at all
    # -- both show up identically as GAVE_UP there (run_tasks's own
    # _finalise_or_retry collapses every non-OK status to GAVE_UP once
    # retries are exhausted, which is immediate here since max_retries=0 by
    # default). `collected[row_id]` -- filled from on_progress as steps
    # arrive, independent of the task's eventual status -- is what actually
    # says whether usable data came back; `.error`'s own text (the one place
    # the original TIMEOUT/STALLED/CRASHED distinction survives) says why a
    # non-OK row stopped. STALLED is grouped with TIMEOUT here (both are an
    # external cutoff imposed by the supervisor, not an engine-side failure).
    n_completed = 0
    n_partial_crash = 0
    n_partial_timeout = 0
    n_no_results = 0

    # `out` is keyed by the _PathRowTask objects passed in as `indices`, not
    # by plain row_id ints -- re-key by .row_id once rather than relying on
    # _PathRowTask.__eq__ matching a bare int (it doesn't: it only compares
    # equal to another _PathRowTask).
    result_by_row = {task.row_id: result for task, result in out.items()}

    failures = []
    for row_id in range(n):
        result = result_by_row.get(row_id)
        has_data = bool(collected.get(row_id))
        if result is not None and result.status is _TaskStatus.OK:
            n_completed += 1
            continue
        error = (result.error or "") if result is not None else ""
        if not has_data:
            n_no_results += 1
            failures.append((row_id, "no results", error))
        elif "timed out" in error or "no activity for over" in error:
            n_partial_timeout += 1
            failures.append((row_id, "partial results", error))
        else:
            n_partial_crash += 1
            failures.append((row_id, "partial results", error))

    if Print_suppress is None:
        print(f"multi_path_v2: {n_completed}/{n} completed, {n_partial_crash}/{n} returned partial "
              f"results after a crash, {n_partial_timeout}/{n} returned partial results after timing "
              f"out, {n_no_results}/{n} returned no results (no retry by default -- see max_retries).")
        max_listed = 10
        for row_id, outcome, error in failures[:max_listed]:
            print(f"  Run {row_id} ({outcome}): {error or 'no reason reported'}")
        if len(failures) > max_listed:
            print(f"  ... and {len(failures) - max_listed} more runs that did not complete.")

    full_results = {}
    for row_id, steps in collected.items():
        if steps:
            full_results[f"Run {row_id}"] = compile_results(steps)

    if n == 1:
        # Match multi_path's own single-run return shape: the one run's own
        # dict of DataFrames, not wrapped in a "Run 0" dict of dicts.
        # No early return here: a single run must still reach the trace_params block below
        if full_results:
            Results = stich(list(full_results.values())[0], Model=Model,
                             Frac_fluid=Frac_fluid, Frac_solid=Frac_solid)
        else:
            print("No calculations returned")
            Results = {}
    else:
        Results = stich(Res=full_results, multi=True, Model=Model, Frac_fluid=Frac_fluid, Frac_solid=Frac_solid)

        for r in Results:
            i = int(r.split(' ')[1].strip())
            row_comp = comp if type(comp) == dict else comp.loc[i].to_dict()
            Results[r]['Input'] = {
                'Model': Model, 'comp': row_comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                'T_C': T_C[i], 'T_path_C': T_path_C[i], 'T_start_C': T_start_C[i], 'T_end_C': T_end_C[i],
                'dt_C': dt_C[i], 'P_bar': P_bar[i], 'P_path_bar': P_path_bar[i], 'P_start_bar': P_start_bar[i],
                'P_end_bar': P_end_bar[i], 'dp_bar': dp_bar[i], 'isenthalpic': isenthalpic,
                'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i], 'fluid_sat': fluid_sat,
                'Crystallinity_limit': Crystallinity_limit, 'MgO_limit': MgO_limit,
                'Suppress': Suppress, 'Suppress_except': Suppress_except,
            }

        # label_results needs the multi-run "Run i" layout with 'Input', so a single run is never labelled
        if label is not None:
            Results = label_results(Results, label)

    if trace_params:
        return _add_trace_results(Results, trace_params)
    return Results


# Plain alias: `multi_path` IS this implementation as of 2026-10-21 (see its
# own docstring), so `multi_path_v2` is now redundant but kept so any
# existing code/notebook/scratch script calling it by that name keeps working.
multi_path_v2 = multi_path


def path_multi(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None,
            T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None,
            P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None,
            isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None,
            fO2_buffer = None, fO2_offset = None, fluid_sat = False, Crystallinity_limit = None,
            Suppress = None, Suppress_except = None, trail = True, shared_dict = None, quiet_engine = True,
            MgO_limit = None):
    """
    Worker function to run a subset of crystallization/decompression models (MELTS or MAGEMin) in parallel.

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
    MgO_limit : float
        Ends run when the MgO content of the liquid (wt%) falls below this value.
    Suppress : list of str
        Phases to exclude from results.
    Suppress_except : bool
        If True, treat `Suppress` as a whitelist.
    trail : bool
        If True, include trailing properties from model output.
    quiet_engine : bool, default=True
        If True, redirects the MELTS engine's own stderr output to the null device for the
        duration of this worker's runs (set to False from `multi_path` for debugging).

    Returns
    -------
    None
        The function returns results using `q.put()`:
        q.put([idx, results])
        where `idx` is the list of completed indices and `results` is a dictionary of output per run.
    """

    idx = []

    if len(index) > 3:
        index = index[:3]
    
    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's own docstring --
        # must run before MELTSdynamic() below, whose own construction is what prints the
        # "alphaMELTS for Python version..." banner (confirmed directly).

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
        try:
            if "MELTS" in Model:
                q.put([i, -1, "Start"])
                if i != 0:
                    melts = melts.addNodeAfter()
                if type(comp) == dict:
                    Results, tr = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
                                            T_initial_C = np.random.random()*200 + 1200, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i],
                                            T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i],
                                            P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i],
                                            isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric,
                                            find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i],
                                            fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress,
                                            Suppress_except=Suppress_except, trail = trail, melts = melts, q = q, index = i,
                                            MgO_limit = MgO_limit)
                else:
                    Results, tr = path_MELTS(Model = Model, comp = comp.loc[i].to_dict(), Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
                                            T_initial_C = np.random.random()*200 + 1200, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i],
                                            T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i],
                                            P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i],
                                            isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric,
                                            find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i],
                                            fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress,
                                            Suppress_except=Suppress_except, trail = trail, melts = melts, q = q, index = i,
                                            MgO_limit = MgO_limit)

            if tr is False:
                break
        except:
            idx.append(i)
            break

    return


def path(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, 
         T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
         P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
         isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, 
         fO2_buffer = None, fO2_offset = None, fluid_sat = None, Crystallinity_limit = None,
         Suppress = None, Suppress_except = False, shared_dict=None, quiet_engine = True, MgO_limit = None):
    '''
    Crystallisation calculations to be performed in parallel. Calculations may be either isobaric or isochoric.

    Parameters:
    ----------
    q: Multiprocessing Queue instance
        Queue instance to record the output variables

    index: int
        index of the calculation in the master code (e.g., position within a for loop) to aid indexing results after calculations are complete.

    Model: string
        "MELTS" or a MAGEMin variant ("Green2025", "Weller2024"). Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict
        Initial compositon for calculations.

    Frac_solid: True/False
        If True, solid phases will be removed from the system at the end of each crystallisation step. Default False.

    Frac_fluid: True/False
        If True, fluid phases will be removed from the system at the end of each crystallisation step. Default False.

    T_C: float
        Calculation temperature - typically used when calculations are performed at a fixed T (e.g.,isothermal degassing).

    T_path_C: np.ndarray
        If a specified temperature path is to be used, T_path_C will be used to determine the T at each step of the model. If 2D, this indicates that multiple calculations with different T_path_C arrays are to be performed.

    T_start_C: float
        Initial temperature used for path calculations.

    T_end_C: float
        Final temperature in crystallisation calculations.

    dt_C: float
        Temperature increment during crystallisation calculations.

    P_bar: float
        Calculation pressure - typically used when calculations are performed at a fixed P (e.g.,isobaric crystallisation).

    P_path_bar: np.ndarray
        If a specified pressure path is to be used, P_path_bar will be used to determine the P at each step of the model. If 2D, this indicates that multiple calculations with different P_path_bar arrays are to be performed.

    P_start_bar: float
        Initial pressure used for path calculations.

    P_end_bar: float
        Final pressure in crystallisation calculations.

    dp_bar: float
        Pressure increment during crystallisation calculations..

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio.

    H2O_Liq: float or np.ndarray
        H2O content of the initial bulk composition.

    isenthalpic: True/False
        If True, calculations will be performed at a constant enthalpy with T treated as a dependent variable.

    isentropic: True/False
        If True, calculations will be performed at a constant entropy with T treated as a dependent variable.

    isochoric: True/False
        If True, the volume of the system will be held constant instead of the pressure. Default is False.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float
        Offset from the buffer spcified in fO2_buffer (log units).

    quiet_engine: bool, default=True
        If True, redirects the MELTS engine's own stderr output to the null device for speed
        (set to False from `multi_path` for debugging; overrides the Jupyter wurlitzer capture).

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    index: int
        index of the calculation

    '''
    Results = {}
    if "MELTS" in Model:
        # try:
            # Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit)
            
        # quiet_engine (default) always wins over the Jupyter wurlitzer pipes() capture below: showing the
        # engine's raw text nicely in a notebook cell and suppressing it for speed are mutually exclusive
        # goals, and the speed one is the default. Set quiet_engine=False to get today's old behaviour back
        # (pipes() in a Jupyter/qtconsole kernel, plain/unsuppressed everywhere else) for debugging.
        #
        # The quiet_engine=True case uses the one-way permanent redirect, not a `with`-scoped one
        # (see _permanently_quiet_engine_output's own docstring): the alphaMELTS engine buffers its
        # own stdout/stderr internally, so a scoped redirect lets banner/status text leak through
        # AFTER the block restores the real fd 1/2. The other two branches (pipes()/nullcontext)
        # exist specifically to SHOW output, not hide it, so they keep their own, separate,
        # properly-scoped `with` -- only this function's own worker process exits right after
        # either branch, so only the quiet_engine=True branch is "one-way forever" by design.
        if quiet_engine:
            _permanently_quiet_engine_output()
            Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
                                    T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C,
                                    dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar,
                                    P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic,
                                    isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus,
                                    fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat,
                                    Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except, q = q, index = 0,
                                    MgO_limit = MgO_limit)
        else:
            output_ctx = pipes() if pipes is not None else contextlib.nullcontext()
            with output_ctx:
                Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid,
                                        T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C,
                                        dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar,
                                        P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic,
                                        isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus,
                                        fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat,
                                        Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except, q = q, index = 0,
                                        MgO_limit = MgO_limit)

        return


