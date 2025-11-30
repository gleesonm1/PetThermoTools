import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.Plotting import *
from petthermotools.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.notebook import tqdm, trange
from pathlib import Path

try:
    shell = get_ipython().__class__.__name__
    # print(shell)
    if shell == 'ZMQInteractiveShell':
        from wurlitzer import sys_pipes as pipes # Jupyter notebook or qtconsole
    elif shell == 'TerminalInteractiveShell':
        pipes = None  # Terminal running IPython
    else:
        pipes = None  # Other type (?)
except:
    pipes = None     # Probably standard Python interpreter (or wurlitzer not available)

def multi_path(cores = None, Model = None, bulk = None, comp = None, Frac_solid = None, Frac_fluid = None, 
               T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
               P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
               Fe3Fet_init = None, Fe3Fet_Liq = None, H2O_init = None, H2O_Liq = None, CO2_init = None, CO2_Liq = None, 
               isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, 
               fO2_buffer = None, fO2_offset = None, 
               Print_suppress = None, fluid_sat = False, Crystallinity_limit = None, Combined = None,
               label = None, timeout = None, print_label = True, Suppress = ['rutile', 'tridymite'], Suppress_except=False,
               multi_processing = True):
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

    Returns
    -------
    Results : dict
        Dictionary with each run's label as key. Values are sub-dictionaries of phase/property DataFrames.
        Includes `Input` key summarizing model configuration per run.
    """
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

    if timeout is None:
        timeout = 180

    timeout_main = timeout

    if Frac_solid is False:
        Frac_solid = None

    if Frac_fluid is False:
        Frac_fluid = None

    if bulk is not None:
        comp = bulk.copy()
        
    # set default values if required
    if Model is None:
        Model == "MELTSv1.0.2"

    if "MELTS" in Model:
        try:
            from meltsdynamic import MELTSdynamic
        except:
            Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

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

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq = CO2_init)
    
    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("Adding small amounts of H$_{2}$O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H$_{2}$O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    if cores is None:
        cores = multiprocessing.cpu_count()

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

    qs = []
    q = Queue()
    
    # perform calculation if only 1 calculation is specified
    if One == 1:
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
                                'Suppress': Suppress, 'Suppress_except': Suppress_except})

            p.start()
            try:
                ret = q.get(timeout = 180)
            except:
                ret = []

            TIMEOUT = 5
            start = time.time()
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

            if Print_suppress is None:
                print(" Complete (time taken = " + str(round(time.time() - s,2)) + " seconds)", end = "", flush = True)

            if len(ret) > 0:
                Results, index = ret
                Results = stich(Results, Model = Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
                return Results
            else:
                Results = {}
                return Results
        else:
            Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                     T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, 
                                     dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, 
                                     P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, 
                                     isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, 
                                     fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, 
                                     Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except)
            
            Results = stich(Results, Model = Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
            return Results

    else: # perform multiple crystallisation calculations
        # first make sure everything is the right length
        L = np.sum(Group)
        if P_bar is None:
            P_bar = [None] * int(L)
        elif type(P_bar) == float or type(P_bar) == int:
            P_bar = np.zeros(int(L)) + P_bar
        if T_C is None:
            T_C = [None] * int(L)
        elif type(T_C) == float or type(T_C) == int:
            T_C = np.zeros(int(L)) + T_C

        if P_start_bar is None:
            P_start_bar = [None] * int(L)
        elif type(P_start_bar) == float or type(P_start_bar) == int:
            P_start_bar = np.zeros(int(L)) + P_start_bar
        if T_start_C is None:
            T_start_C = [None] * int(L)
        elif type(T_start_C) == float or type(T_start_C) == int:
            T_start_C = np.zeros(int(L)) + T_start_C

        if P_end_bar is None:
            P_end_bar = [None] * int(L)
        elif type(P_end_bar) == float or type(P_end_bar) == int:
            P_end_bar = np.zeros(int(L)) + P_end_bar
        if T_end_C is None:
            T_end_C = [None] * int(L)
        elif type(T_end_C) == float or type(T_end_C) == int:
            T_end_C = np.zeros(int(L)) + T_end_C

        if dp_bar is None:
            dp_bar = [None] * int(L)
        elif type(dp_bar) == float or type(dp_bar) == int:
            dp_bar = np.zeros(int(L)) + dp_bar
        if dt_C is None:
            dt_C = [None] * int(L)
        elif type(dt_C) == float or type(dt_C) == int:
            dt_C = np.zeros(int(L)) + dt_C

        if P_path_bar is None:
            P_path_bar = [None] * int(L)
        elif len(np.shape(P_path_bar))  == 1:
            P_path_bar = np.vstack([P_path_bar] * int(L))
        if T_path_C is None:
            T_path_C = [None] * int(L)
        elif len(np.shape(T_path_C))  == 1:
            T_path_C = np.vstack([T_path_C] * int(L))

        if fO2_offset is None:
            fO2_offset = [None] * int(L)
        elif type(fO2_offset) == float or type(fO2_offset) == int:
            fO2_offset = np.zeros(int(L)) + fO2_offset

        index_in = np.arange(int(L))
        combined_results = {}
        index_out = np.array([], dtype=int)

        while len(index_out) < len(index_in):
            index = np.setdiff1d(index_in, index_out)
            groups = np.array_split(index, cores)
            non_empty_groups = [g for g in groups if g.size > 0]
            groups = non_empty_groups

            if len(groups[0])< 15:
                timeout = len(groups[0])*timeout_main
            else:
                timeout = 15*timeout_main
                
            processes = []
            Start = time.time()
            for i in range(len(groups)):
                q = Queue()
                p = Process(target = path_multi, args = (q, groups[i]),
                                kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                        'T_C': T_C, 'T_path_C': T_path_C, 'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                        'P_bar': P_bar, 'P_path_bar': P_path_bar, 'P_start_bar': P_start_bar, 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                                        'isenthalpic': isenthalpic, 'isentropic': isentropic, 'isochoric': isochoric, 'find_liquidus': find_liquidus,
                                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'fluid_sat': fluid_sat, 'Crystallinity_limit': Crystallinity_limit,
                                        'Suppress': Suppress, 'Suppress_except': Suppress_except})
                
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
                    idx_chunks, results = res
                    combined_results.update(results)

                index_out = np.hstack((index_out, idx_chunks))

            print(f"Completed {100*len(index_out)/len(index_in)} %")

        results = combined_results

        Results = stich(Res=results, multi=True, Model=Model, Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
        
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

def path_multi(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None,
            T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None,
            P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None,
            isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None,
            fO2_buffer = None, fO2_offset = None, fluid_sat = False, Crystallinity_limit = None,
            Suppress = None, Suppress_except = None, trail = True):
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
                    if trail is not None:
                        Results, tr = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                                T_initial_C = 1400, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i], 
                                                T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i], 
                                                P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i], 
                                                isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, 
                                                find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                                fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, 
                                                Suppress_except=Suppress_except, trail = trail, melts = melts)
                    else:
                        Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                                T_initial_C = 1400, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i], 
                                                T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i], 
                                                P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i], 
                                                isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, 
                                                find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                                fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, 
                                                Suppress_except=Suppress_except, trail = trail, melts = melts)

                else:
                    if trail is not None:
                        Results, tr = path_MELTS(Model = Model, comp = comp.loc[i].to_dict(), Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                                T_initial_C = 1400, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i], 
                                                T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i], 
                                                P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i], 
                                                isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, 
                                                find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                                fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, 
                                                Suppress_except=Suppress_except, trail = trail, melts = melts)
                    else:
                        Results = path_MELTS(Model = Model, comp = comp.loc[i].to_dict(), Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                                T_initial_C = 1400, T_C = T_C[i], T_path_C = T_path_C[i], T_start_C = T_start_C[i], 
                                                T_end_C = T_end_C[i], dt_C = dt_C[i], P_bar = P_bar[i], P_path_bar = P_path_bar[i], 
                                                P_start_bar = P_start_bar[i], P_end_bar = P_end_bar[i], dp_bar = dp_bar[i], 
                                                isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, 
                                                find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset[i], 
                                                fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, 
                                                Suppress_except=Suppress_except, trail = trail, melts = melts)
            else:
                if fO2_offset[i] is None:
                    fO2_offset[i] = 0.0

                if Frac_solid is None:
                    Frac_solid = False

                if type(comp) == dict:
                    comp_julia = jl.seval("Dict")(comp) 
                else:
                    comp_new = comp.loc[i].to_dict()
                    comp_julia = jl.seval("Dict")(comp_new)

                julia_Suppress = None                
                if Suppress == ['rutile', 'tridymite']:
                    Suppress = None
                    julia_Suppress = None

                if Suppress is not None:
                    jl.Suppress = Suppress
                    jl.seval("Suppress = Vector{String}(Suppress)")
                    julia_Suppress = jl.Suppress


                Results_df = jl.MAGEMinCalc.path(
                                    comp=comp_julia, T_start_C=T_start_C[i], T_end_C=T_end_C[i], dt_C=dt_C[i],
                                    T_C=T_C[i], P_start_bar=P_start_bar[i], P_end_bar=P_end_bar[i], dp_bar=dp_bar[i],
                                    P_bar=P_bar[i], T_path_C=T_path_C[i], P_path_bar=P_path_bar[i], frac_xtal=Frac_solid,
                                    Model=Model, fo2_buffer=fO2_buffer, fo2_offset=fO2_offset[i], find_liquidus=find_liquidus,
                                    suppress = julia_Suppress
                                )
                
                Results = dict(Results_df)
                
                tr = True

            idx.append(i)


            results[f"Run {i}"] = Results

            if tr is False:
                break
        except:
            idx.append(i)
            break

    q.put([idx, results])
    return


def path(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, 
         T_C = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, 
         P_bar = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, 
         isenthalpic = None, isentropic = None, isochoric = None, find_liquidus = None, 
         fO2_buffer = None, fO2_offset = None, fluid_sat = None, Crystallinity_limit = None,
         Suppress = None, Suppress_except = False):
    '''
    Crystallisation calculations to be performed in parallel. Calculations may be either isobaric or isochoric.

    Parameters:
    ----------
    q: Multiprocessing Queue instance
        Queue instance to record the output variables

    index: int
        index of the calculation in the master code (e.g., position within a for loop) to aid indexing results after calculations are complete.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
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

    Returns:
    ----------
    Results: Dict
        Dict containing a series of pandas DataFrames that display the composition and thermodynamic properties of each phase.

    index: int
        index of the calculation

    '''
    Results = {}
    if "MELTS" in Model:
        try:
            # Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, Crystallinity_limit = Crystallinity_limit)
            
            if pipes is not None:
                with pipes():
                    Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                         T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, 
                                         dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, 
                                         P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, 
                                         isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, 
                                         fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, 
                                         Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except)
            else:
                Results = path_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, 
                                     T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, 
                                     dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, 
                                     P_end_bar = P_end_bar, dp_bar = dp_bar, isenthalpic = isenthalpic, 
                                     isentropic = isentropic, isochoric = isochoric, find_liquidus = find_liquidus, 
                                     fO2_buffer = fO2_buffer, fO2_offset = fO2_offset, fluid_sat = fluid_sat, 
                                     Crystallinity_limit = Crystallinity_limit, Suppress = Suppress, Suppress_except = Suppress_except)
            q.put([Results, index])
        except:
            q.put([])

        return

    if "MELTS" not in Model:
        # import julia
        # from julia.api import Julia
        # jl = Julia(compiled_modules=False)
        # from julia import MAGEMinCalc
        # # import pyMAGEMINcalc as MM
        # # try:
        #     # Results = MM.path(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_C = T_C, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_bar = P_bar, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, find_liquidus = find_liquidus, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        # if Frac_solid is None:
        #     Frac_solid = False
            
        # Results = MAGEMinCalc.path(comp = comp, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C,
        #                             T_C = T_C, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar,
        #                             P_bar = P_bar, T_path_C = T_path_C, P_path_bar = P_path_bar, frac_xtal = Frac_solid,
        #                             Model = Model, fo2_buffer = fO2_buffer, fo2_offset = fO2_offset, find_liquidus = find_liquidus)
        # q.put([Results, index])

        # import julia
        from juliacall import Main as jl, convert as jlconvert
        env_dir = Path.home() / ".petthermotools_julia_env"
        jl_env_path = env_dir.as_posix()

        jl.seval(f"""
            import Pkg
            Pkg.activate("{jl_env_path}")
            """)
        # import os, pathlib

        # # 1. Where to install Julia + MAGEMin environment (user home dir, persistent)
        # home = str(pathlib.Path.home())
        # env_path = os.path.join(home, ".MAGEMinEnv")

        # jl.seval(f"""
        # import Pkg
        # Pkg.activate("{env_path}")
        # """)

        jl.seval("using MAGEMinCalc")

        if Frac_solid is None:
            Frac_solid = False

        comp_julia = jl.seval("Dict")(comp)
        # T_path_C = T_path_C.tolist()
        if T_path_C is not None: 
            T_path_C_julia = jlconvert(jl.Vector[jl.Float64], T_path_C)
        else:
            T_path_C_julia = T_path_C

        if P_path_bar is not None:
            P_path_bar_julia = jlconvert(jl.Vector[jl.Float64], P_path_bar)
        else:
            P_path_bar_julia = P_path_bar

        julia_Suppress = None
        if Suppress == ['rutile', 'tridymite']:
            Suppress = None
            julia_Suppress = None

        if Suppress is not None:
            jl.Suppress = Suppress
            jl.seval("Suppress = Vector{String}(Suppress)")
            julia_Suppress = jl.Suppress


        Results = jl.MAGEMinCalc.path(
            comp=comp_julia, T_start_C=T_start_C, T_end_C=T_end_C, dt_C=dt_C,
            T_C=T_C, P_start_bar=P_start_bar, P_end_bar=P_end_bar, dp_bar=dp_bar,
            P_bar=P_bar, T_path_C=T_path_C_julia, P_path_bar=P_path_bar_julia, frac_xtal=Frac_solid,
            Model=Model, fo2_buffer=fO2_buffer, fo2_offset=fO2_offset, find_liquidus=find_liquidus,
            suppress = julia_Suppress
        )
        # Results = jl.pyconvert(dict, Results)
        Results_df = dict(Results)
        q.put([Results_df, index])
        # except:
        #     q.put([])
        # return

