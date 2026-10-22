import numpy as np
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.GenFuncs import _quiet_engine_output
from petthermotools.GenFuncs import _permanently_quiet_engine_output
import contextlib
from petthermotools.Plotting import *
from petthermotools.MELTS import *
from petthermotools.MELTS import _MODEL_CODES
from petthermotools.core_config import MAX_WORKERS, resolve_cores
from petthermotools.mp_utils import run_tasks, TaskStatus
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys, os
import platform
import warnings
from typing import NamedTuple
from tqdm.auto import tqdm, trange
import tempfile # use to create unique temp directories for each specific process

def findSatPressure_multi(cores = MAX_WORKERS, Model = "MELTSv1.2.0", bulk = None, T_fixed_C = None,
                          P_bar_init = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None,
                          fO2_buffer = None, fO2_offset = None, copy_columns = None, quiet_engine = True):
    '''
    Performs multiple volatile saturation pressure ($P_{\text{sat}}$) calculations in parallel 
    for a batch of melt compositions at a fixed temperature.

    The saturation pressure is the minimum pressure at which the melt can hold all its 
    volatile components ($\text{H}_2\text{O}$ and $\text{CO}_2$) without exsolving a fluid phase.

    Parameters:
    ----------
    cores : int, optional
        Number of processes to run in parallel. Defaults to all available CPU cores.
    Model : str, optional, default "MELTSv1.2.0"
        The thermodynamic model to use (e.g., "MELTSv1.2.0" or "pMELTS").
    bulk : pd.DataFrame
        A DataFrame containing the oxide compositions (wt%) and volatile contents 
        ($\text{H}_2\text{O}$ and $\text{CO}_2$) required for the calculations.
    T_fixed_C : float or np.ndarray
        The fixed temperature(s) in Celsius at which to perform the $P_{\text{sat}}$ search. **Required.**
    P_bar_init : float or np.ndarray, optional
        Initial pressure guess (in bars) for the saturation search. Default is $2000$ bar.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or np.ndarray, optional
        Overrides for initial liquid redox state and volatile contents (wt%). These are applied 
        to the input `bulk` composition.
    fO2_buffer : str, optional
        Oxygen fugacity buffer to constrain oxidation states (e.g., "FMQ" or "NNO").
    fO2_offset : float or np.ndarray, optional
        Offset from the specified $\text{f}\text{O}_2$ buffer (in log units).
    copy_columns : list of str, optional
        Placeholder argument (currently not explicitly used in the exposed code).
    quiet_engine : bool, default True
        If True, redirects the MELTS engine's stdout/stderr to the null device for the duration of
        each calculation, which can give a significant speed-up by avoiding large volumes of text
        being written to the console/notebook. Set to False to inspect the raw engine output for
        debugging purposes.

    Returns:
    ----------
    Results : pd.DataFrame
        DataFrame containing the calculated saturation pressure ($P_{\text{sat}}$), the
        corresponding melt composition, and the initial system parameters.
    
    Raises:
    -------
    Warning
        - If $\text{H}_2\text{O}$ is zero in a MELTS model (may affect oxide/apatite saturation).
        - If $\text{Fe}^{3+}/\\Sigma\text{Fe}$ is zero without an $\text{f}\text{O}_2$ buffer in a MELTS model (may cause calculation failure).
    '''
    warnings.warn("findSatPressure_multi() is deprecated and will be removed in a future release; use "
                  "saturation_pressure() instead.", DeprecationWarning, stacklevel=2)
    cores = resolve_cores(cores)
    comp = bulk.copy()

    # if T_fixed_C is None:
    #     raise Warning("Please specify a temperature.")

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("Adding small amounts of H$_{2}$O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H$_{2}$O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    if P_bar_init is None:
        P_bar_init = np.zeros(len(comp['SiO2_Liq'])) + 2000

    comp.loc[:, "Sample_ID_Liq"] = np.linspace(0,len(comp['SiO2_Liq'])-1, len(comp['SiO2_Liq']))
    
    if type(P_bar_init) == int or type(P_bar_init) == float:
        P_bar_init = np.zeros(len(comp['SiO2_Liq'])) + P_bar_init
    comp.loc[:, "P_bar"] = P_bar_init

    if type(T_fixed_C) == int or type(T_fixed_C) == float:
        T_fixed_C = np.zeros(len(comp['SiO2_Liq'])) + T_fixed_C
    comp.loc[:, "T_fixed_C"] = T_fixed_C

    if fO2_offset is None:
        comp.loc[:, "fO2_offset"] = [None] * int(len(comp['SiO2_Liq']))
    else:
        if type(fO2_offset) == float or type(fO2_offset) == int:
            fO2_offset = np.zeros(len(comp['SiO2_Liq'])) + fO2_offset
        comp.loc[:, "fO2_offset"] = fO2_offset    
    
    splits = np.array_split(comp, cores)
    Combined = None

    qs = []
    q = Queue()

    ps = []
    for i in range(cores):
        df = splits[i].reset_index(drop=True)
        p = Process(target = satP_multi, args = (q, i),
                        kwargs = {'Model': Model, 'comp': df,
                                'fO2_buffer': fO2_buffer, 'quiet_engine': quiet_engine})

        ps.append(p)
        p.start()

    TIMEOUT = 20
        
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

    for i in range(len(qs)):
        if len(qs[i]) > 0:
            Res, index = qs[i]
            if Combined is None:
                Combined = Res.copy()
            else:
                Combined = pd.concat([Combined, Res])

    Results = Combined.sort_values(by='Sample_ID_Liq')

    # Reset index if needed
    Results.reset_index(drop=True, inplace=True)

    return Results

def findSatPressure(cores = None, Model = None, bulk = None, T_C_init = None, T_fixed_C = None,
                    P_bar_init = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None,
                    fO2_buffer = None, fO2_offset = None, copy_columns = None, quiet_engine = True):
    """
    Calculates the saturation pressure of a specified composition in the liquid or melt phase. The function will return the saturation pressure of the liquid or melt as a pandas dataframe. If the saturation pressure cannot be calculated, the function will return an empty dataframe.

    Parameters:
    cores (int): The number of CPU cores to use for parallel processing. If left as None, all available CPU cores will be used.
    Model (str): The geochemical model to be used. Options include "MELTSv1.2.0" (default) or "MELTSv1.1.0".
    bulk (dict or pandas.core.frame.DataFrame): A dictionary or pandas dataframe containing the starting composition of the liquid or melt. Elements should be specified in weight percent.
    T_C_init (float, int, or np.ndarray): The starting temperature for the calculation in degrees Celsius. If left as None, the default value of 1200 C will be used.
    P_bar_init (float, int, or np.ndarray): The starting pressure for the calculation in bar. If left as None, the default value of 10000 bar will be used.
    Fe3Fet_Liq (float): The initial concentration of ferric iron in the liquid or melt in weight percent. If left as None, the value will be taken from the bulk composition.
    H2O_Liq (float): The initial concentration of water in the liquid or melt in weight percent. If left as None, the value will be taken from the bulk composition.
    fO2_buffer (float): The  oxygen fugacity buffer to be used in the calculation. If left as None, no buffer will be applied
    fO2_offset (float, int, or np.ndarray): The oxygen fugacity offset (log units) from the specified buffer to be applied to the calculation. If left as None, no offset will be applied. If a list is provided, each offset will be applied to the corresponding element in the bulk composition.
    quiet_engine (bool, optional, default True): If True, redirects the MELTS engine's stdout/stderr to the null device for the duration of each calculation, which can give a significant speed-up by avoiding large volumes of text being written to the console/notebook. Set to False to inspect the raw engine output for debugging purposes.

    Returns:
    Res (dict or pandas.core.frame.DataFrame): A dict or dataframe containing the chemical composition of the melt phase and the saturation pressure and temperature of the melt.

    Examples:

    Calculate the saturation pressure and temperature of a basaltic liquid, with initial conditions set at 1200 C and 10000 bar using all available CPU cores and the MELTSv1.2.0 model.
    >>> findSatPressure(bulk = {'SiO2_Liq': 50, 'TiO2_Liq': 1, 'Al2O3_Liq': 15, 'FeOt_Liq': 10, 'MnO_Liq': 0.5, 'MgO_Liq': 8, 'CaO_Liq': 10, 'Na2O_Liq': 2, 'K2O_Liq': 1, 'P2O5_Liq': 0.5, 'H2O_Liq': 0.5, 'Fe3Fet_Liq': 0.15}, T_C_init = 1200, P_bar_init = 10000)

    """
    warnings.warn("findSatPressure() is deprecated and will be removed in a future release; use "
                  "saturation_pressure() instead.", DeprecationWarning, stacklevel=2)

    # set default values if required
    if Model is None:
        Model == "MELTSv1.2.0"

    comp = bulk.copy()

    if T_C_init is None and T_fixed_C is None:
        T_C_init = 1200

    if P_bar_init is None:
        P_bar_init = 10000

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("Adding small amounts of H$_{2}$O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H$_{2}$O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    cores = resolve_cores(cores)

    qs = []
    q = Queue()

    # specify the number of calculations to be performed in each sequence
    if type(comp) == dict:
        p = Process(target = satP, args = (q, 0),
                    kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                              'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset,
                            'quiet_engine': quiet_engine})

        p.start()
        try:
            ret = q.get(timeout = 60)
        except:
            ret = []

        TIMEOUT = 5
        start = time.time()
        if p.is_alive():
            while time.time() <= TIMEOUT:
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

        if len(ret) > 0:
            Results, index = ret

        Res = Results.copy()

        return Res

    else:
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq'])%cores

        if A > 0:
            Group = np.zeros(A) + cores
            if B > 0:
                Group = np.append(Group, B)
        else:
            Group = np.array([B])

        L = np.sum(Group)
        if type(P_bar_init) == float or type(P_bar_init) == int:
            P_bar_init = np.zeros(int(L)) + P_bar_init
        if T_fixed_C is None:
            if type(T_C_init) == float or type(T_C_init) == int:
                T_C_init = np.zeros(int(L)) + T_C_init
            T_fixed_C = [None]*int(L)
        else:
            if type(T_fixed_C) == float or type(T_fixed_C) == int:
                T_fixed_C = np.zeros(int(L)) + T_fixed_C
            T_C_init = [None] * int(L)
        if fO2_offset is None:
            fO2_offset = [None] * int(L)
        else:
            if type(fO2_offset) == float or type(fO2_offset) == int:
                fO2_offset = np.zeros(int(L)) + fO2_offset

        for j in tqdm(range(len(Group))):
            ps = []

            for i in range(int(cores*j), int(cores*j + Group[j])):
                p = Process(target = satP, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'T_C_init': T_C_init[i],
                                      'T_fixed_C': T_fixed_C[i], 'P_bar_init': P_bar_init[i],
                                    'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i],
                                    'quiet_engine': quiet_engine})

                ps.append(p)
                p.start()

            TIMEOUT = 20
            
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

        Res = pd.DataFrame(data = np.zeros((int(L),15)), columns = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq', 'P_bar', 'T_Liq_C'])
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Results, index = qs[i]
                Res.loc[index] = Results

        if copy_columns is not None:
            if type(copy_columns) == str:
                Res.insert(0, copy_columns, comp[copy_columns])
                # Af_Combined.insert(0, copy_columns, comp[copy_columns])
            elif type(copy_columns) == list:
                j = 0
                for i in copy_columns:
                    Res.insert(j, i, comp[i])
                    # Af_Combined.insert(j, i, comp[i])
                    j = j + 1

        return Res


def saturation_pressure_v1(Model = "MELTSv1.2.0", cores = MAX_WORKERS,
                        bulk = None, T_init_C = None, T_fixed_C = None, P_init_bar = 5000,
                        Fe3Fet_Liq = None, Fe3Fet_init = None, H2O_Liq = None, H2O_init = None,
                        CO2_Liq = None, CO2_init = None, fO2_buffer = None, fO2_offset = None,
                        copy_columns = None, multi_processing = True, timeout = 20, max_retries = 2,
                        quiet_engine = True):
    """
    Estimates the pressure at which volatile saturation occurs for one or more melt compositions
    using MELTS thermodynamic models.

    Supports both single and batched calculations, with optional multiprocessing to accelerate large-scale runs.

    Parameters
    ----------
    Model : str, default 'MELTSv1.2.0'
        The thermodynamic model to use. 
            MELTS:
                "MELTSv1.2.0"
                "MELTSv1.1.0"
                "MELTSv1.0.2" - only H2O
                "pMELTS" - only H2O
    cores : int, default petthermotools.core_config.MAX_WORKERS
        Number of parallel processes to use.
    bulk : dict, pandas.Series, or pandas.DataFrame
        Bulk melt composition(s). If a Series or single dict is passed, a single calculation is run.
        If a DataFrame is passed, a calculation is performed for each row.
    T_init_C : float or np.ndarray, optional
        Initial temperature(s) in degrees Celsius. System will co-solve the liquidus temperature and pressure of volatile saturation. Used if not specifying `T_fixed_C`.
    T_fixed_C : float or np.ndarray, optional
        Fixed temperature(s) for the saturation pressure calculation. Solid phases are disabled in all calculations.
    P_init_bar : float or np.ndarray, default 5000
        Initial pressure guess(es) in bars.
    Fe3Fet_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial Fe³⁺/∑Fe ratio. Prefer `Fe3Fet_init`.
    Fe3Fet_init : float or np.ndarray, optional
        Initial Fe³⁺/∑Fe ratio(s) in the liquid.
    H2O_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial H₂O content. Prefer `H2O_init`.
    H2O_init : float or np.ndarray, optional
        Initial H₂O content(s) in the liquid, in wt%.
    CO2_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial CO₂ content. Prefer `CO2_init`.
    CO2_init : float or np.ndarray, optional
        Initial CO₂ content(s) in the liquid, in wt%.
    fO2_buffer : str, optional
        Oxygen fugacity buffer to use ('FMQ' or 'NNO').
    fO2_offset : float or np.ndarray, optional
        Offset from the fO₂ buffer in log units (e.g., +1.0 = 1 log unit above buffer).
    copy_columns : str or list of str, optional
        If bulk is a DataFrame, copies specified columns from input into output DataFrame for tracking.
    multi_processing : bool, default True
        Whether to parallelize the calculations across multiple cores.
    timeout : int, default 180
        Timeout for individual parallel processes, in seconds.
    quiet_engine : bool, default True
        If True, redirects the MELTS engine's stdout/stderr to the null device for the duration of
        each calculation, which can give a significant speed-up by avoiding large volumes of text
        being written to the console/notebook. Set to False to inspect the raw engine output for
        debugging purposes.

    Returns
    -------
    pandas.DataFrame or dict
        - If multiple calculations are performed, returns a DataFrame with saturation pressures and
          other outputs, indexed by run number.
        - If a single calculation is run (`multi_processing=False` or `bulk` is a single composition),
          returns a dictionary of calculated results.

    Raises
    ------
    Warning
        - If required MELTS packages are not found in the Python path.
        - If unsupported fO₂ buffers are specified.
        - If zero Fe³⁺/∑Fe or H₂O is specified in MELTS without buffer constraints.

    Notes
    -----
    - Missing inputs for temperature, pressure, or fO₂ offset are auto-expanded if multiple runs are inferred.
    - For MELTS-based models, small amounts of H₂O and ferric iron often improve model stability.
    - The function handles multiprocessing with safe timeouts and partial result recovery.
    """
    cores = resolve_cores(cores)

    # converting back to old naming convention, will eventually be fully replaced.
    T_C_init = T_init_C
    P_bar_init = P_init_bar
       

    if bulk is not None:
        comp = bulk.copy()

    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.1.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.1.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.1.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

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
    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq = CO2_init)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0 and "MELTS" in Model:
            raise Warning("Adding small amounts of H$_{2}$O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H$_{2}$O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and "MELTS" in Model and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    # specify the number of calculations to be performed in each sequence
    One = 0
    if type(comp) == pd.core.frame.DataFrame: # simplest scenario - one calculation per bulk composition imported
        A = len(comp['SiO2_Liq'])//cores
        B = len(comp['SiO2_Liq'])%cores
    else:
        if T_fixed_C is not None and type(T_fixed_C) == np.ndarray: # one calculation per T loaded.
            A = len(T_fixed_C)//cores
            B = len(T_fixed_C)%cores
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
    if One == 1:
        if multi_processing:
            p = Process(target = satP, args = (q,0),
                        kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                                  'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                                  'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset,
                                  'quiet_engine': quiet_engine})
            p.start()
            try:
                ret = q.get(timeout = timeout)
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
            
            if len(ret) > 0:
                Results, index = ret
            else:
                Results = {}

            return Results
        else:
            with (_quiet_engine_output() if quiet_engine else contextlib.nullcontext()):
                Results = findSatPressure_MELTS(Model = Model, T_C_init=T_C_init, T_fixed_C=T_fixed_C,
                                                P_bar_init=P_bar_init, comp = comp, fO2_buffer=fO2_buffer,
                                                fO2_offset=fO2_offset)
            return Results
    else: # perform multiple crystallisation calculations
        # first make sure everything is the right length
        L = np.sum(Group)
        if type(P_bar_init) == float or type(P_bar_init) == int:
            P_bar_init = np.zeros(int(L)) + P_bar_init

        if T_C_init is None:
            T_C_init = [None] * int(L)
        elif type(T_C_init) == float or type(T_C_init) == int:
            T_C_init = np.zeros(int(L)) + T_C_init

        if T_fixed_C is None:
            T_fixed_C = [None] * int(L)
        elif type(T_fixed_C) == float or type(T_fixed_C) == int:
            T_fixed_C = np.zeros(int(L)) + T_fixed_C

        if fO2_offset is None:
            fO2_offset = [None] * int(L)
        elif type(fO2_offset) == float or type(fO2_offset) == int:
            fO2_offset = np.zeros(int(L)) + fO2_offset

        index_in = np.arange(int(L))
        combined_results = {}
        index_out = np.array([], dtype=int)
        timeout_main = timeout

        while len(index_out)<len(index_in):
            index = np.setdiff1d(index_in, index_out)
            groups = np.array_split(index, cores)
            non_empty_groups = [g for g in groups if g.size > 0]
            groups = non_empty_groups

            if len(groups[0])< 10:
                timeout = len(groups[0])*timeout_main
            else:
                timeout = 10*timeout_main
                for g in range(len(groups)):
                    groups[g] = groups[g][:10]
            
            processes = []
            Start = time.time()
            for i in range(len(groups)):
                q = Queue()
                p = Process(target = saturationP_multi, args = (q, groups[i]),
                            kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                                      'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                                      'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset,
                                      'quiet_engine': quiet_engine})
                
                processes.append([p, q, groups[i]])
                p.start()

            for p, q, group in processes:
                try:
                    if time.time() - Start < timeout + 5:
                        res = q.get(timeout = (timeout - (time.time()-Start) + 5))
                    else:
                        res = q.get(timeout=1)
                except:
                    res = []
                    idx_chunks = np.array([group[0]], dtype = int)

                p.join(timeout=1)
                p.terminate()

                if len(res) > 0.0:
                    idx_chunks, results = res
                    combined_results.update(results)

                index_out = np.hstack((index_out, idx_chunks))
            
            print(f"Completed {100*len(index_out)/len(index_in)} %")

        results = combined_results

        def results_to_dataframe(results, full_index):
            # Convert 'Run X' → X and create DataFrame
            df = pd.DataFrame.from_dict(results, orient='index')
            df.index = df.index.str.extract(r'Run (\d+)')[0].values.astype(int)

            # Reindex with full range, fill missing with 0.0
            df = df.reindex(full_index, fill_value=np.nan)

            # Optional: sort index if needed
            df = df.sort_index()

            return df
        
        Results = results_to_dataframe(results, index_in)
        
        # --- RETRY LOGIC ---
        # We check if there are any NaNs in the saturation pressure column
        # Assuming the column name is 'P_bar' (adjust if findSatPressure_MELTS uses another name)
        sat_col = 'P_bar' 
        
        if max_retries > 0:
            for retry_count in range(max_retries):
                failed_mask = Results[sat_col].isna()
                failed_indices = Results.index[failed_mask].tolist()

                if T_C_init[0] is not None:
                    if retry_count == 0:
                        T_C_init = T_C_init + 25
                    elif retry_count == 1:
                        T_C_init = T_C_init - 15
                
                if retry_count % 2 == 0:
                    P_bar_init = P_bar_init * 0.7
                else:
                    P_bar_init = P_bar_init * 2

                if len(failed_indices) == 0:
                    break
                    
                print(f"Attempting to resolve {len(failed_indices)} failed calculations (Retry {retry_count + 1})...")

                # comp = bulk.copy()
                # Re-prepare inputs for only the failed indices
                # We increase the timeout slightly for retries to handle edge-case convergence
                retry_timeout = timeout * 1.5 
                

                # Use the same multi-processing logic but on the subset
                # (You can encapsulate your while loop into a helper to keep this clean)
                retry_combined_results = {}
                index_in_retry = np.array(failed_indices)
                index_out_retry = np.array([], dtype=int)

                while len(index_out_retry) < len(index_in_retry):
                    curr_idx = np.setdiff1d(index_in_retry, index_out_retry)
                    groups = np.array_split(curr_idx, cores)
                    groups = [g for g in groups if g.size > 0]

                    processes = []
                    for i in range(len(groups)):
                        q = Queue()
                        p = Process(target=saturationP_multi, args=(q, groups[i]),
                                    kwargs={'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                                            'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset,
                                            'quiet_engine': quiet_engine})
                        processes.append([p, q, groups[i]])
                        p.start()

                    for p, q, group in processes:
                        try:
                            res = q.get(timeout=retry_timeout)
                            if len(res) > 0:
                                idx_chunks, results = res
                                retry_combined_results.update(results)
                                index_out_retry = np.hstack((index_out_retry, idx_chunks))
                        except:
                            index_out_retry = np.hstack((index_out_retry, [group[0]])) # Mark as attempted
                        p.join(timeout=2)
                        p.terminate()
                
                if len(retry_combined_results)>0:
                    # Convert retry results to DF and update the main Results DF
                    retry_df = results_to_dataframe(retry_combined_results, failed_indices)
                    Results.update(retry_df)

        if copy_columns is not None:
            if type(copy_columns) == str:
                Results.insert(0, copy_columns, bulk[copy_columns])
                # Af_Combined.insert(0, copy_columns, comp[copy_columns])
            elif type(copy_columns) == list:
                j = 0
                for i in copy_columns:
                    Results.insert(j, i, bulk[i])
                    # Af_Combined.insert(j, i, comp[i])
                    j = j + 1

        return Results


class _SaturationRowTask(NamedTuple):
    """
    One row's own data for saturation_pressure_v2, carried through
    mp_utils.run_tasks's `indices` -- same rationale as findLiq_multi's
    _RowTask / equilibrate_multi's _EquilibrateRowTask (see either one's own
    docstring): a few small scalars and one composition dict, not an
    O(n)-sized structure passed through task_kwargs, so per-task spawn cost
    does not grow with batch size.

    Carries BOTH `T_C_init` and `T_fixed_C` per row (only one of which is
    normally non-None for a given row) because findSatPressure_MELTS uses
    their presence/absence to pick between its two modes -- see
    saturation_pressure_v2's own docstring for what those modes mean.
    """
    row_id: int
    P_bar_init: float
    T_C_init: float
    T_fixed_C: float
    comp: dict
    fO2_offset: float

    def __hash__(self):
        # see _RowTask.__hash__'s docstring -- same reason (comp is an
        # unhashable dict; row_id alone is a unique identity per task).
        return hash(self.row_id)

    def __eq__(self, other):
        return isinstance(other, _SaturationRowTask) and self.row_id == other.row_id


_cached_sat_melts = None         # per-worker-process cache of the last SUCCESSFULLY-used
_cached_sat_melts_model = None   # MELTSdynamic engine, see _saturation_task's own docstring.
_cached_sat_probed = False       # mode 1 only: has base_melts already run its ONE-TIME
                                 # probe-and-suppress step for this cached engine's lifetime?


def _saturation_task(q, task, *, Model = None, fO2_buffer = None, quiet_engine = True):
    """
    task_fn for saturation_pressure_v2's mp_utils.run_tasks() call. Supports
    cross-task engine reuse (when max_tasks_per_worker > 1), with a MINIMAL
    addNodeAfter() hop per task -- at most one, same spirit as
    _findLiq_task/_equilibrate_task's own reuse, but EXPLICITLY NOT the same
    "always safe" guarantee: this session's investigation (MAINTENANCE_NOTES.md,
    2026-10-14) found that reusing one MELTSdynamic across DIFFERENT,
    independent compositions via addNodeAfter() carries a real, measured risk
    of a crash or a silent wrong answer for a general, multi-phase (solid-rich)
    system. This function is where that's being tried for
    saturation_pressure_v2 specifically BECAUSE mode 1 (T_fixed_C) suppresses
    every solid phase, leaving only liquid + fluid -- a much simpler system
    than equilibrate_multi's general case, with much less phase-assemblage
    state for a reused engine to carry over incorrectly between rows. Still
    needs its own direct verification before being trusted as a default --
    see saturation_pressure_v2's own docstring for the current state of that.

    Engine reuse here follows the SAME convention as _equilibrate_task: cache
    the object reference actually passed to findSatPressure_MELTS/
    findLiq_MELTS (not whatever node either of them may have internally
    reassigned their own local `melts` to -- neither function returns the
    engine object, so there is nothing else to cache), and hand that SAME
    reference back in for the next task, which addNodeAfter()s from it again.

    The hop itself is placed differently per mode, because findSatPressure_MELTS
    only does its OWN internal addNodeAfter() for mode 1:
    - Mode 1 (T_fixed_C given): findSatPressure_MELTS already does exactly one
      `melts = melts.addNodeAfter()` right after its probe-and-suppress step
      -- so the cached engine is passed straight through un-hopped, and that
      existing internal hop is the only one that happens.
      `suppressed` is True for every reused call AFTER the first one this
      cached engine ever sees (tracked by `_cached_sat_probed`, reset
      whenever a fresh engine is constructed), and left at its default (None,
      letting findSatPressure_MELTS run its own probe) only for that first
      call. This mirrors the OLD saturationP_multi's own design exactly
      (confirmed by re-reading it, MAINTENANCE_NOTES.md 2026-10-15): it
      probes-and-suppresses ONCE per spawned process, before its per-item
      loop, then passes suppressed=True for every item in that process
      afterwards -- it never re-probes per row either.

      An EARLIER version of this function instead left `suppressed=None` on
      every call, reasoning the extra suppressions were harmless ("can only
      suppress more solids, never fewer"). That was wrong on two counts, both
      confirmed live against the user's own 149-row T_fixed_C batch
      (MAINTENANCE_NOTES.md 2026-10-15): (1) SPEED -- the probe is itself a
      real equilibrium-adjacent solve (`calcSaturationState()`), not a cheap
      lookup, re-paid on every single reused row instead of once per worker;
      (2) CORRECTNESS -- the probe's `setBulkComposition()` +
      `calcSaturationState()` + `setSystemProperties("Suppress", ...)` calls
      run DIRECTLY on `base_melts` itself (it is never hopped before the
      probe), so repeating them every row means mutating ONE shared ancestor
      node with a DIFFERENT row's composition over and over, across its
      entire reused lifetime -- exactly the "reuse across independent
      compositions" pattern this same session's broader equilibrate_multi
      investigation already found causes crashes and silent wrong answers,
      independently rediscovered here. Probing only once freezes
      `base_melts`'s own engine state after that first call: every later row
      only ever adds a FRESH `addNodeAfter()` child (with its own
      `setBulkComposition()`) off that same, now-immutable ancestor, which is
      the same "star" pattern already proven safe elsewhere in this file.

      Known, accepted limitation (same one OLD code already had, not a new
      risk introduced by matching it): the suppress list for a given cached
      engine's whole lifetime is fixed by whichever row that engine happened
      to probe FIRST. A later row in the same worker whose mineralogy would
      want to suppress a different/additional phase near T_fixed_C+500 does
      not get that chance. This is strictly better than old code's own
      version of the same limitation (which hard-coded row 0 of the WHOLE
      input table regardless of which rows a given worker actually got), but
      is not perfect per-row suppression. In T_fixed_C mode specifically this
      is normally moot: at T_fixed_C+500 the system is usually pure liquid +
      fluid with nothing left to suppress at all (the probe's loop finds
      `PL` containing only "fluid"/"water"), which is exactly why reuse is
      being tried here and not for the general multi-phase case.
    - Mode 2 (T_C_init given, T_fixed_C left None): findSatPressure_MELTS has
      no such internal hop before its findLiq_MELTS calls, so THIS function
      does the one hop itself, only when reusing a cached engine (a freshly
      constructed one needs no hop, same as elsewhere in this codebase).
      `suppressed` has no effect in this mode (confirmed by reading
      findSatPressure_MELTS: the `suppressed`-aware block is inside the
      `if T_fixed_C is not None:` branch only), so it is always left at its
      default here.

    findSatPressure_MELTS's own failure sentinel is a dict with every value
    zeroed (P_bar=0.0, T_Liq_C=0.0, ...), not an empty dict/exception -- this
    wrapper raises explicitly on that so mp_utils.run_tasks sees a genuine
    CRASHED status, mirroring _findLiq_task's/_equilibrate_task's own
    "if not Results: raise" pattern. Only a melts reference that just
    produced a genuine success is cached, for the same reason _findLiq_task's
    own docstring gives: _worker_loop retires the whole process on ANY
    exception, so this is never read after a failure.
    """
    global _cached_sat_melts, _cached_sat_melts_model, _cached_sat_probed

    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's docstring for why --
        # must run before MELTSdynamic() is constructed below, same reasoning as
        # _findLiq_task/_equilibrate_task.

    reused = _cached_sat_melts is not None and _cached_sat_melts_model == Model
    if reused:
        base_melts = _cached_sat_melts
    else:
        from meltsdynamic import MELTSdynamic
        base_melts = MELTSdynamic(_MODEL_CODES.get(Model, 1))
        base_melts.engine.setSystemProperties("Output", "none")   # old saturationP_multi and
        # _findLiq_task_v2's trunk both set this once per fresh engine; this task_fn never did --
        # confirmed missing by reading the code directly (grepped for "Output" in this file and
        # found zero matches inside this function), not previously noticed. Without it the engine
        # writes its full, verbose *_tbl.txt / liquid-model-batch.inp output to the CURRENT WORKING
        # DIRECTORY on (seemingly) every calcEquilibriumState() call -- with cores>1 ALL concurrent
        # workers share that directory, so they're writing to/overwriting the SAME filenames at
        # once, inside a Drive-synced folder besides. Mode 2's nested P-bisection + per-iteration
        # findLiq_MELTS T-search makes far more calcEquilibriumState() calls per row than mode 1's
        # simple, all-solids-suppressed P-only bisection, so this gap is suspected (not yet proven)
        # to disproportionately hurt mode 2 -- see MAINTENANCE_NOTES.md 2026-10-16.
        _cached_sat_probed = False   # fresh engine: its first mode-1 call must probe

    # Mode 2 (T_fixed_C is None) needs its one hop done explicitly here, every time,
    # from the SAME stable base_melts reference -- not chained from wherever the
    # previous task's own hop ended up. Mode 1 passes base_melts straight through
    # un-hopped: findSatPressure_MELTS does its own single internal addNodeAfter()
    # unconditionally for that mode, also always starting from this same reference
    # (it never gets reassigned here), which is what keeps this a "star" pattern
    # for both modes -- same convention _equilibrate_task uses, see its own
    # docstring for why that, not chaining, is what's been tested as safe.
    call_melts = base_melts.addNodeAfter() if task.T_fixed_C is None else base_melts

    # suppressed=True skips findSatPressure_MELTS's own probe-and-suppress step
    # entirely for every mode-1 call after the first one this cached engine ever
    # sees -- see the docstring above for why that is both the fast AND the
    # correct choice, matching old saturationP_multi's own design.
    suppressed = True if (task.T_fixed_C is not None and _cached_sat_probed) else None

    out, tr = findSatPressure_MELTS(Model = Model, T_C_init = task.T_C_init, T_fixed_C = task.T_fixed_C,
                                     P_bar_init = task.P_bar_init, comp = task.comp, fO2_buffer = fO2_buffer,
                                     fO2_offset = task.fO2_offset, trial = True, melts = call_melts,
                                     suppressed = suppressed)

    if out['P_bar'] == 0.0 and out['T_Liq_C'] == 0.0:
        raise RuntimeError(f"findSatPressure_MELTS did not converge at P_bar_init={task.P_bar_init}, "
                            f"T_C_init={task.T_C_init}, T_fixed_C={task.T_fixed_C} (zeroed result)")

    _cached_sat_melts = base_melts
    _cached_sat_melts_model = Model
    if task.T_fixed_C is not None:
        _cached_sat_probed = True

    q.put([out, task.row_id])


def saturation_pressure(Model = "MELTSv1.2.0", cores = MAX_WORKERS,
                        bulk = None, T_init_C = None, T_fixed_C = None, P_init_bar = 5000,
                        Fe3Fet_Liq = None, Fe3Fet_init = None, H2O_Liq = None, H2O_init = None,
                        CO2_Liq = None, CO2_init = None, fO2_buffer = None, fO2_offset = None,
                        copy_columns = None, timeout = None, max_domain_retries = 1,
                        retry_timeout_multiplier = 2.0, max_tasks_per_worker = 20, quiet_engine = True):
    """
    2026-10-21: this IS the v2 implementation, now the default under the
    `saturation_pressure` name (previously available as
    `saturation_pressure_v2`, which remains as a plain alias). UNLIKE
    multi_path/equilibrate_multi/findLiq_multi, reuse is ON by default here
    (max_tasks_per_worker=20, not 1) -- a deliberate choice: mode 1
    (T_fixed_C) suppresses every solid phase, leaving far less phase-
    assemblage state for a reused engine to carry over incorrectly between
    independent compositions than equilibrate_multi's general case, and
    the cross-row hop here is a single, minimal addNodeAfter() per task,
    not the deeper chains tried and abandoned for equilibrate_multi -- see
    the reuse discussion a few paragraphs below, which predates this
    rename and still applies verbatim: this has NOT had the same scale of
    adversarial testing as equilibrate_multi's own reuse investigation, so
    treat the specific default (20) as a starting point to revisit, not a
    settled number; the old, always-reuse-free implementation is kept,
    fully intact, as `saturation_pressure_v1` for anyone who wants it (or
    pass `max_tasks_per_worker=1` here for the same effect without
    switching functions).

    Alternative to saturation_pressure_v1() using the mp_utils.py persistent-worker
    supervisor (the same one findLiq_multi/equilibrate_multi use), instead of the
    original's per-round, fresh-process-per-group-per-round scheme. Same two modes,
    chosen the same way (by whether T_fixed_C or T_init_C is supplied):

    1. T_fixed_C given -- fixed temperature. All solid phases are suppressed (via
       findSatPressure_MELTS's own one-time "probe at T+500C, suppress everything
       except fluid/water" step), and the function searches only for the pressure
       at which a fluid phase is saturated in the melt at that fixed T.
    2. T_init_C given (T_fixed_C left as None) -- unknown temperature. No phases
       are suppressed; the liquidus and the fluid saturation pressure are solved
       for TOGETHER (findLiq_MELTS is called at each trial pressure, so every
       step of the pressure search re-solves the liquidus at that pressure too,
       rather than assuming a fixed T).

    This is a straight engine-supervisor swap -- it does not change
    findSatPressure_MELTS's own search algorithm (same bracket-then-bisect in
    either mode) or its outputs, just how rows are distributed across worker
    processes. Returns the same output shape as saturation_pressure(): a dict for
    a single composition, a DataFrame (one row per input) for a batch.

    UNLIKE equilibrate_multi's final posture, this DOES support cross-row engine
    reuse when max_tasks_per_worker > 1 (default 1, i.e. still no reuse unless
    asked for) -- with a MINIMAL addNodeAfter() hop per task (at most one), not
    the deeper chains tried and abandoned for equilibrate_multi. This is
    deliberately being tried here specifically because mode 1 (T_fixed_C) is a
    much simpler system than equilibrate_multi's general case: every solid phase
    is suppressed, leaving only liquid + fluid, so there is far less
    phase-assemblage state for a reused engine to carry over incorrectly between
    independent compositions. See _saturation_task's own docstring for exactly
    how the hop is placed for each mode, and run the correctness checks described
    there (identical-row and varied-row comparisons against max_tasks_per_worker=1)
    before relying on reuse for anything that matters -- this has NOT yet had the
    same scale of adversarial testing that led equilibrate_multi to abandon reuse
    by default, so treat it as an experiment to validate per-use-case, not a
    proven-safe default.

    A failed row is retried at most max_domain_retries times, for ANY failure --
    unlike equilibrate_multi/findLiq_multi, a genuine calculation error (the
    zeroed-result sentinel) IS retried here too, not just a timeout, matching the
    original saturation_pressure()'s own retry loop (which can't even tell a
    timeout apart from a calculation error, and retries anything still
    unconverged regardless). Every retry -- whether triggered by a timeout or a
    calculation error -- perturbs the starting guess before trying again, same
    spirit as the original's own retry loop: P_bar_init is alternated by x0.7/x2
    each round, and (only for a row using mode 2, i.e. T_C_init was actually
    given) T_C_init is nudged by +25 C then -15 C, compounding across rounds.
    T_fixed_C is NEVER perturbed -- mode 1 requires it stay fixed by definition.
    This is workable here, unlike for a bare equilibrate call, because
    findSatPressure_MELTS's own bracket+bisect search is genuinely sensitive to
    its starting P/T guess -- a calculation error from a bad bracket is exactly
    the kind of failure a different starting point can plausibly fix, so
    retrying blindly isn't just paying the same failure again.

    Parameters
    ----------
    Model : str, default 'MELTSv1.2.0'
        The thermodynamic model to use.
            MELTS:
                "MELTSv1.2.0"
                "MELTSv1.1.0"
                "MELTSv1.0.2" - only H2O
                "pMELTS" - only H2O
    cores : int, default petthermotools.core_config.MAX_WORKERS
        Number of parallel worker processes to use.
    bulk : dict, pandas.Series, or pandas.DataFrame
        Bulk melt composition(s). If a Series or single dict is passed, a single
        calculation is run. If a DataFrame is passed, a calculation is performed
        for each row.
    T_init_C : float or np.ndarray, optional
        Initial temperature(s) in degrees Celsius for mode 2 (unknown T -- the
        liquidus and saturation pressure are co-solved). Ignored for a row where
        `T_fixed_C` is also given.
    T_fixed_C : float or np.ndarray, optional
        Fixed temperature(s) in degrees Celsius for mode 1 (solid phases
        suppressed, only the saturation pressure is solved for).
    P_init_bar : float or np.ndarray, default 5000
        Initial pressure guess(es) in bars.
    Fe3Fet_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial Fe3+/sum(Fe) ratio. Prefer `Fe3Fet_init`.
    Fe3Fet_init : float or np.ndarray, optional
        Initial Fe3+/sum(Fe) ratio(s) in the liquid.
    H2O_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial H2O content. Prefer `H2O_init`.
    H2O_init : float or np.ndarray, optional
        Initial H2O content(s) in the liquid, in wt%.
    CO2_Liq : float or np.ndarray, optional [DEPRECATED]
        Legacy alias for initial CO2 content. Prefer `CO2_init`.
    CO2_init : float or np.ndarray, optional
        Initial CO2 content(s) in the liquid, in wt%.
    fO2_buffer : str, optional
        Oxygen fugacity buffer ('FMQ', 'NNO' or 'IW').
    fO2_offset : float or np.ndarray, optional
        Offset from the fO2 buffer in log units.
    copy_columns : str or list of str, optional
        If bulk is a DataFrame, copies specified columns from the input into the
        output DataFrame for tracking.
    timeout : float or None, default None
        Per-calculation budget in seconds for the first attempt, before a stuck
        worker is killed and (if it failed by timing out, not by a calculation
        error) retried once at `retry_timeout_multiplier` times this value. None
        (the default) resolves from the platform, matching findLiq_multi/
        equilibrate_multi: 20s on macOS/Linux, 45s on Windows -- this search does
        several bracket-expansion and bisection steps per row (more work than a
        single equilibrate call), so raise this if you see rows timing out on
        genuinely slow compositions rather than assume they're unrecoverable.
    max_domain_retries : int, default 1
        How many extra rounds to attempt for rows that timed out. Every round is
        a single batch (one run_tasks call over everything still pending).
    retry_timeout_multiplier : float, default 2.0
        Each retry round's timeout is the previous round's multiplied by this.
    max_tasks_per_worker : int, default 1
        How many successful calculations one worker process runs before it is
        retired and replaced with a fresh one, reusing the same MELTSdynamic
        engine via a single addNodeAfter() hop in between (see _saturation_task's
        own docstring for exactly how). Defaulted to 1 (no reuse) because this
        has not yet had the same scale of adversarial testing that led
        equilibrate_multi to abandon reuse by default -- raise it deliberately,
        as an experiment, especially for mode 1 (T_fixed_C), where the much
        simpler liquid+fluid-only system is expected to be more tolerant of
        reuse than equilibrate_multi's general case ever was.
    quiet_engine : bool, default True
        If True, the MELTS engine's own progress text is permanently redirected
        to the null device for each worker process.

    Returns
    -------
    pandas.DataFrame or dict
        A dict of saturation-pressure results for a single composition, or a
        DataFrame (one row per input, in the original order; rows that never
        converged are NaN) for a batch.
    """
    if timeout is None:
        timeout = 45.0 if platform.system() == "Windows" else 20.0

    comp = bulk.copy()

    if H2O_Liq is not None:
        print('Warning - the kwarg "H2O_Liq" will be removed from v1.1.0 onwards. Please use "H2O_init" instead.')
        if H2O_init is None:
            H2O_init = H2O_Liq

    if CO2_Liq is not None:
        print('Warning - the kwarg "CO2_Liq" will be removed from v1.1.0 onwards. Please use "CO2_init" instead.')
        if CO2_init is None:
            CO2_init = CO2_Liq

    if Fe3Fet_Liq is not None:
        print('Warning - the kwarg "Fe3Fet_Liq" will be removed from v1.1.0 onwards. Please use "Fe3Fet_init" instead.')
        if Fe3Fet_init is None:
            Fe3Fet_init = Fe3Fet_Liq

    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    if fO2_buffer is not None:
        if fO2_buffer not in ("NNO", "FMQ", "IW"):
            raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet_init, H2O_Liq = H2O_init, CO2_Liq = CO2_init)

    if type(comp) == dict:
        if comp['H2O_Liq'] == 0.0:
            raise Warning("Adding small amounts of H2O may improve the ability of MELTS to accurately reproduce the saturation of oxide minerals. Additionally, sufficient H2O is required in the model for MELTS to predict the crystallisation of apatite, rather than whitlockite.")

        if comp['Fe3Fet_Liq'] == 0.0 and fO2_buffer is None:
            raise Warning("MELTS often fails to produce any results when no ferric Fe is included in the starting composition and an fO2 buffer is not set.")

    # --- normalise n, P_bar_init/T_C_init/T_fixed_C/fO2_offset, same precedence as
    # the original saturation_pressure() (DataFrame > T_fixed_C array > fO2_offset
    # array > a single calculation) ---
    if type(comp) == pd.core.frame.DataFrame:
        n = len(comp['SiO2_Liq'])
    elif T_fixed_C is not None and type(T_fixed_C) == np.ndarray:
        n = len(T_fixed_C)
    elif fO2_offset is not None and type(fO2_offset) == np.ndarray:
        n = len(fO2_offset)
    else:
        n = 1

    if type(P_init_bar) in (float, int):
        P_init_bar = np.zeros(n) + P_init_bar
    if T_init_C is None:
        T_init_C = [None] * n
    elif type(T_init_C) in (float, int):
        T_init_C = np.zeros(n) + T_init_C
    if T_fixed_C is None:
        T_fixed_C = [None] * n
    elif type(T_fixed_C) in (float, int):
        T_fixed_C = np.zeros(n) + T_fixed_C
    if fO2_offset is None:
        fO2_offset = [None] * n
    elif type(fO2_offset) in (float, int):
        fO2_offset = np.zeros(n) + fO2_offset

    if type(comp) == dict:
        comps = [comp] * n
    else:
        comps = [comp.loc[i].to_dict() for i in range(n)]

    cores = resolve_cores(cores)

    pending_tasks = [_SaturationRowTask(i, P_init_bar[i], T_init_C[i], T_fixed_C[i], comps[i], fO2_offset[i])
                      for i in range(n)]
    collected: dict = {}   # row_id -> out dict

    pbar = tqdm(total = n, desc = "saturation_pressure_v2")
    ticked: set = set()

    def _capture(item):
        out, idx = item
        collected[idx] = out
        if idx not in ticked:
            ticked.add(idx)
            pbar.update(1)

    current_timeout = timeout
    try:
        for round_num in range(max_domain_retries + 1):
            if not pending_tasks:
                break

            out = run_tasks(
                _saturation_task,
                pending_tasks,
                task_kwargs = dict(Model = Model, fO2_buffer = fO2_buffer, quiet_engine = quiet_engine),
                cores = cores,
                timeout = current_timeout,
                max_retries = 0,   # the selective (TIMEOUT-only) retry below replaces
                                   # run_tasks's own blanket, status-agnostic one.
                max_tasks_per_worker = max_tasks_per_worker,
                on_progress = _capture,
            )

            next_round = []
            for task, res in out.items():
                if res.status is TaskStatus.OK:
                    continue
                # Unlike equilibrate_multi/findLiq_multi, EVERY failure is retried
                # here (up to max_domain_retries), not just TaskStatus.TIMEOUT --
                # matching the original saturation_pressure()'s own retry loop,
                # which can't even tell the two apart (findSatPressure_MELTS
                # returns the identical zeroed-out dict either way) and retries
                # anything still unconverged regardless. That's workable here,
                # unlike for a single equilibrate call, because the starting GUESS
                # itself (not just elapsed time) is genuinely perturbed below --
                # a calculation error from a bad bracket is exactly the kind of
                # failure a different P/T starting point can plausibly fix.
                if round_num < max_domain_retries:
                    # Perturb the starting guess for the retry, same spirit as the
                    # original saturation_pressure()'s own retry loop: alternate
                    # P_bar_init by x0.7/x2, and (only when this row is in mode 2 --
                    # T_C_init was actually given) nudge T_C_init by +25 then -15.
                    # Applied to whatever this task currently carries (not the
                    # original, un-retried value), so a second retry round would
                    # compound off the first's already-perturbed guess, matching
                    # the original's own compounding behaviour.
                    # T_fixed_C is NEVER perturbed -- mode 1 requires it stay fixed.
                    new_P_bar_init = task.P_bar_init * (0.7 if round_num % 2 == 0 else 2.0)
                    if task.T_C_init is not None:
                        new_T_C_init = task.T_C_init + (25.0 if round_num % 2 == 0 else -15.0)
                    else:
                        new_T_C_init = None
                    next_round.append(task._replace(P_bar_init = new_P_bar_init, T_C_init = new_T_C_init))
                elif task.row_id not in ticked:
                    ticked.add(task.row_id)
                    pbar.update(1)

            pending_tasks = next_round
            current_timeout = current_timeout * retry_timeout_multiplier
    finally:
        pbar.close()

    n_ok = len(collected)
    n_failed = n - n_ok
    print(f"saturation_pressure_v2: {n_ok}/{n} converged, {n_failed} did not "
          f"(any failure is retried up to {max_domain_retries} time(s), at {retry_timeout_multiplier}x "
          f"timeout and a perturbed starting guess each round).")

    if n == 1:
        return collected.get(0, {k: np.nan for k in
                                  ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq',
                                   'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq',
                                   'Fe3Fet_Liq', 'P_bar', 'T_Liq_C', 'X_H2O', 'X_CO2']})

    Results = pd.DataFrame([collected.get(i, {}) for i in range(n)])
    Results = Results.reindex(range(n))

    if copy_columns is not None:
        if type(copy_columns) == str:
            Results.insert(0, copy_columns, bulk[copy_columns])
        elif type(copy_columns) == list:
            for j, col in enumerate(copy_columns):
                Results.insert(j, col, bulk[col])

    return Results


# Plain alias: `saturation_pressure` IS this implementation as of 2026-10-21
# (see its own docstring), so `saturation_pressure_v2` is now redundant but
# kept so any existing code/notebook/scratch script calling it by that name
# keeps working.
saturation_pressure_v2 = saturation_pressure


def saturationP_multi(q, index, *, Model = None, comp = None, T_C_init = None, T_fixed_C = None,
               P_bar_init = None, fO2_buffer = None, fO2_offset = None, quiet_engine = True):
    f = open(os.devnull, 'w')
    sys.stdout = f

    # with tempfile.TemporaryDirectory() as tempdir:
    #     original_dir = os.getcwd()
    #     os.chdir(tempdir)

    #     try:
    results = {}
    idx = []
    if quiet_engine:
        _permanently_quiet_engine_output()   # one-way; see that function's own docstring
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

        melts.engine.setSystemProperties("Output",  "none")

        melts.engine.pressure = P_bar_init[0]
        if T_fixed_C[0] is not None:
            melts.engine.temperature = T_fixed_C[0] + 500
        else:
            melts.engine.temperature = T_C_init[0] + 500

        if type(comp) == dict:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'],
                    comp['Fe3Fet_Liq']*((159.69/2)/71.844)*comp['FeOt_Liq'],
                    comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'],
                    comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'],
                    comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'],
                    comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
        else:
            bulk = [comp.loc[0,'SiO2_Liq'], comp.loc[0,'TiO2_Liq'], comp.loc[0,'Al2O3_Liq'],
                    comp.loc[0,'Fe3Fet_Liq']*((159.69/2)/71.844)*comp.loc[0,'FeOt_Liq'],
                    comp.loc[0,'Cr2O3_Liq'], (1- comp.loc[0,'Fe3Fet_Liq'])*comp.loc[0,'FeOt_Liq'],
                    comp.loc[0,'MnO_Liq'], comp.loc[0,'MgO_Liq'], 0.0, 0.0, comp.loc[0,'CaO_Liq'],
                    comp.loc[0,'Na2O_Liq'], comp.loc[0,'K2O_Liq'], comp.loc[0,'P2O5_Liq'],
                    comp.loc[0,'H2O_Liq'], comp.loc[0,'CO2_Liq'], 0.0, 0.0, 0.0]

        melts.engine.setBulkComposition(bulk)
        PL = melts.engine.calcSaturationState()
        if T_fixed_C[0] is not None:
            for p in PL:
                if p != "fluid":
                    if p != "water":
                        melts.engine.setSystemProperties("Suppress", p)
        else:
            Suppress = ['rutile', 'tridymite']
            for p in Suppress:
                melts.engine.setSystemProperties("Suppress", p)

    for i in index:
        melts = melts.addNodeAfter()
        try:
            if type(comp) == dict:
                Results, tr = findSatPressure_MELTS(Model = Model, T_C_init=T_C_init[i], T_fixed_C=T_fixed_C[i],
                                            P_bar_init=P_bar_init[i], comp = comp, fO2_buffer=fO2_buffer,
                                            fO2_offset=fO2_offset[i], trial = True, melts = melts, suppressed = True)
            else:
                Results, tr = findSatPressure_MELTS(Model = Model, T_C_init=T_C_init[i], T_fixed_C=T_fixed_C[i],
                                            P_bar_init=P_bar_init[i], comp = comp.loc[i].to_dict(), fO2_buffer=fO2_buffer,
                                            fO2_offset=fO2_offset[i], trial = True, melts = melts, suppressed = True)

            results[f"Run {i}"] = Results
            idx.append(i)
            if tr is False:
                break
        except:
            idx.append(i)
            break

    q.put([idx, results])

        # finally:# Always return to the original directory before the process exits
        #     os.chdir(original_dir)
        #     q.put([idx, results])

def satP(q, index, *, Model = None, comp = None, T_C_init = None, T_fixed_C = None,
         P_bar_init = None, fO2_buffer = None, fO2_offset = None, quiet_engine = True):
    """Find the saturation pressure for a given composition and temperature.

    This function calculates the volatile saturation pressure for a given composition using the MELTS models. The results are returned in the form of a tuple and added to the specified queue.

    Args:
        q (Queue): The queue to which the results should be added.
        index (int): The index of the calculation.
        Model (str, optional): The model to use for the calculation. Must include "MELTS".
        comp (dict): A dictionary of oxide names and their concentration (in wt%) in the composition.
        T_C_init (float, optional): The initial temperature in degrees Celsius.
        P_bar_init (float, optional): The initial pressure in bar.
        fO2_buffer (float, optional): The fO2 buffer.
        fO2_offset (float, optional): The fO2 offset.
        quiet_engine (bool, optional, default True): If True, redirects the MELTS engine's stdout/stderr
            to the null device for the duration of the calculation. Set to False to inspect the raw
            engine output for debugging purposes.

    Returns:
        None
    """

    if "MELTS" in Model:
        if quiet_engine:
            _permanently_quiet_engine_output()   # one-way; see that function's own docstring
        Results = findSatPressure_MELTS(Model = Model, T_C_init = T_C_init, T_fixed_C=T_fixed_C, P_bar_init = P_bar_init, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index])
        return

def satP_multi(q, index, *, Model = None, comp = None, fO2_buffer = None, quiet_engine = True):
    warnings.warn("satP_multi() has no callers left in petthermotools (it only ever served the deprecated "
                  "findSatPressure_multi()) and is a candidate for removal in a future release.",
                  DeprecationWarning, stacklevel=2)
    if "MELTS" in Model:
        if quiet_engine:
            _permanently_quiet_engine_output()   # one-way; see that function's own docstring
        Results = findSatPressure_MELTS_multi(Model = Model, comp = comp, fO2_buffer=fO2_buffer,
                                              fO2_offset=comp['fO2_offset'],
                                              P_bar = comp['P_bar'], T_fixed_C=comp['T_fixed_C'])
        q.put([Results, index])
        return
    