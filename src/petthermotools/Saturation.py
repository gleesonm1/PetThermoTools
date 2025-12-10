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

def findSatPressure_multi(cores = multiprocessing.cpu_count(), Model = "MELTSv1.2.0", bulk = None, T_fixed_C = None, 
                          P_bar_init = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, 
                          fO2_buffer = None, fO2_offset = None, copy_columns = None):
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

    Returns:
    ----------
    Results : pd.DataFrame
        DataFrame containing the calculated saturation pressure ($P_{\text{sat}}$), the 
        corresponding melt composition, and the initial system parameters.
    
    Raises:
    -------
    Warning
        - If $\text{H}_2\text{O}$ is zero in a MELTS model (may affect oxide/apatite saturation).
        - If $\text{Fe}^{3+}/\Sigma\text{Fe}$ is zero without an $\text{f}\text{O}_2$ buffer in a MELTS model (may cause calculation failure).
    '''
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
                                'fO2_buffer': fO2_buffer})

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
                    fO2_buffer = None, fO2_offset = None, copy_columns = None):
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

    Returns:
    Res (dict or pandas.core.frame.DataFrame): A dict or dataframe containing the chemical composition of the melt phase and the saturation pressure and temperature of the melt.

    Examples:

    Calculate the saturation pressure and temperature of a basaltic liquid, with initial conditions set at 1200 C and 10000 bar using all available CPU cores and the MELTSv1.2.0 model.
    >>> findSatPressure(bulk = {'SiO2_Liq': 50, 'TiO2_Liq': 1, 'Al2O3_Liq': 15, 'FeOt_Liq': 10, 'MnO_Liq': 0.5, 'MgO_Liq': 8, 'CaO_Liq': 10, 'Na2O_Liq': 2, 'K2O_Liq': 1, 'P2O5_Liq': 0.5, 'H2O_Liq': 0.5, 'Fe3Fet_Liq': 0.15}, T_C_init = 1200, P_bar_init = 10000)

    """

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

    if cores is None:
        cores = multiprocessing.cpu_count()

    qs = []
    q = Queue()

    # specify the number of calculations to be performed in each sequence
    if type(comp) == dict:
        p = Process(target = satP, args = (q, 0),
                    kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                              'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

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
                                    'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset[i]})

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

        Res = pd.DataFrame(data = np.zeros((int(L),15)), columns = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq', 'P_bar', 'T_Liq'])
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


def saturation_pressure(Model = "MELTSv1.2.0", cores = multiprocessing.cpu_count(),
                        bulk = None, T_C_init = None, T_fixed_C = None, P_bar_init = 5000,
                        Fe3Fet_Liq = None, Fe3Fet_init = None, H2O_Liq = None, H2O_init = None,
                        CO2_Liq = None, CO2_init = None, fO2_buffer = None, fO2_offset = None, 
                        copy_columns = None, multi_processing = True, timeout = 180):
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
    cores : int, default multiprocessing.cpu_count()
        Number of parallel processes to use.
    bulk : dict, pandas.Series, or pandas.DataFrame
        Bulk melt composition(s). If a Series or single dict is passed, a single calculation is run.
        If a DataFrame is passed, a calculation is performed for each row.
    T_C_init : float or np.ndarray, optional
        Initial temperature(s) in degrees Celsius. System will co-solve the liquidus temperature and pressure of volatile saturation. Used if not specifying `T_fixed_C`.
    T_fixed_C : float or np.ndarray, optional
        Fixed temperature(s) for the saturation pressure calculation. Solid phases are disabled in all calculations.
    P_bar_init : float or np.ndarray, default 5000
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
       
    if "MELTS" in Model:
        try:
            from meltsdynamic import MELTSdynamic
        except:
            Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

    if bulk is not None:
        comp = bulk.copy()

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

    qs = []
    q = Queue()
    if One == 1:
        if multi_processing:
            p = Process(target = satP, args = (q,0),
                        kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                                  'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init,
                                  'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
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

            if len(groups[0])< 15:
                timeout = len(groups[0])*timeout_main
            else:
                timeout = 15*timeout_main
            
            processes = []
            Start = time.time()
            for i in range(len(groups)):
                q = Queue()
                p = Process(target = saturationP_multi, args = (q, groups[i]),
                            kwargs = {'Model': Model, 'comp': comp, 'T_C_init': T_C_init,
                                      'T_fixed_C': T_fixed_C, 'P_bar_init': P_bar_init, 
                                      'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                
                processes.append([p, q, groups[i]])
                p.start()

            for p, q, group in processes:
                try:
                    if time.time() - Start < timeout + 5:
                        res = q.get(timeout = timeout)
                    else:
                        res = q.get(timeout=10)
                except:
                    res = []
                    idx_chunks = np.array([group[0]], dtype = int)

                p.join(timeout=2)
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
            df = df.reindex(full_index, fill_value=0.0)

            # Optional: sort index if needed
            df = df.sort_index()

            return df
        Results = results_to_dataframe(results, index_in)
        
        if copy_columns is not None:
            if type(copy_columns) == str:
                Results.insert(0, copy_columns, comp[copy_columns])
                # Af_Combined.insert(0, copy_columns, comp[copy_columns])
            elif type(copy_columns) == list:
                j = 0
                for i in copy_columns:
                    Results.insert(j, i, comp[i])
                    # Af_Combined.insert(j, i, comp[i])
                    j = j + 1

        return Results

def saturationP_multi(q, index, *, Model = None, comp = None, T_C_init = None, T_fixed_C = None,
               P_bar_init = None, fO2_buffer = None, fO2_offset = None):
    results = {}
    idx = []
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

        
        melts.engine.pressure = P_bar_init[0]
        if T_fixed_C[0] is not None:
            melts.engine.temperature = T_fixed_C[0] + 500
        else:
            melts.engine.temperature = T_C_init[0] + 500

        if type(comp) == dict:
            bulk = [comp['SiO2_Liq'], comp['TiO2_Liq'], comp['Al2O3_Liq'], 
                    comp['Fe3Fet_Liq']*((159.59/2)/71.844)*comp['FeOt_Liq'], 
                    comp['Cr2O3_Liq'], (1- comp['Fe3Fet_Liq'])*comp['FeOt_Liq'], 
                    comp['MnO_Liq'], comp['MgO_Liq'], 0.0, 0.0, comp['CaO_Liq'], 
                    comp['Na2O_Liq'], comp['K2O_Liq'], comp['P2O5_Liq'], 
                    comp['H2O_Liq'], comp['CO2_Liq'], 0.0, 0.0, 0.0]
        else:
            bulk = [comp.loc[0,'SiO2_Liq'], comp.loc[0,'TiO2_Liq'], comp.loc[0,'Al2O3_Liq'], 
                    comp.loc[0,'Fe3Fet_Liq']*((159.59/2)/71.844)*comp.loc[0,'FeOt_Liq'], 
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

def satP(q, index, *, Model = None, comp = None, T_C_init = None, T_fixed_C = None, 
         P_bar_init = None, fO2_buffer = None, fO2_offset = None):
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

    Returns:
        None
    """

    if "MELTS" in Model:
        Results = findSatPressure_MELTS(Model = Model, T_C_init = T_C_init, T_fixed_C=T_fixed_C, P_bar_init = P_bar_init, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        q.put([Results, index])
        return
    
def satP_multi(q, index, *, Model = None, comp = None, fO2_buffer = None):
    if "MELTS" in Model:
        Results = findSatPressure_MELTS_multi(Model = Model, comp = comp, fO2_buffer=fO2_buffer, 
                                              fO2_offset=comp['fO2_offset'],
                                              P_bar = comp['P_bar'], T_fixed_C=comp['T_fixed_C'])
        q.put([Results, index])
        return
    