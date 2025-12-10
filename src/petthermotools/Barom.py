import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    import Thermobar as pt
except:
    print('Thermobar cannot be imported, please check your numpy version')

from petthermotools.GenFuncs import *
from petthermotools.Plotting import *
from petthermotools.Liq import *
from petthermotools.MELTS import *
from petthermotools.Path import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
from scipy import interpolate
from shapely.geometry import MultiPoint, Point, Polygon

import numpy as np
import pandas as pd
import multiprocessing
import time
from pathlib import Path

def path_4_saturation_multi(q, index, *, Model = None, P_bar = None, comp = None, T_maxdrop_C = None, dt_C = None, T_initial_C = None, fO2_buffer = None,
                      fO2_offset = 0.0, H2O_Sat = None, phases = None):
    """
    Worker function to run crystallization simulations at a subset of pressures in parallel
    using alphaMELTS for Python or MAGEMinCalc to determine mineral co-saturations.

    Parameters
    ----------
    q : multiprocessing.Queue
        Queue object for inter-process communication (used to return results).
    index : list of int or numpy.ndarray
        Chunk of indices from the full pressure array (`P_bar`) assigned to this worker process.
    Model : str
        Thermodynamic model to use: e.g., 'MELTSv1.0.2', 'pMELTS', or 'MAGEMinCalc'.
    P_bar : numpy.ndarray
        The full array of pressures (in bar) that is being iterated over by all workers.
    comp : dict
        Bulk composition (oxides in wt%) of the starting material.
    T_maxdrop_C : float
        The minimum temperature (°C) below the liquidus to continue the simulation path down to.
    dt_C : float
        Step size for temperature change (in °C).
    T_initial_C : float
        Starting temperature (°C) for liquidus calculations.
    fO2_buffer : str, optional
        Oxygen fugacity buffer ('FMQ' or 'NNO' for MELTS, 'qfm' or 'nno' for MAGEMinCalc).
    fO2_offset : float, default 0.0
        Offset in log units from the specified fO2 buffer.
    H2O_Sat : bool
        If True, the system is simulated under water-saturated conditions.
    phases : list of str
        List of phases (e.g., ['quartz1', 'alkali-feldspar1']) to monitor for saturation.

    Returns
    -------
    None
        The results are placed into the multiprocessing queue `q` as a two-element list:
        `[list of indices processed (int), dict of raw simulation results]`.
    """

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
    else:
        from juliacall import Main as jl, convert as jlconvert
        env_dir = Path.home() / ".petthermotools_julia_env"
        jl_env_path = env_dir.as_posix()

        jl.seval(f"""
            import Pkg
            Pkg.activate("{jl_env_path}")
            """)

        jl.seval("using MAGEMinCalc")

        comp_julia = jl.seval("Dict")(comp)   

    for i in index:
        # try:
        if "MELTS" in Model:
            Results, tr = path_MELTS(P_bar = P_bar[i],
                    Model=Model,
                    comp=comp,
                    T_maxdrop_C=T_maxdrop_C,
                    dt_C=dt_C,
                    T_initial_C=T_initial_C,
                    find_liquidus=True,
                    fO2_buffer=fO2_buffer,
                    fO2_offset=fO2_offset,
                    fluid_sat=H2O_Sat,
                    Suppress=['rutile', 'tridymite'],
                    phases=phases,
                    trail = True,
                    melts = melts
                )
        else:
            if fO2_offset is None:
                fO2_offset = 0.0

            Results_df = jl.MAGEMinCalc.path(
                                    comp=comp_julia, dt_C=dt_C,
                                    P_bar=P_bar[i], T_min_C = T_maxdrop_C,            
                                    Model=Model, fo2_buffer=fO2_buffer, 
                                    fo2_offset=fO2_offset, find_liquidus=True,
                                    frac_xtal = False
                                    )
            
            Results = dict(Results_df)
            
            tr = True

        idx.append(i)


        results[f"Run {i}"] = Results

        if tr is False:
            break
        # except:
        #     idx.append(i)

    q.put([idx, results])
    return

def mineral_cosaturation(Model="MELTSv1.0.2", cores=int(np.floor(multiprocessing.cpu_count())), bulk=None,
                         phases=['quartz1', 'alkali-feldspar1', 'plagioclase1'],
                         P_bar=np.linspace(250, 5000, 32), Fe3Fet_init=None, H2O_init=None,
                         CO2_init=None, H2O_Sat=False, T_initial_C=None, dt_C=2,
                         T_maxdrop_C=50, T_cut_C=20, find_range=True,
                         find_min=True, fO2_buffer=None, fO2_offset=0.0, timeout = 90, multi_processing = True):
    """
    Determines the pressure at which two or more minerals co-saturate (following the method of Gualda et al. 2014).

    The analysis involves running multiple isobaric crystallization paths over a range of
    pressures, calculating the saturation temperature for the specified phases at each pressure,
    and finding the pressure that minimizes the temperature difference between them.

    Parameters
    ----------
    Model : str
        Thermodynamic model to use.
    cores : int
        Number of parallel processes to use.
    bulk : dict
        Bulk composition of the starting material, reported as oxide concentrations in wt%.
    phases : list
        Mineral phases to track for saturation. Must be a 2 or 3 item list ['quartz1', 'plagioclase1']
    P_bar : array-like
        Array of pressures used for the crystallization calculations.
    Fe3Fet_init : float, optional
        Initial Fe³⁺/Fe_total ratio.
    H2O_init : float, optional
        Initial H₂O content in the system (wt%).
    CO2_init : float, optional
        Initial CO₂ content in the system (wt%).
    H2O_Sat : bool
        If True, run at H₂O saturation. System must have enough H2O to be saturated.
    T_initial_C : float, optional
        Starting temperature in °C for liquidus calculations.
    dt_C : float
        Temperature step size (°C).
    T_maxdrop_C : float
        Maximum temperature drop below the liquidus to search for satuartion (°C).
    T_cut_C : float
        Maximum acceptable difference in phase saturation temperatures for "co-saturation" (°C). Used when trying to fit a polynomial to the results to find the true minimum point.
    find_range : bool
        If True, analyze phase co-saturation ranges. Currently disabled.
    find_min : bool
        If True, determine pressure where phase saturation temperature difference is minimized.
    fO2_buffer : str, optional
        Oxygen fugacity buffer (e.g., "FMQ", "NNO").
    fO2_offset : float, default 0.0
        Offset from the fO2 buffer in log units.
    timeout : int
        Timeout (in seconds) for each process.
    multi_processing : bool
        If True, run simulations in parallel using multiprocessing.

    Returns
    -------
    out : dict
        Contains either 'Output' (DataFrame of phase saturation conditions) or both 'Output' and 'CurveMin' (results from the minimization calculations).
    Results : dict
        Raw simulation outputs for each pressure step.
    """
    ## make sure everything is a float
    T_initial_C = to_float(T_initial_C)
    T_maxdrop_C = to_float(T_maxdrop_C)
    T_cut_C = to_float(T_cut_C)

    P_bar      = to_float(P_bar)

    Fe3Fet_init= to_float(Fe3Fet_init)
    H2O_init   = to_float(H2O_init)
    CO2_init   = to_float(CO2_init)
    fO2_offset = to_float(fO2_offset)
    
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
            
    if H2O_Sat:
        comp['H2O_Liq'] = 20

    comp = comp_fix(Model=Model, comp=comp, Fe3Fet_Liq=Fe3Fet_init,
                    H2O_Liq=H2O_init, CO2_Liq=CO2_init)

    combined_results = {}
    if multi_processing:
        index_in = np.arange(len(P_bar))
        index_out = np.array([], dtype=int)

        while len(index_out) < len(index_in):
            Start = time.time()
            index = np.setdiff1d(index_in, index_out)
            groups = np.array_split(index, cores)
            non_empty_groups = [g for g in groups if g.size > 0]
            groups = non_empty_groups

            processes = []
            for i in range(len(groups)):
                q = Queue()
                p = Process(target = path_4_saturation_multi, args = (q, groups[i]),
                            kwargs = {'Model': Model, 'comp': comp,
                            'T_initial_C': T_initial_C, 'T_maxdrop_C': T_maxdrop_C,
                            'dt_C': dt_C, 'P_bar': P_bar, 'phases': phases,
                            'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
                
                processes.append((p, q, groups[i]))
                p.start()

            for p, q, group in processes:
                try:
                    if time.time() - Start > timeout + 10:
                        res = q.get(timeout = 10)
                    else:
                        res = q.get(timeout=timeout)
                except:
                    if "MELTS" not in Model:
                        print(f"Timeout warning reached. Calculation P_bar = {P_bar[group[0]]} will not be returned. Try increasing the timeout.")
                    res = []
                    idx_chunks = np.array([group[0]], dtype = int)
                    mask = np.ones(len(P_bar), dtype=bool)
                    mask[idx_chunks] = False
                    P_bar = P_bar[mask]

                p.join(timeout = 2)
                p.terminate()

                if len(res) > 0.0:
                    idx_chunks, results = res
                    combined_results.update(results)

                index_out = np.hstack((index_out, idx_chunks))
    else:
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
            from juliacall import Main as jl, convert as jlconvert
            env_dir = Path.home() / ".petthermotools_julia_env"
            jl_env_path = env_dir.as_posix()

            jl.seval(f"""
                import Pkg
                Pkg.activate("{jl_env_path}")
                """)

            jl.seval("using MAGEMinCalc")

            comp_julia = jl.seval("Dict")(comp) 

        if "MELTS" in Model:
            for i in range(len(P_bar)):
                results = {}
                Results = path_MELTS(P_bar = P_bar[i],
                        Model=Model,
                        comp=comp,
                        T_maxdrop_C=T_maxdrop_C,
                        dt_C=dt_C,
                        T_initial_C=T_initial_C,
                        find_liquidus=True,
                        fO2_buffer=fO2_buffer,
                        fO2_offset=fO2_offset,
                        fluid_sat=H2O_Sat,
                        Suppress=['rutile', 'tridymite'],
                        phases=phases,
                        melts = melts
                    )
                results[f"Run {i}"] = Results
                combined_results.update(results)
        else:
            if fO2_offset is None:
                fO2_offset = 0.0

            for i in range(len(P_bar)):
                results = {}
                Results_df = jl.MAGEMinCalc.path(
                                        comp=comp_julia, dt_C=dt_C,
                                        P_bar=P_bar[i], T_min_C = T_maxdrop_C,            
                                        Model=Model, fo2_buffer=fO2_buffer, 
                                        fo2_offset=fO2_offset, find_liquidus=True,
                                        frac_xtal = False
                                        )
                
                Results = dict(Results_df)
                results[f"Run {i}"] = Results
                combined_results.update(results)
            
    results = combined_results

    Results = stich(Res=results, multi=True, Model=Model)

    ## determine the offset between the phases
    if len(phases)  == 3:
        arr = np.zeros((len(Results.keys()), 4))
        arr2 = np.zeros((len(Results.keys()), 4))
        columns = ['P_bar'] + phases + [phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2], '3 Phase Saturation']
        for i, r in enumerate(Results.keys()):
            for idx, p in enumerate(phases):
                if p in Results[r]['mass_g'].keys():
                    arr[i, idx+1] = Results[r]['Conditions'].loc[Results[r]['mass_g'][p] > 0.0, 'T_C'].values[0]
                else:
                    arr[i, idx + 1] = np.nan

            arr[i,0] = Results[r]['Conditions']['P_bar'].loc[0]

        arr2[:, 0] = np.abs(arr[:,1] - arr[:,2])
        arr2[:, 1] = np.abs(arr[:,1] - arr[:,3])
        arr2[:, 2] = np.abs(arr[:,2] - arr[:,3])
        arr2[:, 3] = np.max(arr2[:,0:2], axis = 1)

        r_arr = np.hstack((arr,arr2))

        out = pd.DataFrame(data = r_arr, columns = columns)
        out = out.sort_values('P_bar').reset_index(drop=True)

    else:
        arr = np.zeros((len(Results.keys()),4))
        columns = ['P_bar'] + phases + [phases[0]+' - '+phases[1]]
        for i, r in enumerate(Results.keys()):
            for idx, p in enumerate(phases):
                if p in Results[r]['mass_g'].keys():
                    arr[i, idx+1] = Results[r]['Conditions'].loc[Results[r]['mass_g'][p] > 0.0, 'T_C'].values[0]
                else:
                    arr[i, idx + 1] = np.nan

            arr[i,0] = Results[r]['Conditions']['P_bar'].loc[0]

        arr[:,3] = np.abs(arr[:,1] - arr[:,2])

        out = pd.DataFrame(data = arr, columns = columns)
        out = out.sort_values('P_bar').reset_index(drop=True)

    if find_min:
        res = findmin(out = out, P_bar = P_bar, T_cut_C = T_cut_C)
        out = {'CurveMin': res, 'Output': out}
    else:
        out = {'Output': out}

    return out, Results

def findmin(out = None, P_bar = None, T_cut_C = None):
    """
    Finds the minimum temperature offset between mineral saturation temperatures as a function of pressure.

    Parameters
    ----------
    out : pd.DataFrame
        DataFrame containing phase saturation temperatures.
    P_bar : array-like
        Pressure values (bar) corresponding to the rows in `out`.
    T_cut_C : float
        Maximum acceptable temperature difference (°C) for co-saturation to be considered valid.

    Returns
    -------
    CurveMin : dict
        Dictionary containing:
            - 'P_min': Pressure of minimum saturation temperature difference.
            - 'Res_min': Minimum temperature difference.
            - 'y_new': Interpolated curve of temperature differences.
            - 'P_new': Pressure array for interpolated curve.
            - 'test': 'Pass' or 'Fail' depending on whether the result meets `T_cut_C`.
    """
    Res = out.copy()

    if '3 Phase Saturation' in list(Res.keys()):
        T_cut_old = T_cut_C
        T_cut_C = T_cut_C + np.nanmin(Res[4:])
        Minimum = list(Res.keys())[4:]
        CurveMin = {}
        for m in Minimum:
            if len(Res[m][~np.isnan(Res[m].values)]) > 2:
                y = Res[m][(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)].values
                x = P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)]

                try:
                    y_new = interpolate.UnivariateSpline(x, y, k = 3)

                    P_new = np.linspace(P_bar[P_bar == np.nanmin(P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)])], 
                                        P_bar[P_bar == np.nanmax(P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)])], 200)

                    NewMin = np.nanmin(y_new(P_new))
                    P_min = P_new[np.where(y_new(P_new) == NewMin)][0]
                    if NewMin < T_cut_old:
                        Test = 'Pass'
                    else:
                        Test = 'Fail'

                    CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new(P_new), 'P_new': P_new, 'test': Test}
                except:
                    try:
                        y_new = interpolate.UnivariateSpline(x, y, k = 2)

                        P_new = np.linspace(P_bar[P_bar == np.nanmin(P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)])], 
                                            P_bar[P_bar == np.nanmax(P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C)])], 200)

                        NewMin = np.nanmin(y_new(P_new))
                        P_min = P_new[np.where(y_new(P_new) == NewMin)][0]
                        if NewMin < T_cut_old:
                            Test = 'Pass'
                        else:
                            Test = 'Fail'

                        CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new(P_new), 'P_new': P_new, 'test': Test}
                    except:
                        CurveMin[m] = {'P_min': np.nan, 'Res_min': np.nan, 'y_new': np.nan, 'P_new': np.nan, 'test': 'Fail'}
            else:
                y_new = np.nan
                P_new = np.nan
                NewMin = np.nan
                P_min = np.nan
                Test = 'Fail'
                CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new, 'P_new': P_new, 'test': Test}

    else:
        CurveMin = {}
        m = Res.keys()[3]
        if len(Res[m][~np.isnan(Res[m])]) > 2:
            y = Res[m][(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C*2)].values
            x = P_bar[(~np.isnan(Res[m].values)) & (Res[m].values < T_cut_C*2)]

            try:
                y_new = interpolate.UnivariateSpline(x, y, k = 3)
            except:
                y_new = interpolate.UnivariateSpline(x, y, k = 2)

            P_new = np.linspace(P_bar[P_bar == np.nanmin(P_bar[(~np.isnan(Res[m])) & (Res[m] < T_cut_C*2)])], 
                                    P_bar[P_bar == np.nanmax(P_bar[(~np.isnan(Res[m])) & (Res[m] < T_cut_C*2)])], 200)
            
            NewMin = np.nanmin(y_new(P_new))
            P_min = P_new[np.where(y_new(P_new) == NewMin)][0]
            if NewMin < T_cut_C:
                Test = 'Pass'
            else:
                Test = 'Fail'
        else:
            y_new = np.nan
            P_new = np.nan
            NewMin = np.nan
            P_min = np.nan
            Test = 'Fail'

        CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new, 'P_new': P_new, 'test': Test}

    return CurveMin

def find_mineral_cosaturation(cores = None, Model = None, bulk = None, phases = None, P_bar = None, Fe3Fet_Liq = None, 
                              H2O_Liq = None, H2O_Sat = False, T_initial_C = None, dt_C = None, T_maxdrop_C = None, 
                              T_cut_C = None, find_range = None, find_min = None, fO2_buffer = None, fO2_offset = None):
    print("This function has been removed following update to v0.2.40. Please switch to using the mineral_cosaturation() function")

    return "This function has been removed following update to v0.2.40. Please switch to using the mineral_cosaturation() function"