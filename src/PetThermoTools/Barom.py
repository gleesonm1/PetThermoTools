import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    import Thermobar as pt
except:
    print('Thermobar cannot be imported, please check your numpy version')

from PetThermoTools.GenFuncs import *
from PetThermoTools.Plotting import *
from PetThermoTools.Liq import *
from PetThermoTools.MELTS import *
from PetThermoTools.Path import *
# try:
#     from PetThermoTools.Holland import *
# except:
#     pass
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.notebook import tqdm, trange
from scipy import interpolate
from shapely.geometry import MultiPoint, Point, Polygon

from functools import partial
import concurrent.futures

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import time
        
def path_4_saturation(q, index, *, Model = None, P_bar = None, comp = None, T_maxdrop_C = None, dt_C = None, T_initial_C = None, fO2_buffer = None,
                      fO2_offset = 0.0, H2O_Sat = None, phases = None):
    if "MELTS" in Model:
        Results = path_MELTS(P_bar = P_bar,
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
                phases=phases
            )
        q.put([index, Results])
    return

def path_4_saturation_multi(q, index, *, Model = None, P_bar = None, comp = None, T_maxdrop_C = None, dt_C = None, T_initial_C = None, fO2_buffer = None,
                      fO2_offset = 0.0, H2O_Sat = None, phases = None):
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


        results[f"index = {i}"] = Results

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
                         find_min=True, fO2_buffer=None, fO2_offset=0.0, timeout = 30):

    comp = bulk.copy()
    if H2O_Sat:
        comp['H2O_Liq'] = 20

    comp = comp_fix(Model=Model, comp=comp, Fe3Fet_Liq=Fe3Fet_init,
                    H2O_Liq=H2O_init, CO2_Liq=CO2_init)

    index_in = np.arange(len(P_bar))
    # index_in = index
    combined_results = {}
    index_out = np.array([], dtype=int)

    while len(index_out) < len(index_in):
        index = np.setdiff1d(index_in, index_out)
        groups = np.array_split(index, cores)

        processes = []
        for i in range(len(groups)):
            q = Queue()
            p = Process(target = path_4_saturation_multi, args = (q, groups[i]),
                        kwargs = {'Model': Model, 'comp': comp,
                        'T_initial_C': T_initial_C, 'T_maxdrop_C': T_maxdrop_C,
                        'dt_C': dt_C, 'P_bar': P_bar, 'phases': phases,
                        'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
            
            processes.append((p, q))
            p.start()

        for p, q in processes:
            res = q.get(timeout=timeout)
            p.join(timeout = 2)
            p.terminate()

            idx_chunks, results = res

            index_out = np.hstack((index_out, idx_chunks))
            combined_results.update(results)

    results = combined_results
    # results = {}
    # if len(P_bar) > 1:
    #     A = len(P_bar)//cores
    #     B = len(P_bar) % cores

    # if A > 0:
    #     Group = np.zeros(A) + cores
    #     if B > 0:
    #         Group = np.append(Group, B)
    # else:
    #     Group = np.array([B])

    # # initialise queue
    # qs = []
    # q = Queue()


    # # run calculations
    # for k in range(len(Group)):
    #     ps = []

    #     for kk in range(int(cores*k), int(cores*k + Group[k])):
    #         p = Process(target = path_4_saturation, args = (q, kk),
    #                     kwargs = {'Model': Model, 'comp': comp,
    #                     'T_initial_C': T_initial_C, 'T_maxdrop_C': T_maxdrop_C,
    #                     'dt_C': dt_C, 'P_bar': P_bar[kk], 'phases': phases,
    #                     'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})
            
    #         ps.append(p)
    #         p.start()

    #     start = time.time()
    #     for p in ps:
    #         if time.time() - start < timeout - 10:
    #             try:
    #                 ret = q.get(timeout = timeout - (time.time()-start) + 10)
    #             except:
    #                 ret = []
    #         else:
    #             try:
    #                 ret = q.get(timeout = 10)
    #             except:
    #                 ret = []

    #         qs.append(ret)

    #     TIMEOUT = 5
    #     start = time.time()
    #     for p in ps:
    #         if p.is_alive():
    #             while time.time() - start <= TIMEOUT:
    #                 if not p.is_alive():
    #                     p.join()
    #                     p.terminate()
    #                     break
    #                 time.sleep(.1)
    #             else:
    #                 p.terminate()
    #                 p.join(5)
    #         else:
    #             p.join()
    #             p.terminate()

    # for kk in range(len(qs)):
    #     if len(qs[kk]) > 0:
    #         index, res = qs[kk]
    #         results[f"P = {P_bar[index]:.2f} bars"] = res.copy()

    Results = stich(Res=results, multi=True, Model=Model)

    ## determine the offset between the phases
    if len(phases)  == 3:
        arr = np.zeros((len(Results.keys()), 4))
        arr2 = np.zeros((len(Results.keys()), 4))
        columns = ['P_bar'] + phases + [phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2], '3 Phase Saturation']
        for i, r in enumerate(Results.keys()):
            for idx, p in enumerate(phases):
                if p in Results[r]['Mass'].keys():
                    arr[i, idx+1] = Results[r]['Conditions'].loc[Results[r]['Mass'][p] > 0.0, 'T_C'].values[0]
                else:
                    arr[i, idx + 1] = np.nan

            arr[i,0] = Results[r]['Conditions']['P_bar'].loc[0]

        arr2[:, 0] = np.abs(arr[:,1] - arr[:,2])
        arr2[:, 1] = np.abs(arr[:,1] - arr[:,3])
        arr2[:, 2] = np.abs(arr[:,2] - arr[:,3])
        arr2[:,3] = np.max(arr2[:,0:2], axis = 1)

        r_arr = np.hstack((arr,arr2))

        out = pd.DataFrame(data = r_arr, columns = columns)
        out = out.sort_values('P_bar').reset_index(drop=True)

    else:
        arr = np.zeros((len(Results.keys()),4))
        columns = ['P_bar'] + phases + [phases[0]+' - '+phases[1]]
        for i, r in enumerate(Results.keys()):
            for idx, p in enumerate(phases):
                if p in Results[r]['Mass'].keys():
                    arr[i, idx+1] = Results[r]['Conditions'].loc[Results[r]['Mass'][p] > 0.0, 'T_C'].values[0]
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
    '''
    Carry out multiple calculations in parallel. Allows isobaric, polybaric and isochoric crystallisation to be performed as well as isothermal, isenthalpic or isentropic decompression. All temperature inputs/outputs are reported in degrees celcius and pressure is reported in bars.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default will be determined using Multiprocessing.cpu_count().

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    bulk: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    phases: list
        length 2 or 3, contains the phases of the co-saturated magma. Default = ['quartz1', 'plagioclase1', 'k-feldspar1'].

    P_bar: np.ndarray
        Calculation pressure. Length determines the number of calculations to be performed.

    Fe3Fet_Liq: float or np.ndarray
        Initial Fe 3+/total ratio.

    H2O_Liq: float
        H2O content of the initial melt phase.

    T_initial_C: float
        Starting temperature for the liquidus calculations. Default = 1200

    T_maxdrop_C: float
        Max temperature drop of the calculations. Default = 25

    dt_C: float
        Temperature change at each model step. Default = 1

    T_cut_C: float
        Temperature offset used to indicate whether the model has succeeded or failed in finding a match. Default = 10.

    find_range: True/False
        If True a new DataFrame will be included in Results, indicating all cases where the minimum offset is less than T_cut_C.

    find_min: True/False
        If True, a spline fit will be applied to the data to find the minimum point.

    fO2_buffer: string
        If the oxygen fugacity of the system is to be buffered during crystallisation/decompression, then an offset to a known buffer must be specified. Here the user can define the known buffer as either "FMQ" or "NNO".

    fO2_offset: float
        Offset from the buffer spcified in fO2_buffer (log units).

    Returns:
    ----------
    Results: Dict
        Dictionary containing information regarding the saturation temperature of each phase and the residuals between the different phases
    '''
    try:
        from meltsdynamic import MELTSdynamic
    except:
        Warning('alphaMELTS for Python files are not on the python path. \n Please add these files to the path running \n import sys \n sys.path.append(r"insert_your_path_to_melts_here") \n You are looking for the location of the meltsdynamic.py file')

    T_step_C = dt_C
    dt_C = T_maxdrop_C 

    comp = bulk.copy()
    if H2O_Sat is True:
        comp['H2O_Liq'] = 20

    # set default values if required
    if Model is None:
        Model == "MELTSv1.0.2"

    if cores is None:
        cores = multiprocessing.cpu_count()

    # if comp is entered as a pandas series, it must first be converted to a dict
    if type(comp) == pd.core.series.Series:
        comp = comp.to_dict()

    # ensure the bulk composition has the correct headers etc.
    comp = comp_fix(Model = Model, comp = comp)

    # create base array to be filled
    if P_bar is None:
        P_bar = np.array([-1])
    elif type(P_bar) != np.ndarray:
        P_bar = np.array([P_bar])

    if Fe3Fet_Liq is None:
        Fe3Fet_Liq = np.array([-1])
    elif type(Fe3Fet_Liq) != np.ndarray:
        Fe3Fet_Liq = np.array([Fe3Fet_Liq])

    if H2O_Liq is None:
        H2O_Liq = np.array([-1])
    elif type(H2O_Liq) != np.ndarray:
        H2O_Liq = np.array([H2O_Liq])

    base_array = np.zeros((len(H2O_Liq),len(Fe3Fet_Liq),len(P_bar)))

    if P_bar[0] == -1:
        P_bar = 1000
    if Fe3Fet_Liq[0] == -1:
        Fe3Fet_Liq = None
    if H2O_Liq[0] == -1:
        H2O_Liq = None

    # set default values for remaining parameters
    if T_initial_C is None:
        T_initial_C = 1200

    if dt_C is None:
        dt_C = 25

    if T_step_C is None:
        T_step_C = 1

    if T_cut_C is None:
        T_cut_C = 10

    if phases is None:
        phases = ['quartz1', 'alkali-feldspar1', 'plagioclase1']

    # create main output dictionary
    Results = {}
    if len(phases) == 3:
        List = [phases[0], phases[1], phases[2], 'T_Liq', 'H2O_melt', '3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]
        for l in List:
            Results[l] = base_array.copy()
    else:
        List = [phases[0], phases[1], 'T_Liq', 'H2O_melt', phases[0] + ' - ' + phases[1]]
        for l in List:
            Results[l] = base_array.copy()

    # run calculations if only one initial composition provided
    if type(comp) == dict:
        for i in range(np.shape(base_array)[0]):
            # if H2O specified set H2O on each iteration
            if H2O_Liq is not None:
                comp['H2O_Liq'] = H2O_Liq[i]

            for j in range(np.shape(base_array)[1]):
                #if Fe3Fet specified set Fe3Fet on each iteration
                if Fe3Fet_Liq is not None:
                    comp['Fe3Fet_Liq'] = Fe3Fet_Liq[j]

                # determine how many processes to run in parallel
                if len(P_bar) > 1:
                    A = len(P_bar)//cores
                    B = len(P_bar) % cores

                if A > 0:
                    Group = np.zeros(A) + cores
                    if B > 0:
                        Group = np.append(Group, B)
                else:
                    Group = np.array([B])

                # initialise queue
                qs = []
                q = Queue()


                # run calculations
                for k in range(len(Group)):
                    ps = []

                    for kk in range(int(cores*k), int(cores*k + Group[k])):
                        p = Process(target = satTemperature, args = (q, kk),
                                    kwargs = {'Model': Model, 'comp': comp,
                                    'T_initial_C': T_initial_C, 'T_step_C': T_step_C,
                                    'dt_C': dt_C, 'P_bar': P_bar[kk], 'phases': phases,
                                    'H2O_Liq': H2O_Liq, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

                        ps.append(p)
                        p.start()

                    TIMEOUT = 300

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


                    # for p in ps:
                    #     try:
                    #         ret = q.get(timeout = 180)
                    #     except:
                    #         ret = []
                    #
                    #     qs.append(ret)
                    #
                    # TIMEOUT = 20
                    # start = time.time()
                    # for p in ps:
                    #     if p.is_alive():
                    #         time.sleep(.1)
                    #         while time.time() - start <= TIMEOUT:
                    #             if not p.is_alive():
                    #                 p.join()
                    #                 p.terminate()
                    #                 break
                    #             time.sleep(.1)
                    #         else:
                    #             p.terminate()
                    #             p.join(5)
                    #     else:
                    #         p.join()
                    #         p.terminate()

                # # extract results
                for kk in range(len(qs)):
                    if len(qs[kk]) > 0:
                        Res, index = qs[kk]
                        for l in Results:
                            if l != 'sat_surface':
                                Results[l][i,j,index] = Res[l]

        # covert any empty values to nan
        for l in Results:
            if l != 'sat_surface':
                Results[l][np.where(Results[l] == 0.0)] = np.nan

        if find_min is not None:
            if H2O_Liq is not None:
                Results = findMinimum(Results = Results, P_bar = P_bar, T_cut_C = T_cut_C, H2O_Liq = H2O_Liq, Fe3Fet_Liq = Fe3Fet_Liq, phases = phases)
            else:
                Results = findMinimum(Results = Results, P_bar = P_bar, T_cut_C = T_cut_C, H2O_Liq = H2O_Liq, Fe3Fet_Liq = Fe3Fet_Liq, phases = phases)

        if find_range is not None:
            Results['range'] = np.zeros(np.shape(Results[l]))
            if len(phases) == 3:
                Results['range'][np.where(Results['3 Phase Saturation'] <= T_cut_C)] = True
                Results['range'][np.where(Results['3 Phase Saturation'] > T_cut_C)] = False
            else:
                Results['range'][np.where(Results[phases[0] + ' - ' + phases[1]] <= T_cut_C)] = True
                Results['range'][np.where(Results[phases[0] + ' - ' + phases[1]] > T_cut_C)] = False

    

    return Results

def findMinimum(Results = None, P_bar = None, T_cut_C = None, H2O_Liq = None, Fe3Fet_Liq = None, phases = None):
    '''
    Take the results of SatPress and search for the minimum point using a spline fit.
    '''
    if T_cut_C is None:
        T_cut_C = 10

    if H2O_Liq is None and Fe3Fet_Liq is None:
        if '3 Phase Saturation' in list(Results.keys()):
            Minimum = {'3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]}
            CurveMin = {}
            for m in Minimum:
                if len(Results[m][0,0,:][~np.isnan(Results[m][0,0,:])]) > 2:
                    y = Results[m][0,0,:][(np.where(~np.isnan(Results[m][0,0,:]))) and (np.where(Results[m][0,0,:] < T_cut_C))]
                    x = P_bar[(np.where(~np.isnan(Results[m][0,0,:]))) and (np.where(Results[m][0,0,:] < T_cut_C))]

                    try:
                        y_new = interpolate.UnivariateSpline(x, y, k = 5)
                    except:
                        y_new = interpolate.UnivariateSpline(x, y, k = 2)

                    P_new = np.linspace(P_bar[np.where(P_bar == np.nanmin(P_bar[(np.where(~np.isnan(Results[m][0,0,:]))) and (np.where(Results[m][0,0,:] < T_cut_C))]))], P_bar[np.where(P_bar == np.nanmax(P_bar[(np.where(~np.isnan(Results[m][0,0,:]))) and (np.where(Results[m][0,0,:] < T_cut_C))]))], 200)

                    NewMin = np.nanmin(y_new(P_new))
                    P_min = P_new[np.where(y_new(P_new) == NewMin)][0]
                    if NewMin < T_cut_C:
                        Test = 'Pass'
                    else:
                        Test = 'Fail'

                    CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new(P_new), 'P_new': P_new, 'test': Test}
                else:
                    y_new = np.nan
                    P_new = np.nan
                    NewMin = np.nan
                    P_min = np.nan
                    Test = 'Fail'
                    CurveMin[m] = {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new, 'P_new': P_new, 'test': Test}

            Results['CurveMin'] = CurveMin
        else:
            m = phases[0] + ' - ' + phases[1]
            if len(Results[m][0,0,:][~np.isnan(Results[m][0,0,:])]) > 2:
                y = Results[m][0,0,:][np.where(~np.isnan(Results[m][0,0,:]))]
                x = P_bar[np.where(~np.isnan(Results[m][0,0,:]))]

                try:
                    y_new = interpolate.UnivariateSpline(x, y, k=5)
                except:
                    y_new = interpolate.UnivariateSpline(x, y, k = 2)

                P_new = np.linspace(P_bar[np.where(P_bar == np.nanmin(P_bar[np.where(~np.isnan(Results[m][0,0,:]))]))], P_bar[np.where(P_bar == np.nanmax(P_bar[np.where(~np.isnan(Results[m][0,0,:]))]))], 200)

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

            Results['CurveMin'] = {phases[0] + ' - ' + phases[1]: {'P_min': P_min, 'Res_min': NewMin, 'y_new': y_new, 'P_new': P_new, 'test': Test}}

    if H2O_Liq is not None and Fe3Fet_Liq is None:
        if '3 Phase Saturation' in list(Results.keys()):
            X, Y = np.meshgrid(P_bar, H2O_Liq)
            Y = Results['H2O_melt'][:,0,:].copy()

            Minimum = {'3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]}
            CurveMin = {}
            for m in Minimum:
                if len(Results[m][:,0,:][~np.isnan(Results[m][:,0,:])]) > 4:
                    Res = Results[m][:,0,:].copy()
                    Res[np.where(Res > T_cut_C*2)] = np.nan
                    for i in range(len(H2O_Liq)):
                        Res[i, :][np.where(Results['H2O_melt'][i,0,:] < 0.99*np.nanmax(Results['H2O_melt'][i,0,:]))] = np.nan

                    A = Res.copy()
                    Res[np.where(Res > T_cut_C)] = np.nan
                    X, Y = np.meshgrid(P_bar, H2O_Liq)
                    Y = Results['H2O_melt'][:,0,:]

                    try:
                        z_new = interpolate.SmoothBivariateSpline(X[np.where(~np.isnan(A) & ~np.isnan(Y))], Y[np.where(~np.isnan(A) & ~np.isnan(Y))], A[np.where(~np.isnan(A) & ~np.isnan(Y))], kx = 3, ky = 3)
                    except:
                        z_new = interpolate.SmoothBivariateSpline(X[np.where(~np.isnan(A) & ~np.isnan(Y))], Y[np.where(~np.isnan(A) & ~np.isnan(Y))], A[np.where(~np.isnan(A) & ~np.isnan(Y))], kx = 2, ky = 2)

                    H2O_new = np.linspace(Y[np.where(Res == np.nanmin(Res))] - (H2O_Liq[1]-H2O_Liq[0]),
                                        Y[np.where(Res == np.nanmin(Res))] + (H2O_Liq[1]-H2O_Liq[0]), 20)
                    P_new = np.linspace(X[np.where(Res == np.nanmin(Res))] - (P_bar[1]-P_bar[0]),
                                        X[np.where(Res == np.nanmin(Res))] + (P_bar[1]-P_bar[0]), 20)

                    X_new, Y_new = np.meshgrid(P_new, H2O_new)
                    x = X[~np.isnan(A)].flatten()
                    y = Y[~np.isnan(A)].flatten()

                    MyPoly = MultiPoint(list(zip(x, y))).convex_hull

                    points = list(zip(X_new.flatten(), Y_new.flatten()))
                    Include = np.zeros(len(X_new.flatten()))
                    for i in range(len(points)):
                        p = Point(points[i])
                        Include[i] = p.within(MyPoly)

                    YayNay = Include.reshape(X_new.shape)
                    x_new = X_new[np.where(YayNay == True)].flatten()
                    y_new = Y_new[np.where(YayNay == True)].flatten()
                    Res_min = np.nanmin(z_new(x_new, y_new, grid = False))
                    P_min = x_new[np.where(z_new(x_new, y_new, grid = False) == Res_min)]
                    H2O_min = y_new[np.where(z_new(x_new, y_new, grid = False) == Res_min)]
                    if Res_min < T_cut_C:
                        Test = 'Pass'
                    else:
                        Test = 'Fail'

                    CurveMin[m] = {'Res_min': Res_min, 'P_min': P_min[0], 'H2O_min': H2O_min[0], 'z_new': z_new, 'test': Test}
                else:
                    CurveMin[m] = {'Res_min': np.nan, 'P_min': np.nan, 'H2O_min': np.nan, 'z_new': np.nan, 'test': 'Fail'}

            Results['CurveMin'] = CurveMin

    if H2O_Liq is not None and Fe3Fet_Liq is not None:
        if '3 Phase Saturation' in list(Results.keys()):
            X, Y = np.meshgrid(P_bar, H2O_Liq)
            Y = Results['H2O_melt'][:,0,:].copy()

            Minimum = {'3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]}
            CurveMin = {}
            for m in Minimum:
                Res_min_save = 50
                for w in range(len(Fe3Fet_Liq)):
                    if len(Results[m][:,w,:][~np.isnan(Results[m][:,w,:])]) > 4:
                        Res = Results[m][:,w,:].copy()
                        Res[np.where(Res > T_cut_C*2)] = np.nan
                        for i in range(len(H2O_Liq)):
                            Res[i, :][np.where(Results['H2O_melt'][i,w,:] < 0.99*np.nanmax(Results['H2O_melt'][i,w,:]))] = np.nan

                        A = Res.copy()
                        X, Y = np.meshgrid(P_bar, H2O_Liq)
                        Y = Results['H2O_melt'][:,w,:].copy()

                        try:
                            z_new = interpolate.SmoothBivariateSpline(X[np.where(~np.isnan(A) & ~np.isnan(Y))], Y[np.where(~np.isnan(A) & ~np.isnan(Y))], A[np.where(~np.isnan(A) & ~np.isnan(Y))], kx = 3, ky = 3)
                        except:
                            z_new = interpolate.SmoothBivariateSpline(X[np.where(~np.isnan(A) & ~np.isnan(Y))], Y[np.where(~np.isnan(A) & ~np.isnan(Y))], A[np.where(~np.isnan(A) & ~np.isnan(Y))], kx = 2, ky = 2)

                        H2O_new = np.linspace(Y[np.where(Res == np.nanmin(Res))] - (H2O_Liq[1]-H2O_Liq[0]),
                                            Y[np.where(Res == np.nanmin(Res))] + (H2O_Liq[1]-H2O_Liq[0]), 20)
                        P_new = np.linspace(X[np.where(Res == np.nanmin(Res))] - (P_bar[1]-P_bar[0]),
                                            X[np.where(Res == np.nanmin(Res))] + (P_bar[1]-P_bar[0]), 20)

                        X_new, Y_new = np.meshgrid(P_new, H2O_new)
                        x = X[~np.isnan(A)].flatten()
                        y = Y[~np.isnan(A)].flatten()

                        MyPoly = MultiPoint(list(zip(x, y))).convex_hull

                        points = list(zip(X_new.flatten(), Y_new.flatten()))
                        Include = np.zeros(len(X_new.flatten()))
                        for i in range(len(points)):
                            p = Point(points[i])
                            Include[i] = p.within(MyPoly)

                        YayNay = Include.reshape(X_new.shape)
                        x_new = X_new[np.where(YayNay == True)].flatten()
                        y_new = Y_new[np.where(YayNay == True)].flatten()
                        Res_min = np.nanmin(z_new(x_new, y_new, grid = False))
                        P_min = x_new[np.where(z_new(x_new, y_new, grid = False) == Res_min)]
                        H2O_min = y_new[np.where(z_new(x_new, y_new, grid = False) == Res_min)]
                        if Res_min < T_cut_C:
                            Test = 'Pass'
                        else:
                            Test = 'Fail'

                        if Res_min < Res_min_save:
                            Res_min_save = Res_min
                            CurveMin[m] = {'Res_min': Res_min, 'P_min': P_min[0], 'H2O_min': H2O_min[0], 'z_new': z_new, 'Fe3Fet_Liq': Fe3Fet_Liq[w], 'test': Test}

            Results['CurveMin'] = CurveMin


    return Results

def polymin(P_bar = None, Res = None):
    '''
    Finds the minimum residual temperature using a 2nd degree polynomial.
    '''
    arr = np.sort(Res)
    Ind = np.where(Res == arr[0])[0][0]

    if P_bar[Ind] == np.nanmax(P_bar):
        p = np.array([0,1,0])
        p_min = np.array([P_bar[Ind]])
    elif P_bar[Ind] == np.nanmin(P_bar):
        p = np.array([0,1,0])
        p_min = np.array([P_bar[Ind]])
    else:
        p = np.polyfit(P_bar[np.array([Ind-1, Ind, Ind+1])],Res[np.array([Ind-1, Ind, Ind+1])],2)

        x = np.linspace(np.nanmin(P_bar),np.nanmax(P_bar),501)
        y = p[0]*x**2 + p[1]*x + p[2]

        p_min = x[np.where(y == np.nanmin(y))]

    return p, p_min

def satTemperature(q, index, *, Model = None, comp = None, phases = None, T_initial_C = None, T_step_C = None, dt_C = None, P_bar = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):
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

    T_initial_C: float
        Initial guess for the liquidus temperature.

    T_step_C: float
        The temperature drop at each step of the calculation.

    dt_C: float
        Total temperature drop allowed duringmodel runs.

    P_bar: float
         Specifies the pressure of calculation (bar).

    Returns:
    ----------
    Results: Dict
        Dict containing a series of floats that represent the saturation temperature and residual temperature for each calculation.

    index: int
        index of the calculation

    '''

    Results = {}
    if "MELTS" in Model:

        try:
            Results = phaseSat_MELTS(Model = Model, comp = comp, phases = phases, T_initial_C = T_initial_C, T_step_C = T_step_C, dt_C = dt_C, P_bar = P_bar, H2O_Liq = H2O_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
        except:
            Results = {phases[0]: np.nan, phases[1]: np.nan, phases[2]: np.nan, 'T_Liq': np.nan, 'H2O_melt': np.nan}
            if len(phases) == 2:
                del Results[phases[2]]

        if len(phases) == 3:
            Res = ['3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]
            for R in Res:
                Results[R] = np.nan

            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]) and ~np.isnan(Results[phases[2]]):
                Results['3 Phase Saturation'] = np.nanmax(np.array([abs(Results[phases[0]] - Results[phases[1]]), abs(Results[phases[0]] - Results[phases[2]]), abs(Results[phases[1]] - Results[phases[2]])]))
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]):
                Results[phases[0] + ' - ' + phases[1]] = abs(Results[phases[0]] - Results[phases[1]])
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[2]]):
                Results[phases[0] + ' - ' + phases[2]] = abs(Results[phases[0]] - Results[phases[2]])
            if ~np.isnan(Results[phases[1]]) and ~np.isnan(Results[phases[2]]):
                Results[phases[1] + ' - ' + phases[2]] = abs(Results[phases[1]] - Results[phases[2]])
        else:
            Results[phases[0] + ' - ' + phases[1]] = np.nan
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]):
                Results[phases[0] + ' - ' + phases[1]] = abs(Results[phases[0]] - Results[phases[1]])

        q.put([Results, index])
        return

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        Results = {phases[0]: np.nan, phases[1]: np.nan, phases[2]: np.nan, 'T_Liq': np.nan, 'H2O_melt': np.nan}
        if len(phases) == 2:
            del Results[phases[2]]

        #try:
        Result = MM.path(comp = comp, phases = phases, T_min_C = dt_C, dt_C = T_step_C, P_bar = P_bar, find_liquidus = True)
        #Result = stich(Result, Model = Model)

        for i in range(len(phases)):
            try:
                Results[phases[i]] = Result['Conditions']['T_C'][Result[phases[i]+'_prop']['mass'] > 0.0].values[0]
                print(Results[phases[i]])
            except:
                Results[phases[i]] = np.nan

        Results['T_Liq'] = Result['Conditions']['T_C'].values[0]
        Results['H2O_melt'] = Result['liq']['H2O'].values[0]

        if len(phases) == 3:
            Res = ['3 Phase Saturation', phases[0] + ' - ' + phases[1], phases[0] + ' - ' + phases[2], phases[1] + ' - ' + phases[2]]
            for R in Res:
                Results[R] = np.nan

            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]) and ~np.isnan(Results[phases[2]]):
                Results['3 Phase Saturation'] = np.nanmax(np.array([abs(Results[phases[0]] - Results[phases[1]]), abs(Results[phases[0]] - Results[phases[2]]), abs(Results[phases[1]] - Results[phases[2]])]))
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]):
                Results[phases[0] + ' - ' + phases[1]] = abs(Results[phases[0]] - Results[phases[1]])
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[2]]):
                Results[phases[0] + ' - ' + phases[2]] = abs(Results[phases[0]] - Results[phases[2]])
            if ~np.isnan(Results[phases[1]]) and ~np.isnan(Results[phases[2]]):
                Results[phases[1] + ' - ' + phases[2]] = abs(Results[phases[1]] - Results[phases[2]])
        else:
            Results[phases[0] + ' - ' + phases[1]] = np.nan
            if ~np.isnan(Results[phases[0]]) and ~np.isnan(Results[phases[1]]):
                Results[phases[0] + ' - ' + phases[1]] = abs(Results[phases[0]] - Results[phases[1]])
        
        q.put([Results, index])
        #except:
        #    q.put([Results, index])
        return


