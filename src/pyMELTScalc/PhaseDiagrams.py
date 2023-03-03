import numpy as np
import pandas as pd
from pyMELTScalc.GenFuncs import *
from pyMELTScalc.Plotting import *
from pyMELTScalc.Path import *
from pyMELTScalc.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
from tqdm.notebook import tqdm, trange
import pickle

def phaseDiagram_calc(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None, i_max = 25):

    comp = bulk.copy()

    if cores is None:
        cores = multiprocessing.cpu_count()

    if Model is None:
        Model = "MELTSv1.0.2"

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

    A = len(P_flat)//cores + 1

    # B = len(P_flat) % cores
    #
    # if A > 0:
    #     Group = np.zeros(A) + cores
    #     if B > 0:
    #         Group = np.append(Group, B)
    # else:
    #     Group = np.array([B])

    c = 0
    j = 0
    while len(T_flat)>1:
        print('Attempt ' + str(j))
        if j > i_max:
            break

        qs = []
        q = Queue()
        ps = []
        j = j + 1
        for i in range(cores):
            T_path_C = T_flat[i*A:(i+1)*A]
            P_path_bar = P_flat[i*A:(i+1)*A]

            if len(T_path_C) > 100:
                T_path_C = T_path_C[:99]
                P_path_bar = P_path_bar[:99]

            p = Process(target = path, args = (q,i), kwargs = {'Model': Model, 'comp': comp, 'T_path_C': T_path_C, 'P_path_bar': P_path_bar, 'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset})

            ps.append(p)
            p.start()

        TIMEOUT = 300 #+ #0.5*len(T_flat)
        start = time.time()
        for p in ps:
            if p.is_alive():
                while time.time() - start <= TIMEOUT:
                    if not p.is_alive():
                        p.terminate()
                        p.join()
                        break
                    time.sleep(.1)
                else:
                    p.terminate()
                    p.join()
            else:
                p.terminate()
                p.join()

            try:
                ret = q.get(timeout = 2)
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
                Combined = Results[i]['All'].copy()
                c = 1
            else:
                Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)

        flat = np.round(np.array([T.flatten(), P.flatten()]).T,2)
        Res_flat = np.round(np.array([Combined['T_C'], Combined['P_bar']]).T,2)
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

    return Combined


def phaseDiagram_eq(cores = None, Model = None, bulk = None, T_C = None, P_bar = None, T_min_C = None, T_max_C = None, T_num = None, P_min_bar = None, P_max_bar = None, P_num = None, Fe3Fet_Liq = None, H2O_Liq = None, fO2_buffer = None, fO2_offset = None):

    comp = bulk.copy()

    if cores is None:
        cores = multiprocessing.cpu_count()

    if Model is None:
        Model = "MELTSv1.0.2"

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

        TIMEOUT = 30 #+ #0.5*len(T_flat)
        start = time.time()
        for p in ps:
            if TIMEOUT  - (time.time() - start) > 5:
                try:
                    ret = q.get(timeout = TIMEOUT - (time.time() - start))
                except:
                    ret = []

            else:
                try:
                    ret = q.get(timeout = 5)
                except:
                    ret = []

            qs.append(ret)

        # TIMEOUT = 5 #+ #0.5*len(T_flat)
        # start = time.time()
        # for p in ps:
        #     if p.is_alive():
        #         while time.time() - start <= TIMEOUT:
        #             if not p.is_alive():
        #                 p.terminate()
        #                 p.join()
        #                 break
        #             time.sleep(.1)
        #         else:
        #             p.terminate()
        #             p.join()
        #     else:
        #         p.terminate()
        #         p.join()

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
            Combined = Results[i]['All'].copy()
            c = 1
        else:
            Combined = pd.concat([Combined, Results[i]['All']], axis = 0, ignore_index = True)

    return Combined


def equilibrate(q, index,*, Model = None, P_bar = None, T_C = None, comp = None, fO2_buffer = None, fO2_offset = None):

    Results = equilibrate_MELTS(Model = Model, P_bar = P_bar, T_C = T_C, comp = comp, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
    Phases = {}
    # for pp in PhaseList:
    #     Phases[pp] = pd.DataFrame(columns = ['T_C', 'P_bar', 'SiO2', 'TiO2', 'Al2O3', 'FeOt', 'Cr2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'mass', 'rho'], data = np.array([T_C, P_bar, PhaseComp[pp][0], PhaseComp[pp][1], PhaseComp[pp][2], (71.844/(159.69/2))*PhaseComp[pp][3]+PhaseComp[pp][5], PhaseComp[pp][4], PhaseComp[pp][6], PhaseComp[pp][7], PhaseComp[pp][10], PhaseComp[pp][11], PhaseComp[pp][12], PhaseComp[pp][13], PhaseComp[pp][14], PhaseComp[pp][15], PhaseProp[pp]['mass'], PhaseProp[pp]['rho']]).reshape(1,17))

    q.put([Results, index])
    return
    # except:
    #     q.put([PhaseList, PhaseComp, PhaseProp, index])
    #     return