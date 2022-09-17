import numpy as np
import pandas as pd
import sys
import os
import pickle
import matplotlib.pyplot as plt
import subprocess
import glob
from subprocess import Popen, PIPE
import Thermobar as pt
from pyMELTScalc.Barom import *
from pyMELTScalc.Liq import *
from pyMELTScalc.Holland import *
from pyMELTScalc.MELTS import *
import asyncio
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time

def multi_crystallise(cores = None, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, isochoric = None, find_liquidus = None):

    if Model is None:
        Model == "MELTS"

    if cores is None:
        cores = 4

    if type(P_bar) == np.ndarray:
        A = len(P_bar)//cores
        B = len(P_bar) % cores
    else:
        A = len(comp['SiO2'])//cores
        B = len(comp['SiO2']) % cores

    Group = np.zeros(A) + cores
    Group = np.append(Group, B)

    qs = []
    q = Queue()

    for j in range(len(Group)):
        ps = []

        for i in range(int(cores*j), int(cores*j + Group[j])):
            if type(comp) == dict:
                p = Process(target = crystallise, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_path_bar': P_bar[i], 'isochoric': isochoric, 'find_liquidus': find_liquidus})
            elif type(P_bar) != np.ndarray:
                p = Process(target = crystallise, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_path_bar': P_bar, 'isochoric': isochoric, 'find_liquidus': find_liquidus})
            else:
                p = Process(target = crystallise, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_path_bar': P_bar[i], 'isochoric': isochoric, 'find_liquidus': find_liquidus})

            ps.append(p)
            p.start()

        for p in ps:
            try:
                ret = q.get(timeout = 180)
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

    if type(P_bar) == np.ndarray:
        Results = {}
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Res, index = qs[i]
                Results['P = ' + str(round(P_bar[index],2)) + ' bars'] = Res
    else:
        Results = {}
        for i in range(len(qs)):
            if len(qs[i]) > 0:
                Res, index = qs[i]
                Results['index = ' + str(index)] = Res

    return Results

def crystallise(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_path_C = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, P_start_bar = None, P_end_bar = None, dp_bar = None, isochoric = None, find_liquidus = None):

    Results = {}

    if Model == "MELTS":
        try:
            Results = crystallise_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, isochoric = isochoric, find_liquidus = find_liquidus)
            q.put([Results, index])
            return
        except:
            q.put([Results, index])
            return

    if Model == "Holland":
        #try:
        Results = crystallise_holland(comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_path_C = T_path_C, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, isochoric = isochoric, find_liquidus = find_liquidus)
        q.put([Results, index])
        return
        #except:
        #    q.put([Results, index])
        #    return

