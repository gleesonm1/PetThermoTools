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

def multi_iso_crystallise(cores = None, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_bar = None, isochoric = None, find_liquidus = None):
    '''
    Carry out multiple crystallisation calculations in parallel. Allows isobaric or isochoric calculations to be performed. All temperature inputs/outputs are reported in degrees celcius and pressure is reported in bars.

    Parameters:
    ----------
    cores: int
        number of processes to run in parallel. Default is 4.

    Model: string
        "MELTS" or "Holland". Dictates whether MELTS or MAGEMin calculations are performed. Default "MELTS".
        Version of melts can be specified "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", or "pMELTS". Default "v.1.0.2".

    comp: Dict or pd.DataFrame
        Initial compositon for calculations. If type == Dict, the same initial composition will be used in all calculations.

    Frac_solid:
        If True, solid phases will be removed from the system at the end of each crystallisation step. Default False.

    Frac_fluid:
        If True, fluid phases will be removed from the system at the end of each crystallisation step. Default False.

    T_start_C: float
        Initial temperature used for crystallisation calculations.

    T_end_C: float
        Final temperature in crystallisation calculations. Default = 750.

    dt_C: float
        Temperature increment during crystallisation calculations. Default = 2.

    P_bar: float or np.ndarray
        Single pressure or an array of pressures equal to the number of calculations to be performed. Specifies the pressure of calculation (bar).

    isochoric: True/False
        If True, the volume of the system will be held constant instead of the pressure. Default is False.

    find_liquidus: True/False
        If True, the calculations will start with a search for the melt liquidus temperature. Default is False.

    Returns:
    ----------
    Results: Dict
        Dictionary where each entry represents the results of a single calculation. Within the dictionary each single calculation is reported as a series of pandas DataFrames, displaying the composition and thermodynamic properties of each phase.

    '''

    if Model is None:
        Model == "MELTSv1.0.2"

    if cores is None:
        cores = 4

    if dt_C is None and T_path_C is None:
        dt_C = 2

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
                p = Process(target = iso_crystallise, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp, 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_path_bar': P_bar[i], 'isochoric': isochoric, 'find_liquidus': find_liquidus})
            elif type(P_bar) != np.ndarray:
                p = Process(target = iso_crystallise, args = (q, i),
                            kwargs = {'Model': Model, 'comp': comp.loc[i].to_dict(), 'Frac_solid': Frac_solid, 'Frac_fluid': Frac_fluid,
                                    'T_start_C': T_start_C, 'T_end_C': T_end_C, 'dt_C': dt_C,
                                    'P_path_bar': P_bar, 'isochoric': isochoric, 'find_liquidus': find_liquidus})
            else:
                p = Process(target = iso_crystallise, args = (q, i),
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

def iso_crystallise(q, index, *, Model = None, comp = None, Frac_solid = None, Frac_fluid = None, T_start_C = None, T_end_C = None, dt_C = None, P_path_bar = None, isochoric = None, find_liquidus = None):

    Results = {}

    if "MELTS" in Model:
        try:
            Results = crystallise_MELTS(Model = Model, comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, isochoric = isochoric, find_liquidus = find_liquidus)
            q.put([Results, index])
            return
        except:
            q.put([Results, index])
            return

    if Model == "Holland":
        try:
            Results = crystallise_holland(comp = comp, Frac_solid = Frac_solid, Frac_fluid = Frac_fluid, T_start_C = T_start_C, T_end_C = T_end_C, dt_C = dt_C, P_path_bar = P_path_bar, isochoric = isochoric, find_liquidus = find_liquidus)
            q.put([Results, index])
            return
        except:
            q.put([Results, index])
            return

