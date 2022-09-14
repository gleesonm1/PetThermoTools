import numpy as np
import julia
import time
from pyMELTScalc.Liq import *
from pyMELTScalc.Crystallise import *
from pyMELTScalc.Barom import *
from pyMELTScalc.MELTS import *
#from julia import Main
from julia import MAGEMinCalc
from os import path
from pathlib import Path
import pickle

#with open('directory', 'rb') as f:
#    this_directory = pickle.load(f)

#file_dir = Path(__file__).parent
#Main.using("MAGEMinCalc")

#Main.include("functions.jl")

def findLiq_holland(P, T_initial, bulk):
    T_Liq = 0

    T = T_initial

    bulk_in = bulk.copy()

    Liq = ["liq","fl"]

    start = time.time()

    PhaseList = MAGEMinCalc.satPhase(P, T, bulk)

    i = set.intersection(set(Liq),set(PhaseList))

    Step = np.array([3,1,0.1])
    for k in range(len(Step)):
        if len(i) == len(PhaseList):
            while len(i) == len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T - Step[k]
                    Ret = MAGEMinCalc.satPhase(P, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

        if len(i) < len(PhaseList):
            while len(i) < len(PhaseList):
                bulk = bulk_in.copy()
                if time.time() - start > 60:
                    return T_Liq
                else:
                    T = T + Step[k]
                    Ret = MAGEMinCalc.satPhase(P, T, bulk)
                    PhaseList = Ret['Phase']
                    i = set.intersection(set(Liq),set(PhaseList))

    if "liq" in Ret['Phase']:
        if Ret['Liq_Frac'] < 0.9:
            return T_Liq
        else:
            T_Liq = T
            return T_Liq
    else:
        return T_Liq