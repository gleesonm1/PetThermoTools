import numpy as np
import pandas as pd
from pyMELTScalc.Barom import *
from pyMELTScalc.Crystallise import *
from pyMELTScalc.Holland import *
from pyMELTScalc.MELTS import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process

def findLiq(q, Model, P, T_initial, bulk, index):

    T_Liq = 0
    T_in = T_initial
    H2O_Melt = 0

    if Model == "MELTS":
        try:
            T_Liq, H2O_Melt = findLiq_MELTS(P, Model, T_initial, bulk)
            q.put([T_Liq, H2O_Melt, index, T_in])
        except:
            q.put([T_Liq, H2O_Melt, index])

    if Model == "Holland":
        T_Liq = findLiq_holland(P, T_initial, bulk)
        q.put([T_Liq, H2O_Melt, index, T_in])

