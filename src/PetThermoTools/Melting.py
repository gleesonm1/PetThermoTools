import numpy as np
import pandas as pd
from PetThermoTools.GenFuncs import *
from PetThermoTools.Plotting import *
from PetThermoTools.MELTS import *
from PetThermoTools.Compositions import *
import multiprocessing
from multiprocessing import Queue
from multiprocessing import Process
import time
import sys
from tqdm.notebook import tqdm, trange
# import pyMelt as m

def comp_check(comp_lith, Model, MELTS_filter, Fe3Fet):
    if type(comp_lith) == str:
        if Model != "pyMelt":
            comp = Compositions[comp_lith]
        else:
            comp = comp_lith
    else:
        comp = comp_lith.copy()

    # if comp is entered as a pandas series, it must first be converted to a dict
    if Model != "pyMelt":
        if type(comp) == pd.core.series.Series:
            comp = comp.to_dict()

        comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet)

    if "MELTS" in Model and MELTS_filter == True:
        if type(comp) == pd.core.frame.DataFrame:
            comp['K2O_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['P2O5_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['H2O_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['CO2_Liq'] = np.zeros(len(comp['SiO2_Liq']))
        else:
            comp['K2O_Liq'] = 0
            comp['P2O5_Liq'] = 0
            comp['H2O_Liq'] = 0
            comp['CO2_Liq'] = 0
    
    return comp

def AdiabaticDecompressionMelting(cores = multiprocessing.cpu_count(), 
                                  Model = "pMELTS", bulk = "KLB-1", comp_lith_1 = None, 
                                  comp_lith_2 = None, comp_lith_3 = None, Tp_C = 1350, Tp_Method = "pyMelt",
                                  P_start_bar = 30000, P_end_bar = 5000, dp_bar = 200, 
                                  P_path_bar = None, Frac = False, prop = None, 
                                  fO2_buffer = None, fO2_offset = None, Fe3Fet = None, MELTS_filter = True):
    
    try:
        import pyMelt as m
        Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
                    'KG1': m.lithologies.matthews.kg1(),
                    'G2': m.lithologies.matthews.eclogite(),
                    'hz': m.lithologies.shorttle.harzburgite()}
    except ImportError:
        raise RuntimeError('You havent installed pyMelt or there is an error when importing pyMelt. pyMelt is currently required to estimate the starting point for the melting calculations.')

    if bulk is not None and comp_lith_1 is None:
        comp_lith_1 = bulk

    comp_1 = comp_check(comp_lith_1, Model, MELTS_filter, Fe3Fet)

    if comp_lith_2 is not None:
        comp_2 = comp_check(comp_lith_2, Model, MELTS_filter, Fe3Fet)
    else:
        comp_2 = None

    if comp_lith_3 is not None:
        comp_3 = comp_check(comp_lith_3, Model, MELTS_filter, Fe3Fet)
    else:
        comp_3 = None

    One = 0
    if Model != "pyMelt":
        if type(comp_1) == pd.core.frame.DataFrame: # simplest scenario - one calculation per bulk composition imported
            A = len(comp_1['SiO2_Liq'])//cores
            B = len(comp_1['SiO2_Liq'])%cores
        if P_start_bar is not None and type(P_start_bar) == np.ndarray: #one calculation per P loaded
            A = len(P_start_bar)//cores
            B = len(P_start_bar)%cores
        elif Tp_C is not None and type(Tp_C) == np.ndarray: # one calculation per T loaded.
            A = len(Tp_C)//cores
            B = len(Tp_C)%cores

        else: # just one calculation
            One = 1
            A = 1
            B = 0
    else:
        if P_start_bar is not None and type(P_start_bar) == np.ndarray: #one calculation per P loaded
            A = len(P_start_bar)//cores
            B = len(P_start_bar)%cores
        elif Tp_C is not None and type(Tp_C) == np.ndarray: # one calculation per T loaded.
            A = len(Tp_C)//cores
            B = len(Tp_C)%cores

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
        p = Process(target = AdiabaticMelt, args = (q, 1),
                    kwargs = {'Model': Model, 'comp_1': comp_1, 'comp_2': comp_2, 'comp_3': comp_3,
                              'Tp_C': Tp_C, 'P_path_bar': P_path_bar, 
                              'P_start_bar': P_start_bar, 'P_end_bar': P_end_bar, 'dp_bar': dp_bar,
                              'fO2_buffer': fO2_buffer, 'fO2_offset': fO2_offset, 'Frac': Frac, 'prop': prop})

        p.start()
        try:
            ret = q.get(timeout = 600)
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
            if Model != "pyMelt":
                Results = stich(Results, Model = Model)

                # make mass relative to 1 at the start of the model
                Tot = Results['Mass'].sum(axis = 1)[0]
                Results['Mass'] = Results['Mass']/Tot
                Results['All'].loc[:,Results['All'].columns.str.contains('mass')] = Results['All'].loc[:,Results['All'].columns.str.contains('mass')]/Tot
            return Results
        else:
            Results = {}
            return Results

    return Results

def AdiabaticMelt(q, index, *, Model = None, comp_1 = None, comp_2 = None, comp_3 = None, Tp_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, P_path_bar = None, Frac = None, fO2_buffer = None, fO2_offset = None, prop = None):
    '''
    Melting calculations to be performed in parallel.

    '''
    Results = {}
    if "MELTS" in Model:
        try:
            Results = AdiabaticDecompressionMelting_MELTS(Model = Model, comp = comp_1, Tp_C = Tp_C, P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
            q.put([Results, index])
        except:
            q.put([])

        return

    if Model == "Holland":
        import pyMAGEMINcalc as MM
        Results = MM.AdiabaticDecompressionMelting(comp = comp_1, T_p_C = Tp_C, P_start_kbar = P_start_bar/1000, P_end_kbar = P_end_bar/1000, dp_kbar = dp_bar/1000, Frac = 0)
        q.put([Results, index])
        return
    
    if Model == "pyMelt":
        try:
            import pyMelt as m
            Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
                        'KG1': m.lithologies.matthews.kg1(),
                        'G2': m.lithologies.matthews.eclogite(),
                        'hz': m.lithologies.shorttle.harzburgite()}
        except ImportError:
            raise RuntimeError('You havent installed pyMelt or there is an error when importing pyMelt. pyMelt is currently required to estimate the starting point for the melting calculations.')

        lith_1 = Lithologies[comp_1]

        if prop is None:
            mantle = m.mantle([lith_1],[1],[comp_1])
            column = mantle.adiabaticMelt(Tp_C, Pstart = P_start_bar/10000, Pend = P_end_bar/10000, dP = -dp_bar/10000)
            Results['All'] = pd.DataFrame(data = np.zeros((len(column.P), 3)), columns = ['T_C', 'P_bar', 'mass_Liq'])
            Results['All']['T_C'] = column.T
            Results['All']['P_bar'] = column.P*10000
            Results['All']['mass_Liq'] = column.F
        
        if prop is not None:
            if type(prop) == float or type(prop) == int:
                mantle = m.mantle([lith_1],[1],[comp_1])
                column = mantle.adiabaticMelt(Tp_C, Pstart = P_start_bar/10000, Pend = P_end_bar/10000, dP = -dp_bar/10000)
                Results['All'] = pd.DataFrame(data = np.zeros((len(column.P), 3)), columns = ['T_C', 'P_bar', 'mass_Liq'])
                Results['All']['T_C'] = column.T
                Results['All']['P_bar'] = column.P*10000
                Results['All']['mass_Liq'] = column.F
            else:
                if len(prop) == 2:
                    lith_2 = Lithologies[comp_2]

                    mantle = m.mantle([lith_1, lith_2], [prop[0], prop[1]], [comp_1, comp_2])
                    column = mantle.adiabaticMelt(Tp_C, Pstart = P_start_bar/10000, Pend = P_end_bar/10000, dP = -dp_bar/10000)
                    Results['All'] = pd.DataFrame(data = np.zeros((len(column.P), 5)), columns = ['T_C', 'P_bar', 'mass_Liq_tot',
                                                                                                  'mass_Liq_'+comp_1,
                                                                                                  'mass_Liq_'+comp_2])
                    Results['All']['T_C'] = column.T
                    Results['All']['P_bar'] = column.P*10000
                    Results['All']['mass_Liq_tot'] = column.F
                    Results['All']['mass_Liq_'+comp_1] = column.lithologies[comp_1].F
                    Results['All']['mass_Liq_'+comp_2] = column.lithologies[comp_2].F

                if len(prop) == 3:
                    lith_2 = Lithologies[comp_2]
                    lith_3 = Lithologies[comp_3]
                        
                    mantle = m.mantle([lith_1, lith_2, lith_3], [prop[0], prop[1], prop[2]], [comp_1, comp_2, comp_3])
                    column = mantle.adiabaticMelt(Tp_C, Pstart = P_start_bar/10000, Pend = P_end_bar/10000, dP = -dp_bar/10000)
                    Results['All'] = pd.DataFrame(data = np.zeros((len(column.P), 6)), columns = ['T_C', 'P_bar', 'mass_Liq_tot',
                                                                                                  'mass_Liq_'+comp_1,
                                                                                                  'mass_Liq_'+comp_2,
                                                                                                  'mass_Liq_'+comp_3])
                    Results['All']['T_C'] = column.T
                    Results['All']['P_bar'] = column.P*10000
                    Results['All']['mass_Liq_tot'] = column.F
                    Results['All']['mass_Liq_'+comp_1] = column.lithologies[comp_1].F
                    Results['All']['mass_Liq_'+comp_2] = column.lithologies[comp_2].F
                    Results['All']['mass_Liq_'+comp_3] = column.lithologies[comp_3].F
        
        q.put([Results, index])
        return
        

