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

def AdiabaticDecompressionMelting(cores = multiprocessing.cpu_count(), 
                                  Model = "pMELTS", bulk = "KLB-1", comp_lith_1 = None, 
                                  comp_lith_2 = None, comp_lith_3 = None, Tp_C = 1350, Tp_Method = "pyMelt",
                                  P_start_bar = 30000, P_end_bar = 2000, dp_bar = 200, 
                                  P_path_bar = None, Frac = False, prop = None, 
                                  fO2_buffer = None, fO2_offset = None, Fe3Fet = None, MELTS_filter = True):
    """
    Perform adiabatic decompression melting calculations using MELTS, MAGEMin, or pyMelt.

    Simulates mantle melting along an adiabatic upwelling path (e.g., ridge or plume) with user-defined 
    starting potential temperature, pressure range, and step size. Supports single-lithology mantle sources 
    (e.g., KLB-1) at present, with expansion to multi-lithology systems in development.

    Parameters
    ----------
    cores : int, optional
        Number of CPU cores to use for multiprocessing. Defaults to total available.
    Model : str, optional
        Thermodynamic model. MELTS variants: "MELTSv1.0.2", "MELTSv1.1.0", "MELTSv1.2.0", "pMELTS";
        or MAGEMin: "Green2025", "Weller2024". Alternatively calculations can be performed using 
        pyMelt (Matthews et al. 2020): "pyMelt"
    bulk : dict ot str, optional
        Bulk composition name or composition dictionary.
        Default is "KLB-1".
    Tp_C : float or np.ndarray, optional
        Mantle potential temperature(s) in °C. Default is 1350.
    Tp_Method : str, optional
        Method to calculate the starting pressure for adiabatic melting. Default is "pyMelt".
    P_start_bar, P_end_bar, dp_bar : float or array, optional
        Starting, ending, and step size pressures (in bar) for adiabatic decompression. 
        Defaults: 30000, 2000, and 200, respectively.
    P_path_bar : np.ndarray, optional
        User-specified pressure path (in bar). If given, overrides `P_start_bar`, `P_end_bar`, and `dp_bar`.
    fO2_buffer : {"FMQ", "NNO"}, optional
        Redox buffer for constraining oxygen fugacity.
    fO2_offset : float, optional
        Offset (log units) from the chosen fO2 buffer.
    Fe3Fet : float, optional
        Initial Fe³⁺/ΣFe ratio for the bulk composition. If None, values is taken from the "bulk" variable or set according to fO2 buffer positions.
    MELTS_filter : bool, default=True
        If True, filters oxide components to avoid issues in MELTS calculations (e.g., K2O content set to 0.0).

    Returns
    -------
    Results : dict
        Dictionary containing DataFrames for the system and phase compositions and properties.

    Notes
    -----
    - Currently limited to single-lithology melting.
    - Normalizes output mass so that total initial mass = 1.

    Examples
    --------
    Run a single adiabatic decompression path from 3.0 GPa to 0.2 GPa:

    >>> results = AdiabaticDecompressionMelting(Model="pMELTS", bulk="KLB-1",
    ...                                        Tp_C=1350, P_start_bar=30000,
    ...                                        P_end_bar=2000, dp_bar=200)

    Run with an explicit pressure path:

    >>> import numpy as np
    >>> P_path = np.linspace(30000, 2000, 20)
    >>> results = AdiabaticDecompressionMelting(Model="pMELTS", comp_lith_1="KLB-1",
    ...                                        P_path_bar=P_path, Tp_C=1400)
    """
    
    Tp_C        = to_float(Tp_C)

    P_path_bar = to_float(P_path_bar)
    P_start_bar= to_float(P_start_bar)
    P_end_bar  = to_float(P_end_bar)
    dp_bar     = to_float(dp_bar)

    Fe3Fet = to_float(Fe3Fet)
    fO2_offset = to_float(fO2_offset)

    if fO2_buffer is not None:
        if fO2_buffer != "NNO":
            if fO2_buffer != "FMQ":
                raise Warning("fO2 buffer specified is not an allowed input. This argument can only be 'FMQ' or 'NNO' \n if you want to offset from these buffers use the 'fO2_offset' argument.")

    if "MELTS" not in Model:
        if fO2_buffer == "FMQ":
            fO2_buffer = "qfm"
        if fO2_buffer == "NNO":
            fO2_buffer = "nno"

    # if Tp_Method == "pyMelt":
    #     try:
    #         import pyMelt as m
    #         Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
    #                     'KG1': m.lithologies.matthews.kg1(),
    #                     'G2': m.lithologies.matthews.eclogite(),
    #                     'hz': m.lithologies.shorttle.harzburgite()}
    #     except ImportError:
    #         raise RuntimeError('You havent installed pyMelt or there is an error when importing pyMelt. pyMelt is currently required to estimate the starting point for the melting calculations.')

    if bulk is not None and comp_lith_1 is None:
        comp_lith_1 = bulk

    comp_1 = comp_check(comp_lith_1, Model, MELTS_filter, Fe3Fet)

    # place holders for when code is expanded to account for multi-lithology mantle
    if comp_lith_2 is not None:
        comp_2 = comp_check(comp_lith_2, Model, MELTS_filter, Fe3Fet)
    else:
        comp_2 = None

    if comp_lith_3 is not None:
        comp_3 = comp_check(comp_lith_3, Model, MELTS_filter, Fe3Fet)
    else:
        comp_3 = None


    # At present calculations only work for a single simulation - this represents a placeholder for when the code is expanded to account for multiple simulations
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
                              'Tp_C': Tp_C, 'Tp_Method': Tp_Method, 'P_path_bar': P_path_bar, 
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

def AdiabaticMelt(q, index, *, Model = None, comp_1 = None, comp_2 = None, comp_3 = None, Tp_Method = "pyMelt",
                  Tp_C = None, P_start_bar = None, P_end_bar = None, dp_bar = None, P_path_bar = None, 
                  Frac = None, fO2_buffer = None, fO2_offset = None, prop = None):
    '''
    Melting calculations to be performed in parallel.

    '''
    Results = {}
    if "MELTS" in Model:
        try:
            Results = AdiabaticDecompressionMelting_MELTS(Model = Model, comp = comp_1, Tp_C = Tp_C, Tp_Method = "pyMelt",
                                                          P_path_bar = P_path_bar, P_start_bar = P_start_bar, P_end_bar = P_end_bar, dp_bar = dp_bar, 
                                                          fO2_buffer = fO2_buffer, fO2_offset = fO2_offset)
            q.put([Results, index])
        except:
            q.put([])

        return
    
    elif Model == "pyMelt":
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

    else:
        try:
            import pyMelt as m
            Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
                        'KG1': m.lithologies.matthews.kg1(),
                        'G2': m.lithologies.matthews.eclogite(),
                        'hz': m.lithologies.shorttle.harzburgite()}
        except ImportError:
            raise RuntimeError('You havent installed pyMelt or there is an error when importing pyMelt. pyMelt is currently required to estimate the starting point for the melting calculations.')

        if Tp_Method == "pyMelt":
            lz = m.lithologies.matthews.klb1()
            mantle = m.mantle([lz], [1], ['Lz'])
            T_start_C = mantle.adiabat(P_start_bar/10000.0, Tp_C)
        else:
            T_start_C = None

        from juliacall import Main as jl, convert as jlconvert

        jl.seval("using MAGEMinCalc")

        comp_1['O'] = comp_1['Fe3Fet_Liq']*(((159.59/2)/71.844)*comp_1['FeOt_Liq'] - comp_1['FeOt_Liq'])

        if type(comp_1) == dict:
            comp_julia = jl.seval("Dict")(comp_1) 
        else:
            comp_new = comp_1.loc[0].to_dict()
            comp_julia = jl.seval("Dict")(comp_new)

        Output_jl = jl.MAGEMinCalc.AdiabaticDecompressionMelting(comp = comp_julia, P_start_kbar = P_start_bar/1000.0, 
                                                              P_end_kbar = P_end_bar/1000.0, dp_kbar = dp_bar/1000.0,
                                                              T_start_C = T_start_C, fo2_buffer = fO2_buffer, fo2_offset = fO2_offset, Model = Model)
        Results = dict(Output_jl)
        # Results = stich(Results, Model = Model)
        
        q.put([Results, index])
        return
        

