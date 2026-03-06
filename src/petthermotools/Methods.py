### functions to allow specific methods in MELTS
import numpy as np 
import pandas as pd
from petthermotools.GenFuncs import *
from petthermotools.Liq import *
from petthermotools.MELTS import *

def spinel_Fe3(spl):
    oxides = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MnO', 'MgO', 'CaO', 'Cr2O3', 'V2O5', 'ZnO']

    molar_mass = [60.08, 79.87, 101.96, 71.84, 70.94, 40.30, 56.08, 151.99, 149.88, 81.38]
    mm = pd.Series(data = molar_mass, index = oxides)
    
    cat = pd.Series(data = [1,1,2,1,1,1,1,2,2,1], index = oxides)
    X = ((spl.div(mm))).div(((spl.div(mm))).sum(axis=1),axis=0) ## Is the issue here?
    ox = pd.Series(data = [2,2,3,1,1,1,1,3,3,1], index = oxides)

    at_prop = spl.div(mm).multiply(ox)
    ox_no = at_prop.multiply(32/at_prop.sum(axis = 1), axis = 0)
    cat_no = ox_no.multiply(cat/ox)
    cat_no_norm = 24*cat_no.div(cat_no.sum(axis=1),axis=0)

    cat_no_norm['Fe2O3'] = (2*32)*(1-(24/cat_no.sum(axis=1))).values
    cat_no_norm["FeO"] = cat_no_norm['FeOt'].values - cat_no_norm['Fe2O3']
    cat_no_norm.loc[cat_no_norm['Fe2O3']<0.0, "FeO"] = cat_no_norm.loc[cat_no_norm['Fe2O3']<0.0, "FeOt"]
    cat_no_norm.loc[cat_no_norm['Fe2O3']<0.0, "Fe2O3"] = 0.0
    cat_no_norm = cat_no_norm.drop(columns = "FeOt")

    Fe3Fe2 = cat_no_norm['Fe2O3']/(cat_no_norm['FeO'] + cat_no_norm['Fe2O3'])

    spl['Fe2O3'] = 1.1113*(Fe3Fe2*spl['FeOt'])
    spl['FeO'] = (1-Fe3Fe2)*spl['FeOt']
    spl = spl.drop(columns = "FeOt")

    return spl, X 

def MELTS_OSaS(ol = None, spl = None, liq = None, P_bar = None, T_C = None):
    spl, X = spinel_Fe3(spl)
    
    liq_data = equilibrate_multi(Model = "MELTSv1.0.2", Suppress="All", bulk=liq,
                      P_bar = P_bar, T_C = T_C, fO2_buffer = "FMQ")
    ol_data = supCalc(phase="olivine", bulk=ol,
                      P_bar = P_bar, T_C = T_C)
    spl_data = supCalc(phase="spinel", bulk=spl,
                      P_bar = P_bar, T_C = T_C)
    
    muSiO2 = liq_data['liquid1_prop']['mu_sio2_Liq']
    muFe2SiO4=ol_data['Thermodynamics']['mu_fe2sio4']
    muFe3O4 = spl_data['Thermodynamics']['mu_fe3o4']

    ## correction for muFe3O4 taken from Bell et al. 2025
    muFe3O4_correction  = (-88317.116658736 + 0.968070202*muFe3O4+ 
                       66581.745702629*X['Cr2O3'] + 75747.288514812*X['MgO'] + 53625.743547049*X['FeOt'])
    
    muO2 = 3*muSiO2+2*muFe3O4_correction-3*muFe2SiO4

    ## calculations for hs, ss and mu0_O2 taken from Bell et al. 2025
    t  =  (T_C + 273.15).astype(float)
    tr = 298.15
    hs = 23.10248*(t-tr) + 2.0*804.8876*(np.sqrt(t)-np.sqrt(tr)) \
        - 1762835.0*(1.0/t-1.0/tr) - 18172.91960*np.log(t/tr) \
        + 0.5*0.002676*(t*t-tr*tr)

    ss = 205.15 + 23.10248*np.log(t/tr)  \
        - 2.0*804.8876*(1.0/np.sqrt(t)-1.0/np.sqrt(tr)) \
        - 0.5*1762835.0*(1.0/(t*t)-1.0/(tr*tr)) \
        + 18172.91960*(1.0/t-1.0/tr) + 0.002676*(t-tr)

    mu0_O2= hs - t*ss
    
    logfO2 = np.log10(np.exp((muO2-mu0_O2)/(8.314*t)))

    return logfO2