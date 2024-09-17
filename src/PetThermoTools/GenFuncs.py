import numpy as np
import pandas as pd
# from PetThermoTools.Barom import *
# from PetThermoTools.Liq import *
# from PetThermoTools.Crystallise import *
# from PetThermoTools.MELTS import *
# try:
#     from PetThermoTools.Holland import *
# except:
#     pass

Names = {'liquid1': '_Liq',
        'liquid2': '_Liq2',
        'olivine1': '_Ol',
        'olivine2': '_Ol2',
        'clinopyroxene1': '_Cpx',
        'clinopyroxene2': '_Cpx2',
        'plagioclase1': '_Plag',
        'plagioclase2': '_Plag2',
        'spinel1': '_Sp',
        'spinel2': '_Sp2',
        'k-feldspar1': '_Kspar',
        'k-feldspar2': '_Kspar2',
        'garnet1': '_Grt',
        'garnet2': '_Grt2',
        'rhm-oxide1': '_Rhm',
        'rhm-oxide2': '_Rhm2',
        'quartz1': '_Qtz',
        'quartz2': '_Qtz2',
        'orthopyroxene1': '_Opx',
        'orthopyroxene2': '_Opx2',
        'apatite1': '_Apa',
        'apatite2': '_Apa2'}

Names_MM = {'liq': '_Liq',
            'ol': '_Ol',
            'cpx': '_Cpx',
            'opx': '_Opx',
            'g': '_Grt',
            'pl4t': '_Plag',
            'spn': '_Sp'}

def comp_fix(Model = None, comp = None, Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None):
    '''
    Ensure that the input variables contain the correct column headers for the following variables.

    Parameters:
    ----------
    Model: string
        "MELTSvx.x.x" or "Holland" determines which function list is followed.

    comp: dict or DataFrame
        inputed composition for calculations

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, and type(Fe3Fet_Liq) == np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their Fe3Fet_Liq value. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable or an oxygen fugacity buffer must be chosen.

    H2O_Liq: float or np.ndarray
        H2O content of the initial melt phase. If type(comp) == dict, and type(H2O_Liq) = np.ndarray a new DataFrame will be constructed with bulk compositions varying only in their H2O_Liq value. If comp is a pd.DataFrame, a single H2O_Liq value may be passes (float) and will be used as the initial melt H2O content for all starting compositions. Alternatively, if an array of H2O_Liq values is passed, equal to the number of compositions specified in comp, a different initial melt H2O value will be passed for each sample. If None, H2O_Liq must be specified in the comp variable.

    Returns:
    ---------
    comp: dict or DataFrame
        new composition file with correct headers.
    '''
    if Model is None:
        Model = "MELTSv1.0.2"

    if "MELTS" in Model:
        # check all required columns are present with appropriate suffix
        Columns_bad = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeOt', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'Fe3Fet']
        Columns_ideal = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq']

        Comp_start = comp.copy()
        if "FeO_Liq" in list(Comp_start.keys()) and "Fe2O3_Liq" in list(Comp_start.keys()):
            if "FeOt_Liq" not in list(Comp_start.keys()):
                comp['FeOt_Liq'] = comp['FeO_Liq'] + 71.844/(159.69/2)*comp['Fe2O3_Liq']
            if"Fe3Fet_Liq" not in list(Comp_start.keys()):
                comp['Fe3Fet_Liq'] = (1 - comp['FeO_Liq']/(comp['FeO_Liq'] + 71.844/(159.69/2)*comp['Fe2O3_Liq']))
            Comp_start = comp.copy()

        if "FeO" in list(Comp_start.keys()) and "Fe2O3" in list(Comp_start.keys()):
            if "FeOt" not in list(Comp_start.keys()):
                comp['FeOt'] = comp['FeO'] + 71.844/(159.69/2)*comp['Fe2O3']
            if"Fe3Fet" not in list(Comp_start.keys()):
                comp['Fe3Fet'] = 1 - comp['FeO']/(comp['FeO'] + 71.844/(159.69/2)*comp['Fe2O3'])
            Comp_start = comp.copy()

        if type(comp) == pd.core.frame.DataFrame:
            for el in Comp_start:
                if el in Columns_bad:
                    comp = comp.rename(columns = {el:el + '_Liq'})

            for el in Columns_ideal:
                if el not in list(comp.keys()):
                    comp[el] = np.zeros(len(comp.iloc[:,0]))

        elif type(comp) == dict:
            for el in Comp_start:
                if el in Columns_bad:
                    comp[el + '_Liq'] = comp[el]
                    del comp[el]

            for el in Columns_ideal:
                if el not in list(comp.keys()):
                    comp[el] = 0.0
    else:
        # check all required columns are present with appropriate suffix
        Columns_bad = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O', 'H2O', 'Cr2O3' 'Fe3Fet']
        Columns_ideal = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'Cr2O3_Liq', 'H2O_Liq', 'Fe3Fet_Liq']
        Comp_start = comp.copy()
        if type(comp) == pd.core.frame.DataFrame:
            for el in Comp_start:
                if el in Columns_bad:
                    comp = comp.rename(columns = {el:el + '_Liq'})

            for el in Columns_ideal:
                if el not in list(comp.keys()):
                    comp[el] = np.zeros(len(comp.iloc[:,0]))

        elif type(comp) == dict:
            for el in Comp_start:
                if el in Columns_bad:
                    comp[el + '_Liq'] = comp[el]
                    del comp[el]

            for el in Columns_ideal:
                if el not in list(comp.keys()):
                    comp[el] = 0.0

    # set the liquid Fe redox state if specified separate to the bulk composition
    if Fe3Fet_Liq is not None:
        if type(comp) == dict:
            if type(Fe3Fet_Liq) != np.ndarray:
                comp['Fe3Fet_Liq'] = Fe3Fet_Liq
            else:
                Comp = pd.DataFrame.from_dict([comp]*len(Fe3Fet_Liq))
                Comp['Fe3Fet_Liq'] = Fe3Fet_Liq
                comp = Comp.copy()
        else:
            comp['Fe3Fet_Liq'] = np.zeros(len(comp.iloc[:,0])) + Fe3Fet_Liq


    if H2O_Liq is not None:
        if type(comp) == dict:
            if type(H2O_Liq) != np.ndarray:
                comp['H2O_Liq'] = H2O_Liq
            else:
                Comp = pd.DataFrame.from_dict([comp]*len(H2O_Liq))
                Comp['H2O_Liq'] = H2O_Liq
                comp = Comp.copy()
        else:
            comp['H2O_Liq'] = np.zeros(len(comp.iloc[:,0])) + H2O_Liq

    if CO2_Liq is not None:
        if type(comp) == dict:
            if type(CO2_Liq) != np.ndarray:
                comp['CO2_Liq'] = CO2_Liq
            else:
                Comp = pd.DataFrame.from_dict([comp]*len(CO2_Liq))
                Comp['CO2_Liq'] = CO2_Liq
                comp = Comp.copy()
        else:
            comp['CO2_Liq'] = np.zeros(len(comp.iloc[:,0])) + CO2_Liq

    return comp

def stich(Res, multi = None, Model = None, Frac_fluid = None, Frac_solid = None):
    '''
    Takes the outputs from the multiple crystallisation/decompression calculations and stiches them together into a single dataframe. Additionally, it adds the relevant suffix to the composition and properties of each mineral (e.g., SiO2 -> SiO2_Liq for the liquid phase).

    Parameters:
    ----------
    Res: dict
        Final results from the multiple crystallisation/decompression calculations.

    multi: True/False
        If True, Results is composed of multiple dictionaries, each representing a single crystallisation/decompression calculation. Default is False.

    Model: string
        "MELTSvx.x.x" or "Holland" determines which function list is followed.

    Returns:
    ----------
    Results: dict
        A copy of the input dict with a new DataFrame titled 'All' included.
    '''
    Results = Res.copy()
    if "MELTS" in Model:
        Order = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'FeO', 'FeOt', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'Fe3Fet']
        if multi is None:
            Results = stich_work(Results = Results, Order = Order, Model = "MELTS", Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
        else:
            for Ind in Res:
                Result = Res[Ind].copy()
                Result = stich_work(Results = Result, Order = Order, Model = "MELTS", Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
                Results[Ind] = Result.copy()
    else:
        Order = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O', 'H2O', 'Fe3Fet']
        if multi is None:
            Results = stich_work(Results = Results, Order = Order, Model = "Holland", Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
        else:
            for Ind in Res:
                Result = Res[Ind].copy()
                Result = stich_work(Results = Result, Order = Order, Model = "Holland", Frac_fluid = Frac_fluid, Frac_solid = Frac_solid)
                Results[Ind] = Result.copy()

    return Results

def stich_work(Results = None, Order = None, Model = "MELTS", Frac_fluid = None, Frac_solid = None):
    '''
    Does the work required by Stich.
    '''
    Res = Results.copy()
    if "MELTS" in Model:    
        Results['Conditions'] = Results['Conditions'].rename(columns = {'temperature':'T_C'})
        Results['Conditions'] = Results['Conditions'].rename(columns = {'pressure':'P_bar'})
    else:
        Results['Conditions'] = Results['Conditions'].rename(columns = {'P_kbar': 'P_bar'})
        Results['Conditions']['P_bar'] = Results['Conditions']['P_bar']*1000

    SN = []
    for R in Results:
        if '_prop' not in R and R != 'Conditions' and R != "sys":
            SN += [R]

    for R in SN:
        # if "_prop" not in R and R != "Conditions" and R!= "sys":
        if "MELTS" in Model:
            Results[R].loc[:,'FeOt'] = Results[R].loc[:,'FeO'] + 71.844/(158.69/2)*Results[R].loc[:,'Fe2O3']
            Results[R].loc[:,'Fe3Fet'] = (71.844/(159.69/2)*Results[R].loc[:,'Fe2O3'])/Results[R].loc[:,'FeOt']
            try:
                Results[R][Results[R + '_prop']['mass'] == 0.0] = np.nan
            except:
                Results[R][Results[R + '_prop']['Mass'] == 0.0] = np.nan
            Results[R] = Results[R][Order]
            if R == "fluid1":
                El = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'FeOt', 'Fe3Fet',
                      'MnO','MgO', 'CaO', 'Na2O', 'K2O', 'P2O5']
                for e in El:
                    Results[R] = Results[R].drop(columns = e)
                Results[R].loc[:,'X_H2O_mol'] = (Results[R].loc[:, 'H2O']/18)/(Results[R].loc[:, 'H2O']/18 + Results[R].loc[:, 'CO2']/44)
                Results[R].loc[:,'X_CO2_mol'] = 1 - Results[R].loc[:, 'X_H2O_mol']
        else:
            Results[R] = Results[R].rename(columns = {'FeO': 'FeOt'})
            Tot = Results[R].sum(axis = 1)
            for el in Results[R]:
                Results[R][el] = 100*Results[R][el]/Tot
            Results[R]['Fe3Fet'] = Results[R]['O']/(((159.9/2)/71.844)*Results[R]['FeOt'] - Results[R]['FeOt'])
            try:
                Results[R][Results[R + '_prop']['mass'] == 0.0] = np.nan
            except:
                Results[R][Results[R + '_prop']['Mass'] == 0.0] = np.nan
                
            Results[R] = Results[R][Order]

    if "MELTS" in Model:
        Remove = np.where(Results['Conditions']['h'] == 0.0)[0]
    else:
        Remove = np.where(Results['Conditions']['T_C'] == 0.0)[0]

    for R in Results:
        if len(Remove) > 0:
            Results[R] = Results[R].drop(labels = Remove)

    # if "MELTS" in Model:
    #     Remove_2 = np.where(Results['Conditions']['h'] == np.nan)[0]
    #     for R in Results:
    #         if len(Remove_2) > 0:
    #             Results[R] = Results[R].drop(labels = Remove_2)

    Results_Mass = pd.DataFrame(data = np.zeros((len(Results['Conditions']['T_C']), len(SN))), columns = SN)
    Results_Volume = Results_Mass.copy()
    Results_rho = Results_Mass.copy()
    for n in SN:
        try:
            Results_Mass[n] = Results[n + '_prop']['mass']
        except:
            Results_Mass[n] = Results[n + '_prop']['Mass']

        if "MELTS" in Model:
            Results_Volume[n] = Results[n + '_prop']['v']
            Results_rho[n] = Results[n + '_prop']['rho']

    if Frac_solid is True or Frac_fluid is True:
        if Frac_solid is None:
            Results_Mass['fluid1_cumsum'] = Results_Mass['fluid1'].cumsum()
        elif Frac_fluid is None:
            for n in SN:
                if n != 'liquid1' and n!= 'fluid1':
                    Results_Mass[n + '_cumsum'] = Results_Mass[n].cumsum()
        else:
            for n in SN:
                if n != 'liquid1':
                    Results_Mass[n + '_cumsum'] = Results_Mass[n].cumsum()

    Results_All = Results['Conditions'].copy()
    for R in Results:
        if R != "Conditions" and R != "sys":
            if "MELTS" in Model:
                if any(n in R for n in Names):
                    for n in Names:
                        if n in R:
                            Results[R] = Results[R].add_suffix(Names[n])                
                else:
                    if '_prop' in R:
                        Results[R] = Results[R].add_suffix('_' + R[:-5])
                    else:
                        Results[R] = Results[R].add_suffix('_' + R)
            else:
                if any(n in R for n in Names_MM):
                    for n in Names_MM:
                        if n in R:
                            Results[R] = Results[R].add_suffix(Names_MM[n])
                else:       
                    if '_prop' in R:
                        Results[R] = Results[R].add_suffix('_' + R[:-5])
                    else:
                        Results[R] = Results[R].add_suffix('_' + R)

            Results_All = pd.concat([Results_All, Results[R]], axis = 1)

    Results['All'] = Results_All
    Results['Mass'] = Results_Mass
    if "MELTS" in Model:
        Results['Volume'] = Results_Volume
        Results['rho'] = Results_rho

    if Results['Mass'].sum(axis = 1).iloc[-1] == 0.0:
        for R in Results:
            Results[R].drop(Results[R].tail(1).index,inplace=True)

    return Results