import numpy as np
import pandas as pd
# from pyMELTScalc.Barom import *
# from pyMELTScalc.Liq import *
# from pyMELTScalc.Crystallise import *
# from pyMELTScalc.MELTS import *
# try:
#     from pyMELTScalc.Holland import *
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

def comp_fix(Model = None, comp = None, Fe3Fet_Liq = None):
    '''
    Ensure that the input variables contain the correct column headers for the following variables.

    Parameters:
    ----------
    Model: string
        "MELTSvx.x.x" or "Holland" determines which function list is followed.

    comp: dict or DataFrame
        inputed composition for calculations

    Fe3Fet_Liq: float or np.ndarray
        Fe 3+/total ratio. If type(comp) == dict, Fe3Fet_Liq must be a float and will set the Fe redox state in the initial composition. If comp is a pd.DataFrame, a single Fe3Fet_Liq value may be passed (float) and will be used as the Fe redox state for all starting compostions, or an array of Fe3Fet_Liq values, equal to the number of compositions specified in comp can specify a different Fe redox state for each sample. If None, the Fe redox state must be specified in the comp variable.

    Returns:
    ---------
    comp: dict or DataFrame
        new composition file with correct headers.
    '''
    # set the liquid Fe redox state if specified separate to the bulk composition
    if Fe3Fet_Liq is not None:
        if type(comp) == dict:
            comp['Fe3Fet_Liq'] = Fe3Fet_Liq
        else:
            comp['Fe3Fet_Liq'] = np.zeros(len(comp.iloc[:,0])) + Fe3Fet_Liq

    if "MELTS" in Model:
        # check all required columns are present with appropriate suffix
        Columns_bad = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'Fe3Fet']
        Columns_ideal = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq']
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

    return comp

def stich(Res, multi = None, Model = None):
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
    Order = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'FeO', 'FeOt', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'Fe3Fet']
    if multi is None:

        Results['Conditions'] = Results['Conditions'].rename(columns = {'temperature':'T_C'})
        Results['Conditions'] = Results['Conditions'].rename(columns = {'pressure':'P_bar'})

        for R in Results:
            if "_prop" not in R and R != "Conditions":
                Results[R]['FeOt'] = Results[R]['FeO'] + 71.844/(158.69/2)*Results[R]['Fe2O3']
                Results[R]['Fe3Fet'] = (71.844/(159.69/2)*Results[R]['Fe2O3'])/Results[R]['FeOt']
                Results[R][Results[R + '_prop']['mass'] == 0.0] = np.nan
                Results[R] = Results[R][Order]
            if len(np.where(Res['Conditions']['temperature'] == 0.0)[0]) > 0:
                Results[R] = Results[R].drop(labels = np.where(Res['Conditions']['temperature'] == 0.0)[0])

        Results_All = Results['Conditions']

        for R in Results:
            if R != "Conditions":
                if any(n in R for n in Names):
                    for n in Names:
                        if n in R:
                            Results[R] = Results[R].add_suffix(Names[n])
                else:
                    Results[R] = Results[R].add_suffix('_' + R)


                Results_All = pd.concat([Results_All, Results[R]], axis = 1)

        Results['All'] = Results_All


    else:
        for Ind in Res:
            Result = Res[Ind].copy()
            Result['Conditions'] = Result['Conditions'].rename(columns = {'temperature':'T_C'})
            Result['Conditions'] = Result['Conditions'].rename(columns = {'pressure':'P_bar'})

            for R in Result:
                if "_prop" not in R and R != "Conditions":
                    Result[R]['FeOt'] = Result[R]['FeO'] + 71.844/(158.69/2)*Result[R]['Fe2O3']
                    Result[R]['Fe3Fet'] = (71.844/(159.69/2)*Result[R]['Fe2O3'])/Result[R]['FeOt']
                    Result[R][Result[R + '_prop']['mass'] == 0.0] = np.nan
                    Result[R] = Result[R][Order]
                if len(np.where(Res[Ind]['Conditions']['temperature'] == 0.0)[0]) > 0:
                    Result[R] = Result[R].drop(labels = np.where(Res[Ind]['Conditions']['temperature'] == 0.0)[0])

            Results_All = Result['Conditions']

            for R in Result:
                if R != "Conditions":
                    if any(n in R for n in Names):
                        for n in Names:
                            if n in R:
                                Result[R] = Result[R].add_suffix(Names[n])
                    else:
                        Result[R] = Result[R].add_suffix('_' + R)


                    Results_All = pd.concat([Results_All, Result[R]], axis = 1)

            Result['All'] = Results_All

            Results[Ind] = Result.copy()

    return Results