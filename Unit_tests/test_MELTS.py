import unittest
import pandas as pd
import pyMELTScalc as M

## Setting a liquid composition

Liq_test=pd.DataFrame(data={'SiO2_Liq': 52.72,
                            'TiO2_Liq': 2.08,
                            'Al2O3_Liq': 13.26,
                            'FeOt_Liq': 9.21,
                            'MgO_Liq': 9.32,
                            'CaO_Liq': 10.26,
                            'Na2O_Liq': 2.24,
                            'K2O_Liq': 0.43,
                            'P2O5_Liq': 0.22,
                            'Fe3Fet_Liq': 0.095,
                            'H2O_Liq': 0}, index=[0])
decimalPlace=2

class test_compfix(unittest.TestCase):
    def test_compfix_single_df(self):
        self.assertAlmostEqual(M.comp_fix(Model="pMELTS",
                                          comp=Liq_test,
                                          CO2_Liq=0.2)['CO2_Liq'], 0.2,
decimalPlace, "CO2 doesn't match test value")

#     def test_Simple_mantle_S(self):
#         self.assertAlmostEqual(ss.Lee_Wieser_sulfide_melting(Modes=Modes, KDs=KDs_Cu,
#                         N=100, S_Mantle=200,
#                         S_Sulf=S_Sulf, S_Melt_SCSS_2_ppm=1000,
#                         elem_Per=30, Prop_S6=0)['S_Residue'][3], 175.2577319587629,
# decimalPlace, "Residual S doesnt match test value")



# class test_SCSS_SCAS_Total(unittest.TestCase):
#     def test_SCSS_SCAS_Jugo(self):
#         self.assertAlmostEqual(ss.calculate_S_Total_SCSS_SCAS(deltaQFM=0.3,
#             SCSS=1000, SCAS=5000, model='Jugo')['Total_S_ppm'][0], 1031.6227766016839,
# decimalPlace, "ST from Jugo doesnt match test value")

#     def test_SCSS_SCAS_Nash(self):
#         self.assertAlmostEqual(ss.calculate_S_Total_SCSS_SCAS(Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'], T_K=Liq_test['T_K'],
#             SCSS=1000, SCAS=5000, model='Nash')['Total_S_ppm'][0], 1006.9936479729372,
# decimalPlace, "ST from Nash doesnt match test value")

#     def test_SCSS_SCAS_Nash_higher(self):
#         self.assertAlmostEqual(ss.calculate_S_Total_SCSS_SCAS(Fe3Fet_Liq=Liq_test['Fe3Fet_Liq']*5, T_K=Liq_test['T_K'],
#             SCSS=1000, SCAS=5000, model='Nash')['Total_S_ppm'][0], 5000.02347426457,
# decimalPlace, "ST from Nash doesnt match test value")

# class test_S6_corrections(unittest.TestCase):
#     def test_s6_Jugo(self):
#         self.assertAlmostEqual(ss.calculate_S6St_Jugo2010_eq10(deltaQFM=0.3), 0.030653430031715508,
# decimalPlace, "S6 from Jugo doesnt match test value")

#     def test_s6_Nash(self):
#         self.assertAlmostEqual(ss.calculate_S6St_Nash2019(T_K=Liq_test['T_K'], Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])[0],0.006945076552384788,
# decimalPlace, "S6 from Nash doesnt match test value")

#     def test_O2022_O2023(self):
#         self.assertAlmostEqual(ss.calculate_OM2022_S6St(df=Liq_test, T_K=Liq_test['T_K'], logfo2=None,
#                     Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])['S6St_Liq'][0],0.008973,
# decimalPlace, "S6 from ONeill doesnt match test value")

# class test_SCSS_calc_sulfide(unittest.TestCase):
#     def test_Smythe2017_calcsulf_Smythe(self):
#         self.assertAlmostEqual(ss.calculate_S2017_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Ni_Liq=100, Cu_Liq=150,
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Fe_FeNiCu_Sulf='Calc_Smythe')['SCSS2_ppm_ideal_Smythe2017'][0], 968.5364514652154,
# decimalPlace, "SCSS calculated from Smythe2017 using Calc Sulf Smythe not equal to test value")

#     def test_Smythe2017_calcsulf_ONeill(self):
#         self.assertAlmostEqual(ss.calculate_S2017_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Ni_Liq=100, Cu_Liq=150,
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Fe_FeNiCu_Sulf='Calc_ONeill')['SCSS2_ppm_ideal_Smythe2017'][0], 930.3940622962762,
# decimalPlace, "SCSS calculated from Smythe2017 using Calc Sulf Smythe not equal to test value")

#     def test_O2021_calcsulf_ONeill(self):
#         self.assertAlmostEqual(ss.calculate_O2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Ni_Liq=100, Cu_Liq=150,
# Fe_FeNiCu_Sulf='Calc_ONeill')['SCSS2_ppm'][0], 955.4837286314,
# decimalPlace, "SCSS calculated from O2021 with Oneill Sulf not equal to test value")

#     def test_O2021_calcsulf_Smythe(self):
#         self.assertAlmostEqual(ss.calculate_O2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Ni_Liq=100, Cu_Liq=150,
# Fe_FeNiCu_Sulf='Calc_Smythe')['SCSS2_ppm'][0], 994.6546925260972,
# decimalPlace, "SCSS calculated from O2021 with Smythe Sulf not equal to test value")


#     def test_Bl2021_calcsulf_Smythe(self):
#         self.assertAlmostEqual(ss.calculate_B2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Ni_Liq=100, Cu_Liq=150,
# Fe_FeNiCu_Sulf='Calc_Smythe')['SCSS2_ppm_eq11'][0], 1074.7469538219757,
# decimalPlace, "SCSS calculated from B2021 not equal to test value")



#     def test_L2021_calcsulf_Smythe(self):
#         self.assertAlmostEqual(ss.calculate_Liu2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Ni_Liq=100, Cu_Liq=150,
# Fe_FeNiCu_Sulf='Calc_Smythe')['SCSS2_ppm'][0], 1034.4417069163972,
# decimalPlace, "SCSS calculated from L2021 not equal to test value")





# class test_SCSS_fixed_sulfide(unittest.TestCase):
#     def test_Smythe2017(self):
#         self.assertAlmostEqual(ss.calculate_S2017_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Fe_FeNiCu_Sulf=0.65)['SCSS2_ppm_ideal_Smythe2017'][0], 971.1628818423652,
# decimalPlace, "SCSS calculated from Smythe2017 not equal to test value")

#     def test_O2021(self):
#         self.assertAlmostEqual(ss.calculate_O2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Fe_FeNiCu_Sulf=0.65)['SCSS2_ppm'][0], 997.3519490880699,
# decimalPlace, "SCSS calculated from O2021 not equal to test value")


#     def test_Bl2021(self):
#         self.assertAlmostEqual(ss.calculate_B2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'],
# Fe_FeNiCu_Sulf=0.65)['SCSS2_ppm_eq11'][0], 1077.661400710366,
# decimalPlace, "SCSS calculated from B2021 not equal to test value")

#     def test_F2015(self):
#         self.assertAlmostEqual(ss.calculate_F2015_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])['SCSS2_ppm'][0], 1367.444449878481,
# decimalPlace, "SCSS calculated from F2015 not equal to test value")

#     def test_L2021(self):
#         self.assertAlmostEqual(ss.calculate_Liu2021_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'], Fe_FeNiCu_Sulf=0.65)['SCSS2_ppm'][0], 1037.2468559826248,
# decimalPlace, "SCSS calculated from L2021 not equal to test value")

#     def test_OM2022(self):
#         self.assertAlmostEqual(ss.calculate_OM2022_SCSS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'], Fe_FeNiCu_Sulf=0.65)['SCSS2_ppm'][0], 997.3519490880699,
# decimalPlace, "SCSS calculated from OM2022 not equal to test value")

# class test_SCAS_fixed_sulfide(unittest.TestCase):
#     def test_CD2019(self):
#         self.assertAlmostEqual(ss.calculate_CD2019_SCAS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])['SCAS6_ppm'][0], 7043.40006526207,
# decimalPlace, "SCAS calculated from CD2019 not equal to test value")

#     def test_ZT2022(self):
#         self.assertAlmostEqual(ss.calculate_ZT2022_SCAS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])['SCAS6_ppm'][0], 7025.572329292286,
# decimalPlace, "SCAS calculated from ZT2022 not equal to test value")

#     def test_MK2015(self):
#         self.assertAlmostEqual(ss.calculate_MK2015_SCAS(df=Liq_test,
# T_K=Liq_test['T_K'], P_kbar=Liq_test['P_kbar'],
# Fe3Fet_Liq=Liq_test['Fe3Fet_Liq'])['SCAS6_ppm'][0], 0.005101492327453715,
# decimalPlace, "SCAS calculated from MK2015 not equal to test value")






if __name__ == '__main__':
     unittest.main()