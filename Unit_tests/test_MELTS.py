import unittest
import pandas as pd
import petthermotools as M

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
                                          CO2_Liq=0.2)['CO2_Liq'][0], 0.2,
decimalPlace, "CO2 doesn't match test value")


if __name__ == '__main__':
    unittest.main()