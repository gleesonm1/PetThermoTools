import numpy as np
import pyMelt as m

# Saved compositions for use in pyMELTScalc
Compositions = {'KLB-1': {'SiO2_Liq': 44.48,
                        'TiO2_Liq': 0.16,
                        'Al2O3_Liq': 3.59,
                        'FeOt_Liq': 8.10,
                        'MgO_Liq': 39.22,
                        'CaO_Liq': 3.44,
                        'Na2O_Liq': 0.30,
                        'K2O_Liq': 0.02,
                        'Cr2O3_Liq': 0.42,
                        'H2O_Liq': 0.00,
                        'Fe3Fet_Liq': 0.040},
                'KG1': {'SiO2_Liq': 47.0,
                        'TiO2_Liq': 0.78,
                        'Al2O3_Liq': 9.75,
                        'FeOt_Liq': 9.77,
                        'MgO_Liq': 23.6,
                        'CaO_Liq': 7.35,
                        'Na2O_Liq': 1.52,
                        'K2O_Liq': 0.12,
                        'Cr2O3_Liq': 0.20,
                        'H2O_Liq': 0.0,
                        'Fe3Fet_Liq': 0.155},
                'G2': {'SiO2_Liq': 50.05,
                        'TiO2_Liq': 1.97,
                        'Al2O3_Liq': 15.76,
                        'FeOt_Liq': 9.35,
                        'MgO_Liq': 7.90,
                        'CaO_Liq': 11.74,
                        'Na2O_Liq': 3.04,
                        'K2O_Liq': 0.03,
                        'Cr2O3_Liq': 0.0,
                        'H2O_Liq': 0.0,
                        'Fe3Fet_Liq': 0.18}}

Lithologies = {'KLB-1': m.lithologies.matthews.klb1(),
               'KG1': m.lithologies.matthews.kg1(),
               'G2': m.lithologies.matthews.eclogite(),
               'hz': m.lithologies.shorttle.harzburgite()}