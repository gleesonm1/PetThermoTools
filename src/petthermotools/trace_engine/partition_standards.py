import numpy as np 
import pandas as pd

phases = ['olivine', 'orthopyroxene', 'clinopyroxene', 'sphene', 'garnet', 'melilite', 'aenigmatite', 'cummingtonite',
 'clinoamphibole', 'orthoamphibole', 'hornblende', 'biotite', 'muscovite', 'alkali-feldspar', 'plagioclase', 'quartz',
 'tridymite', 'cristobalite', 'nepheline', 'kalsilite', 'leucite', 'corundum', 'sillimanite', 'rutile', 'perovskite',
 'spinel', 'rhm-oxide', 'ortho-oxide', 'whitlockite', 'apatite', 'fluid', 'alloy-solid', 'alloy-liquid', 'calcite', 'aragonite',
 'magnesite', 'siderite', 'dolomite', 'spurrite', 'tilleyite', 'diamond', 'graphite']

# elements = [
#     "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
#     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
#     "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
#     "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
#     "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
#     "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
#     "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
#     "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
#     "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
#     "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cc", "Es", "Fm",
#     "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
#     "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
# ]

common_elements = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
]

coeffs_standard = {
    phase: {element: 0 for element in common_elements} for phase in phases
}

coeffs = {"Default": coeffs_standard,
"D_Gibson_2010": {
    'olivine': {
        "Rb": 0.0003, "Ba": 0.000005, "Th": 0.00005, "U": 0.00038, "Nb": 0.0005,
        "Ta": 0.0005, "La": 0.0005, "Ce": 0.0005, "Pb": 0.003, "Pr": 0.0008,
        "Nd": 0.00042, "Sr": 0.00004, "Zr": 0.0033, "Hf": 0.0022, "Sm": 0.0011,
        "Eu": 0.0016, "Ti": 0.015, "Gd": 0.0011, "Tb": 0.0015, "Dy": 0.0027,
        "Ho": 0.0016, "Y": 0.0099, "Er": 0.013, "Yb": 0.020, "Lu": 0.020
    },
    'orthopyroxene': {
        "Rb": 0.0002, "Ba": 0.000006, "Th": 0.002, "U": 0.002, "Nb": 0.004,
        "Ta": 0.004, "La": 0.0031, "Ce": 0.0040, "Pb": 0.009, "Pr": 0.0048,
        "Nd": 0.01200, "Sr": 0.0007, "Zr": 0.013, "Hf": 0.03, "Sm": 0.0200,
        "Eu": 0.0130, "Ti": 0.086, "Gd": 0.0130, "Tb": 0.0190, "Dy": 0.0110,
        "Ho": 0.0065, "Y": 0.052, "Er": 0.045, "Yb": 0.080, "Lu": 0.120
    },
    'clinopyroxene': {
        "Rb": 0.0004, "Ba": 0.0004, "Th": 0.0059, "U": 0.0094, "Nb": 0.015,
        "Ta": 0.015, "La": 0.0490, "Ce": 0.0800, "Pb": 0.012, "Pr": 0.126,
        "Nd": 0.17800, "Sr": 0.091, "Zr": 0.119, "Hf": 0.284, "Sm": 0.2930,
        "Eu": 0.3350, "Ti": 0.350, "Gd": 0.3500, "Tb": 0.4030, "Dy": 0.4000,
        "Ho": 0.4270, "Y": 0.426, "Er": 0.420, "Yb": 0.400, "Lu": 0.376
    },
    'garnet': {
        "Rb": 0.0002, "Ba": 0.00007, "Th": 0.009, "U": 0.028, "Nb": 0.015,
        "Ta": 0.015, "La": 0.0010, "Ce": 0.0050, "Pb": 0.005, "Pr": 0.014,
        "Nd": 0.05200, "Sr": 0.0007, "Zr": 0.270, "Hf": 0.400, "Sm": 0.2500,
        "Eu": 0.4960, "Ti": 0.600, "Gd": 0.84800, "Tb": 1.4770, "Dy": 2.2000,
        "Ho": 3.3150, "Y": 3.100, "Er": 4.400, "Yb": 6.600, "Lu": 7.100
    },
    'spinel': {
        "Rb": 0.0001, "Ba": 0.0001, "Th": 0.0, "U": 0.0, "Nb": 0.0, "Ta": 0.0,
        "La": 0.0100, "Ce": 0.0100, "Pb": 0.0, "Pr": 0.01, "Nd": 0.0100,
        "Sr": 0.0, "Zr": 0.0, "Hf": 0.0, "Sm": 0.0100, "Eu": 0.0100, "Ti": 0.15,
        "Gd": 0.0100, "Tb": 0.0100, "Dy": 0.0100, "Ho": 0.0100, "Y": 0.01,
        "Er": 0.0100, "Yb": 0.0100, "Lu": 0.0100
    },
    'plagioclase': {
        "Rb": 0.03, "Ba": 0.33, "Th": 0.05, "U": 0.11, "Nb": 0.01, "Ta": 0.0,
        "La": 0.2700, "Ce": 0.200, "Pb": 0.36, "Pr": 0.17, "Nd": 0.1400,
        "Sr": 2.0, "Zr": 0.01, "Hf": 0.01, "Sm": 0.1100, "Eu": 0.7300, "Ti": 0.04,
        "Gd": 0.0660, "Tb": 0.0600, "Dy": 0.0550, "Ho": 0.0480, "Y": 0.03,
        "Er": 0.0100, "Yb": 0.031, "Lu": 0.0250
    }
}}