import numpy as np 
import pandas as pd 
import re
from petthermotools.GenFuncs import *

# Ionic Radii (in Angstroms, Å) and Valencies for supplied elements
element_data = {
    "Na": {"valency": 1, "ionic_radius": 1.18},
    "P": {"valency": 5, "ionic_radius": 0.38},
    "K": {"valency": 1, "ionic_radius": 1.51},
    "Sc": {"valency": 3, "ionic_radius": 0.870},
    "Ti": {"valency": 4, "ionic_radius": 0.74},
    "V": {"valency": 5, "ionic_radius": 0.54},
    "Cr": {"valency": 3, "ionic_radius": 0.615},
    "Mn": {"valency": 2, "ionic_radius": 0.96},
    "Co": {"valency": 2, "ionic_radius": 0.90},
    "Ni": {"valency": 2, "ionic_radius": 0.690},
    "Ga": {"valency": 3, "ionic_radius": 0.620},
    "Rb": {"valency": 1, "ionic_radius": 1.61},
    "Sr": {"valency": 2, "ionic_radius": 1.26},
    "Y": {"valency": 3, "ionic_radius": 1.019},
    "Nb": {"valency": 5, "ionic_radius": 0.64},
    "Cd": {"valency": 2, "ionic_radius": 1.10},
    "In": {"valency": 3, "ionic_radius": 0.92},
    "Cs": {"valency": 1, "ionic_radius": 1.74},
    "Ba": {"valency": 2, "ionic_radius": 1.42},
    "La": {"valency": 3, "ionic_radius": 1.16},
    "Ce": {"valency": 3, "ionic_radius": 1.143},
    "Pr": {"valency": 3, "ionic_radius": 1.126},
    "Nd": {"valency": 3, "ionic_radius": 1.109},
    "Sm": {"valency": 3, "ionic_radius": 1.079},
    "Eu": {"valency": 3, "ionic_radius": 1.066},
    "Gd": {"valency": 3, "ionic_radius": 1.053},
    "Tb": {"valency": 3, "ionic_radius": 1.040},
    "Dy": {"valency": 3, "ionic_radius": 1.027},
    "Ho": {"valency": 3, "ionic_radius": 1.015},
    "Er": {"valency": 3, "ionic_radius": 1.004},
    "Tm": {"valency": 3, "ionic_radius": 0.994},
    "Yb": {"valency": 3, "ionic_radius": 0.985},
    "Lu": {"valency": 3, "ionic_radius": 0.977},
    "Ta": {"valency": 5, "ionic_radius": 0.64},
    "Ra": {"valency": 2, "ionic_radius": 1.48},
    "Ac": {"valency": 3, "ionic_radius": 1.25},
    "Th": {"valency": 4, "ionic_radius": 1.035},
    "Pa": {"valency": 5, "ionic_radius": 0.78},
}

# Simple lookup dictionary for ONLY ionic radii
ionic_radii = {elem: data["ionic_radius"] for elem, data in element_data.items()}

radii = pd.DataFrame.from_dict(ionic_radii, orient = 'index').T

trivalent_radii = {
    elem: data["ionic_radius"]
    for elem, data in element_data.items()
    if data["valency"] == 3
}

trivalent_radii_df = pd.DataFrame.from_dict(trivalent_radii, orient='index').T
trivalent_radii_names = trivalent_radii_df.columns

"""Composition handling for the trace engine: one clean step, one normalise step.

    wt  = oxides(df, phase="px", fe_mode="total")   # clean wt% table
    cat = normalise(wt, "oxygen", 6)                # cations per 6 O

Replaces ensure_oxide_columns, any_cation_fractions, _oxides, px_oxides, _mol,
strip and num. The old names are kept below as one-line wrappers.
"""
import re
import pandas as pd

from petthermotools.GenFuncs import OXIDE_MW, OXIDE_N, CATION_N

# One table of constants instead of three dicts indexed in parallel
OXIDES = pd.DataFrame({"mw": OXIDE_MW, "n_cat": CATION_N, "n_O": OXIDE_N}).dropna()
VOLATILES = ["H2O", "CO2"]
_FE2O3_TO_FEO = 2 * OXIDES.at["FeO", "mw"] / OXIDES.at["Fe2O3", "mw"]

# Column suffixes recognised for each phase (case-insensitive; anything after is dropped,
# so _Cpx, _cpx1, _clinopyroxene1 all work). phase=None strips any "_word" suffix.
SUFFIXES = {
    "px": "cpx|opx|clinopyroxene|orthopyroxene",
    "ol": "ol|olivine",
    "liq": "liq|liquid|melt",
    "grt": "Grt|grt|g|garnet"
}


def oxides(df, phase=None, fe_mode="as_given", volatiles=True):
    """Clean wt% table with one column per oxide in OXIDES and nothing else.

    - suffixes stripped; non-oxide columns (Sample_ID, P_kbar, other phases) ignored
    - missing columns and blank/NaN/text cells -> 0
    - iron, row by row:
        "as_given": FeO/Fe2O3 where given, otherwise FeOt (split by Fe3Fet if present)
        "total":    all Fe as FeO (FeOt if > 0, else FeO + 0.8998*Fe2O3)
    - volatiles=False drops H2O and CO2
    """
    if fe_mode not in ("as_given", "total"):
        raise ValueError("fe_mode must be 'as_given' or 'total'")

    pat = re.compile(rf"_({SUFFIXES[phase] if phase else '[A-Za-z]'}).*$", re.IGNORECASE)
    wanted = [*OXIDES.index, "FeOt", "Fe3Fet"]
    cols = {}
    for c in df.columns:
        name = pat.sub("", c) if isinstance(c, str) else c
        if name in wanted:
            if name in cols:
                raise ValueError(f"'{name}' appears twice ({cols[name]!r}, {c!r}); pass phase= to pick one")
            cols[name] = c

    wt = pd.DataFrame({n: pd.to_numeric(df[c], errors="coerce") for n, c in cols.items()},
                      index=df.index).reindex(columns=wanted).fillna(0.0)

    FeOt, Fe3Fet = wt.pop("FeOt"), wt.pop("Fe3Fet")
    if fe_mode == "total":
        wt["FeO"] = FeOt.where(FeOt > 0, wt["FeO"] + _FE2O3_TO_FEO * wt["Fe2O3"])
        wt["Fe2O3"] = 0.0
    else:
        no_fe = (wt["FeO"] + wt["Fe2O3"]) <= 0
        wt.loc[no_fe, "FeO"] = (FeOt * (1 - Fe3Fet))[no_fe]
        wt.loc[no_fe, "Fe2O3"] = (FeOt * Fe3Fet / _FE2O3_TO_FEO)[no_fe]

    return wt if volatiles else wt.drop(columns=VOLATILES, errors="ignore")


def normalise(wt, basis, total=1.0):
    """Normalise a clean wt% table (from `oxides`).

    basis="oxide":  oxide mole fractions summing to `total`
    basis="cation": cation fractions summing to `total`
    basis="oxygen": cations per `total` oxygens (apfu); columns keep oxide names
    """
    k = OXIDES.loc[wt.columns]
    mol = wt / k["mw"]
    if basis == "oxide":
        return mol.div(mol.sum(axis=1), axis=0) * total
    cat = mol * k["n_cat"]
    if basis == "cation":
        return cat.div(cat.sum(axis=1), axis=0) * total
    if basis == "oxygen":
        return cat.div((mol * k["n_O"]).sum(axis=1), axis=0) * total
    raise ValueError("basis must be 'oxide', 'cation' or 'oxygen'")


# ---- old names, now one line each -------------------------------------------------
def px_apfu(pyx, fe_mode="as_given"):
    return normalise(oxides(pyx, "px", fe_mode, volatiles=False), "oxygen", 6)


def px_cation_fraction(pyx, fe_mode="as_given"):
    return normalise(oxides(pyx, "px", fe_mode, volatiles=False), "cation")


def liq_mole_fractions(liquid, fe_mode="as_given"):
    return normalise(oxides(liquid, "liq", fe_mode), "oxide")


def any_cation_fractions(df, return_total_cations=False, phase=None, fe_mode="as_given"):
    wt = oxides(df, phase, fe_mode)
    frac = normalise(wt, "cation")
    if return_total_cations:
        return frac, (wt / OXIDES.loc[wt.columns, "mw"] * OXIDES.loc[wt.columns, "n_cat"]).sum(axis=1)
    return frac