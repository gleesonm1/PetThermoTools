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

PX_SUFFIX = re.compile(r"_(cpx|opx|clinopyroxene|orthopyroxene).*$", flags=re.IGNORECASE)
LIQ_SUFFIX = re.compile(r"_(liq|liquid|melt).*$", flags=re.IGNORECASE)

def _oxides(df, suffix_pattern, fe_mode):
    if fe_mode not in ("as_given", "total"):
        raise ValueError("fe_mode must be 'as_given' or 'total'")
    px = strip(df, suffix_pattern)
 
    ox = [k for k in OXIDE_MW
          if k in CATION_N and k in OXIDE_N and k not in ("FeOt", "Fe3Fet")]
    out = pd.DataFrame({k: num(px, k) for k in ox}, index=px.index)
 
    FeOt, Fe3Fet = num(px, "FeOt"), num(px, "Fe3Fet")
    f = 2 * OXIDE_MW["FeO"] / OXIDE_MW["Fe2O3"]          # wt% Fe2O3 -> wt% FeO
 
    if fe_mode == "total":
        out["FeO"] = FeOt.where(FeOt > 0, out["FeO"] + f * out["Fe2O3"])
        out["Fe2O3"] = 0.0
    else:
        no_fe = (out["FeO"] + out["Fe2O3"]) <= 0
        out.loc[no_fe, "FeO"] = (FeOt * (1 - Fe3Fet))[no_fe]
        out.loc[no_fe, "Fe2O3"] = (FeOt * Fe3Fet / f)[no_fe]
    return out

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

def ensure_oxide_columns(df):
    """
    Ensures all keys in oxide_dict exist as columns in df.
    Fills missing columns with 0.0 while retaining any additional pre-existing columns.
    """
    for col in df.columns:
        for oxide in OXIDE_MW:
            if col.startswith(f"{oxide}_"):
                # Extract the suffix portion (e.g., '_Liq')
                suffix = (col[len(oxide):])
                break

    expected_cols = [f"{oxide}{suffix}" for oxide in OXIDE_MW]
    missing_cols = [col for col in expected_cols if col not in df.columns]
    
    if missing_cols:
        all_columns = list(df.columns) + missing_cols
        df = df.reindex(columns=all_columns, fill_value=0.0)

        if f"FeO{suffix}" in missing_cols and f"Fe2O3{suffix}" in missing_cols:
            if f"FeOt{suffix}" not in missing_cols and f"Fe3Fet{suffix}" not in missing_cols:
                df[f"FeO{suffix}"] = df[f"FeOt{suffix}"]*(1-df[f"Fe3Fet{suffix}"])
                df[f"Fe2O3{suffix}"] = df[f"FeOt{suffix}"]*((159.69/2)/71.844)*(df[f"Fe3Fet{suffix}"])

    return df

#### liquid ####

def liq_mole_fractions(liquid, fe_mode="as_given"):
    """Oxide mole fractions of a melt; columns are matched by name, so column order does not matter."""
    ox = _oxides(liquid, LIQ_SUFFIX, fe_mode)
    moles = ox / pd.Series(OXIDE_MW)[ox.columns]
    return moles.div(moles.sum(axis=1), axis=0)

def any_cation_fractions(df, return_total_cations=False):
    data = df.copy()

    suffix = None
    for col in data.columns:
        if '_' in col:
            suffix = col[col.find('_'):]
            break

    if suffix is not None:
        if f"FeOt{suffix}" in data.columns:
            data = data.drop(columns=[f'FeOt{suffix}'])
        if f"Fe3Fet{suffix}" in data.columns:
            data = data.drop(columns=[f'Fe3Fet{suffix}'])
        data.columns = data.columns.str.replace(suffix, "", regex=False)
    else:
        if "FeOt" in data.columns:
            data = data.drop(columns=['FeOt'])
        if "Fe3Fet" in data.columns:
            data = data.drop(columns=['Fe3Fet'])

    oxide_mw = pd.Series(OXIDE_MW)
    cation_n = pd.Series(CATION_N)

    moles = data.div(oxide_mw)
    cat = moles.mul(cation_n)
    total_cations = cat.sum(axis=1)
    cat_frac = cat.div(total_cations, axis=0)

    if return_total_cations:
        return cat_frac, total_cations
    return cat_frac

#### clinopyroxene ####

def strip(df, pattern):
    out = df.copy()
    out.columns = [pattern.sub("", c) if isinstance(c, str) else c for c in out.columns]
    return out.loc[:, ~out.columns.duplicated()]
 
 
def num(df, name):
    """Column as float with blanks/NaN/text -> 0; all-zero Series if the column is absent."""
    if name in df.columns:
        return pd.to_numeric(df[name], errors="coerce").fillna(0.0).astype(float)
    return pd.Series(0.0, index=df.index)
 
 
def px_oxides(pyx, fe_mode="as_given"):
    """Clean wt% table with exactly the oxides the normalisation uses.

    * strips _Cpx/_Opx/_clinopyroxene/_orthopyroxene suffixes (any case)
    * ignores everything that isn't an oxide (Sample_ID, P_kbar, ..._Liq columns)
    * missing oxide columns and blank/NaN cells -> 0 (what the Excel macros do)
    * resolves iron row by row:
        fe_mode="as_given": FeO/Fe2O3 used where given; rows with neither fall
                            back to FeOt (split using Fe3Fet if supplied)
        fe_mode="total":    all Fe as FeO (FeOt if > 0, else FeO + 0.8998*Fe2O3)
    """
    return _oxides(pyx, PX_SUFFIX, fe_mode)
 
 
def _mol(ox):
    cols = ox.columns
    moles = ox / pd.Series(OXIDE_MW)[cols]
    return moles * pd.Series(CATION_N)[cols], moles * pd.Series(OXIDE_N)[cols]
 
 
def px_cation_fraction(pyx, fe_mode="as_given"):
    cat, _ = _mol(px_oxides(pyx, fe_mode))
    return cat.div(cat.sum(axis=1), axis=0)
 
 
def px_apfu(pyx, charge_balance=False, fe_mode="as_given"):
    """Cations per 6 oxygens (columns keep oxide names, e.g. apfu['Al2O3'] = Al)."""
    cat, oxy = _mol(px_oxides(pyx, fe_mode))
    return cat.div(oxy.sum(axis=1), axis=0) * 6

# def px_cation_fraction(pyx):
#     px = pyx.copy()

#     suffix = None
#     for col in px.columns:
#         if '_' in col:
#             suffix = col[col.find('_'):]
#             break
#     suffix = suffix or '_Cpx'

#     if f"FeOt{suffix}" in px.columns:
#         px = px.drop(columns = [f"FeOt{suffix}"])
#     if f"FeOt" in px.columns:
#         px = px.drop(columns = ['FeOt'])

#     if f"Fe3Fet{suffix}" in px.columns:
#         px = px.drop(columns = [f'Fe3Fet{suffix}'])
#     if "Fe3Fet" in px.columns:
#         px = px.drop(columns = ['Fe3Fet'])

#     pattern = r"_(Cpx|Opx|clinopyroxene|orthopyroxene).*"
#     px.columns = px.columns.str.replace(pattern, "", regex=True)

#     moles = px.div(pd.DataFrame.from_dict(OXIDE_MW, orient = 'index').T.loc[0], axis = 1)
#     # moles_O = moles.multiply(pd.DataFrame.from_dict(OXIDE_N, orient = 'index').T.loc[0], axis = 1)
#     moles_cat = moles.multiply(pd.DataFrame.from_dict(CATION_N, orient = "index").T.loc[0], axis = 1)

#     cation_frac = moles_cat.div(moles_cat.sum(axis = 1),axis = 0)

#     return cation_frac

# def px_apfu(pyx, charge_balance = False):
#     px = pyx.copy()

#     suffix = None
#     for col in px.columns:
#         if '_' in col:
#             suffix = col[col.find('_'):]
#             break
#     suffix = suffix or '_Cpx'

#     if f"FeOt{suffix}" in px.columns:
#         px = px.drop(columns = [f"FeOt{suffix}"])
#     if f"FeOt" in px.columns:
#         px = px.drop(columns = ['FeOt'])

#     if f"Fe3Fet{suffix}" in px.columns:
#         px = px.drop(columns = [f'Fe3Fet{suffix}'])
#     if "Fe3Fet" in px.columns:
#         px = px.drop(columns = ['Fe3Fet'])

#     pattern = r"_(Cpx|Opx|clinopyroxene|orthopyroxene).*"
#     px.columns = px.columns.str.replace(pattern, "", regex=True)

#     moles = px.div(pd.DataFrame.from_dict(OXIDE_MW, orient = 'index').T.loc[0], axis = 1)
#     moles_O = moles.multiply(pd.DataFrame.from_dict(OXIDE_N, orient = 'index').T.loc[0], axis = 1)
#     moles_cat = moles.multiply(pd.DataFrame.from_dict(CATION_N, orient = "index").T.loc[0], axis = 1)

#     sum_O = moles_O.sum(axis = 1)

#     cat_6_ox = moles_cat.div(sum_O, axis = 0)*6

#     return cat_6_ox