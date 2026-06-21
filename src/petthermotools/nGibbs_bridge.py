"""
nGibbs_bridge.py
----------------
Single instantiation point for nGibbs emulators and the input-table builder
that translates PTT-style multi_path arguments into a 2D batch array for
``emulator.ForwardMB``.

Imported by both Path.py (for the actual dispatch) and __init__.py (so users
can inspect ``ptt._nGibbs_models``).  Neither file contains its own registry,
eliminating duplication and the need for a circular-import workaround.

"""

from __future__ import annotations

import contextlib
import io
import os
import sys

import numpy as np
import pandas as pd
from petthermotools.MELTS import calc_phase_props_MELTS_parallel, calc_phase_props_MELTS # Property Calculators

_FE2O3_TO_FEO = 0.8998

_nGibbs_models = {}

try:
    import ngibbs as _ngibbs_pkg
    from ngibbs.config.constants import ALLOWED_MELTS_OXIDES, REQUIRED_MELTS_OXIDES
    from ngibbs.utils.math_utils import grid_sample_explicit, grid_sample # Helpers. Grid Sample used elsewhere. 
    from ngibbs.engine.API import MELTS102EmulatorCPU, MELTS120EmulatorCPU
    _nGibbs_models['nMELTSv1.0.2'] = MELTS102EmulatorCPU
    _nGibbs_models['nMELTSv1.2.0'] = MELTS120EmulatorCPU

except:
    pass # Installation hint in module __init__.py

def find_oxide_cols(label):
    """
    Identify which oxide is represented by the string label.

    Returns the oxide name (e.g. 'SiO2', 'FeO') if found, else None.

    Matching rules (in priority order):
      1. Exact match: label == oxide  (e.g. 'H2O', 'SiO2')
      2. Fe3Fet check: 'fe3fet' anywhere in label (case-insensitive)
      3. Oxide prefix with underscore: label starts with oxide + '_'
         (e.g. 'H2O_Liq', 'SiO2_Liq', 'H2O_init')
      4. For total iron, also allow oxide prefix with 't' or 'T' suffix: label starts with oxide + 't' or oxide + 'T'

    Deliberately does NOT match columns where the oxide name is merely a
    substring (e.g. 'XH2O_fl_VESIcal' should NOT map to 'H2O').
    """
    # 1. Exact match
    if label in ALLOWED_MELTS_OXIDES:
        return label
    # 2. Fe3Fet (covers 'Fe3Fet_Liq', 'Fe3FeT', 'fe3fet_init', …)
    if 'fe3fet' in label.lower():
        return 'Fe3FeT'
    # 3. Oxide prefix with underscore separator
    for oxide in ALLOWED_MELTS_OXIDES:
        if label.startswith(oxide + '_'):
            return oxide
        if label.startswith(oxide + 't'):
            return oxide
        if label.startswith(oxide + 'T'):
            return oxide
    return None


### FUNCTION TO GET P, T, fO2 vectos from crazy ptt input args. Make it a **kwargs type

def handle_kwargs(modelname="nMELTSv1.0.2", CO2_OK=False, verbose=False, **kwargs):
    """
    processes kwargs to create a dict of compositional overrides. One such override may contain multiple values to trigger a grid calculation
    Warns users about functionalities that are not supported by nGibbs such as phase suppression
    """
    overrides = {}
    for key, val in kwargs.items():
        if val is None:
            continue
        if key == 'fO2_buffer':
            if val.lower() not in ['fmq', 'qfm']:
                if verbose:
                    print(f"[WARNING] Ignoring fO2 buffer '{val}' as it is not recognised/supported by nGibbs. nGibbs only supports 'FMQ' and 'QFM' buffers. Other buffers will be supported in future versions")
                continue
            else:
                overrides['fO2'] = 0 #If offset is not provided, default to FMQ (0 log units from FMQ)
                continue
        if key == 'fO2_offset':
            overrides['fO2'] = val #If buffer is provided, apply offset to get logfO2 relative to FMQ
            continue
        if val is not None:
            if  "suppress" in key.lower():
                if verbose:
                    print(f"Phase suppression not supported by nGibbs, although new models may be trained for your use case. Email: benjthyer@gmail.com")
                    print(f"Phases in the {modelname} isothermal Cr-bearing model:")
                    print(_nGibbs_models[modelname].cr.isothermal_emulator.ml_indexer.all_phases)
                 #This can be specified further when not all emulators share the same phase. But they probably should...
            else:
                col = find_oxide_cols(key)
                if col is not None:
                    if (col != 'CO2') or CO2_OK:
                        overrides[col] = val
                    else:
                        if verbose:
                            print(f"Warning: CO2 is only supported by nMELTSv1.2.0. Ignoring CO2 argument {key}.")
                else:
                    if verbose:
                        print(f"Warning: nGibbs is ignoring Unrecognised keyword argument '{key}'")

    if 'fO2' in overrides and 'Fe3Fet' in overrides:
        raise ValueError("Simultaneous fO2 and Fe3Fet overrides are mutually exclusive!")

    # CO2 is only supported by nMELTSv1.2.0
    if 'CO2' in overrides and 'v1.2.0' not in modelname:
        if verbose:
            print(f"Warning: CO2 is only supported by nMELTSv1.2.0. Ignoring CO2={overrides['CO2']} for {modelname}.")
        del overrides['CO2']

    # Format overrides as lists, check that no more than one is multi-valued
    multi_run = False
    for key, val in overrides.items():
        listval = condition_inputs_as_lists(val)
        if len(listval) > 1 and key != 'fO2': # fO2 is handled differently because sometimes it is needed for parallel calculation mode
            if multi_run:
                raise ValueError(f"Multiple multi-valued overrides found.\n" 
                                 "Only one multi-valued override arg (CO2, H2O, Fe3Fet liq or init) " 
                                 "is allowed to trigger a grid calculation.")
            multi_run = True
        overrides[key] = listval

    return overrides

def extract_composition_from_table(table: pd.DataFrame, overrides={}, CO2_OK=False, verbose=False):

    """
    Convert between PTT-style oxide names (e.g. 'SiO2_Liq') and nGibbs-style
    bare names (e.g. 'SiO2').

    Outputs headers and a numpy array table with their values

    Overrides are taken from arguments and overwrite the values in the table, or are added as new columns if not already present
    """

    for key, val in overrides.items():
        if len(val) > 1:
            raise ValueError(f"Multi-valued override for '{key}' found in arguments.\n"
                             "Multi-valued overrides trigger grid calculations, which are supported only for a single composition.\n"
                             "For parallel, custom calculations, place these free-to-vary values in a composition table")
    keepIDX = []
    headers = []
    FeR = None
    for idx, header in enumerate(table.columns):
        oxide = find_oxide_cols(header)
        if oxide is not None:
            if oxide == 'Fe3Fet':
                if FeR is not None:
                    raise ValueError("Multiple Fe3Fet columns found in composition table.")
                else:
                    if 'Fe3Fet' in overrides:
                        FeR = overrides['Fe3Fet'][0]
                    else:
                        FeR = table.to_numpy(dtype=np.float32)[:, idx]
                    continue
            if (oxide != 'CO2') or CO2_OK:
                keepIDX.append(idx)
                headers.append(oxide)
            else:
                if verbose:
                    print(f"Warning: CO2 is only supported in nGibbs by nMELTSv1.2.0. Ignoring CO2 column in composition table.")
                continue

    for oxide, override_val in overrides.items():
        if oxide != 'Fe3Fet':
            if oxide in headers:
                idx = headers.index(oxide)
                table.iloc[:, idx] = override_val[0]
            else:
                table[oxide] = override_val[0]
                keepIDX.append(table.columns.get_loc(oxide))
                headers.append(oxide)

    out_table = table.iloc[:, keepIDX].to_numpy(dtype=np.float32)

    if FeR is None and 'Fe3Fet' in overrides:
        FeR = overrides['Fe3Fet'][0]

    if FeR is not None:
        if 'Fe2O3' in headers:            
            raise ValueError("Fe2O3 column already exists in composition table. Cannot apply Fe3Fet override.")
        # Recast Fe3Fet into Fe2O3 and append to the table 
        FeOt = out_table[:, headers.index('FeO')]
        Fe2O3 = FeOt * FeR / _FE2O3_TO_FEO
        out_table = np.column_stack([out_table, Fe2O3])
        out_table[:, headers.index('FeO')] = FeOt * (1.0 - FeR)
        headers.append('Fe2O3')

    missing = set(REQUIRED_MELTS_OXIDES) - set(headers)
    assert not missing, f"Missing nGibbs-required oxides in composition table: {missing}"
    assert len(headers) == len(set(headers)), f"Duplicate oxide columns found in composition table: {headers}"

    return headers, out_table


def extract_composition_1D(vector: pd.Series | dict, overrides={}, CO2_OK=False, verbose=False):

    """
    Convert between PTT-style oxide names (e.g. 'SiO2_Liq') and nGibbs-style
    bare names (e.g. 'SiO2').

    Outputs headers and a list of values in the same order as the headers

    Overrides are taken from arguments and overwrite the values in the table, or are added as new columns if not already present
    """

    values = []
    headers = []
    FeR = None
    for key, val in vector.items():
        oxide = find_oxide_cols(key)
        if oxide is not None:
            if oxide == 'Fe3Fet':
                if FeR is not None:
                    raise ValueError("Multiple Fe3Fet columns found in composition table.")
                else:
                    if 'Fe3Fet' in overrides:
                        FeR = overrides['Fe3Fet'][0]
                    else:
                        FeR = val
                    continue
            if (oxide != 'CO2') or CO2_OK:
                values.append(val)
                headers.append(oxide)
            else:
                if verbose:
                    print(f"Warning: CO2 is only supported in nGibbs by nMELTSv1.2.0. Ignoring CO2 column in composition table.")
                continue

    for oxide, override_val in overrides.items():
        if oxide != 'Fe3Fet':
            if oxide in headers:
                idx = headers.index(oxide)
                values[idx] = override_val[0]
            else:
                values.append(override_val[0])
                headers.append(oxide)

    if FeR is None and 'Fe3Fet' in overrides:
        FeR = overrides['Fe3Fet'][0]

    if FeR is not None:
        if 'Fe2O3' in headers:            
            raise ValueError("Fe2O3 column already exists in composition table. Cannot apply Fe3Fet override.")
        # Recast Fe3Fet into Fe2O3 and append to vector
        FeOt = values[headers.index('FeO')]
        Fe2O3 = FeOt * FeR / _FE2O3_TO_FEO
        values.append(Fe2O3)
        values[headers.index('FeO')] = FeOt * (1.0 - FeR)
        headers.append('Fe2O3')

    missing = set(REQUIRED_MELTS_OXIDES) - set(headers)
    assert not missing, f"Missing nGibbs-required oxides in composition: {missing}"
    assert len(headers) == len(set(headers)), f"Duplicate oxide columns found in composition: {headers}"

    return headers, np.array(values, dtype=np.float32).reshape(1, -1)

def build_PT_vectors(
    T_C=None, T_path_C=None, T_start_C=None, T_end_C=None, dt_C=None,
    T_min_C=None, T_max_C=None, T_num=None,
    P_bar=None, P_path_bar=None, P_start_bar=None, P_end_bar=None, dp_bar=None,
    P_min_bar=None, P_max_bar=None, P_num=None,
):
    """
    Resolve PTT-style T and P arguments into 1D numpy arrays and determine
    whether the result is a correlated PT path or an independent P×T grid.

    Resolution priority for each axis:
      *_path   >  *_start / *_end / d*  |  *_min / *_max / *_num  >  scalar (*_C or *_bar)

    Parameters
    ----------
    T_C : float or array-like
        Fixed temperature(s) in °C.
    T_path_C : array-like
        Explicit temperature path in °C (used as-is).
    T_start_C, T_end_C, dt_C : float
        Generate a temperature path with ``np.arange``.  Direction is inferred
        from start vs end; the sign of dt_C is ignored.
    T_min_C, T_max_C, T_num : float, float, int
        Generate a temperature vector with ``np.linspace(T_min_C, T_max_C, T_num)``.
    P_bar : float or array-like
        Fixed pressure(s) in bar.
    P_path_bar : array-like
        Explicit pressure path in bar (used as-is).
    P_start_bar, P_end_bar, dp_bar : float
        Generate a pressure path with ``np.arange``.  Direction is inferred
        from start vs end; the sign of dp_bar is ignored.
    P_min_bar, P_max_bar, P_num : float, float, int
        Generate a pressure vector with ``np.linspace(P_min_bar, P_max_bar, P_num)``.

    Returns
    -------
    P_vec, T_vec : np.ndarray
        1D float32 arrays of pressures and temperatures.
    is_path : bool
        True  → P_vec and T_vec are correlated step-by-step (PT path or single
                point).  run_nGibbs will broadcast a scalar side to match the
                other and build a path input array.
        False → P_vec and T_vec are independent axes (grid mode).
                run_nGibbs will form a P × T cross-product via
                grid_sample_explicit.

    is_path logic
    -------------
    Grid mode requires ALL of the following:
      • P has more than one value (from any input form)
      • T has more than one value (from any input form, including range syntax)
      • The user did NOT supply both T_path_C and P_path_bar simultaneously
        (that would be an explicit correlated path)
    Everything else is path/single-point mode.
    """

    def _any_multi(x):
        """True if x resolves to more than one element."""
        if x is None:
            return False
        try:
            return len(np.atleast_1d(x)) > 1
        except TypeError:
            return False

    # ── is_path determination (uses pre-resolution inputs) ───────────────────
    has_explicit_correlated_path = (T_path_C is not None and P_path_bar is not None)
    _P_will_be_multi = (
        _any_multi(P_bar) or _any_multi(P_path_bar) or _any_multi(P_start_bar) or
        (P_min_bar is not None and P_max_bar is not None and P_num is not None and int(P_num) > 1)
    )
    _T_will_be_multi = (
        _any_multi(T_C) or
        _any_multi(T_path_C) or
        (T_start_C is not None and T_end_C is not None and dt_C is not None) or
        (T_min_C is not None and T_max_C is not None and T_num is not None and int(T_num) > 1)
    )

    if has_explicit_correlated_path:
        is_path = True
    elif _P_will_be_multi and _T_will_be_multi:
        is_path = False  # grid mode: P × T
    else:
        is_path = True   # path/single-point mode

    # ── Temperature ──────────────────────────────────────────────────────────
    if T_path_C is not None:
        T_vec = np.asarray(T_path_C, dtype=np.float32).ravel()
    elif T_start_C is not None and T_end_C is not None and dt_C is not None:
        if dt_C == 0:
            raise ValueError("dt_C must be non-zero.")
        step = abs(float(dt_C)) if T_end_C >= T_start_C else -abs(float(dt_C))
        T_vec = np.arange(T_start_C, T_end_C + step * 0.5, step, dtype=np.float32)
    elif T_min_C is not None and T_max_C is not None and T_num is not None:
        T_vec = np.linspace(T_min_C, T_max_C, int(T_num), dtype=np.float32)
    elif T_C is not None:
        T_vec = np.atleast_1d(np.asarray(T_C, dtype=np.float32)).ravel()
    else:
        raise ValueError(
            "No temperature specified. Provide T_C, T_path_C, T_start_C/T_end_C/dt_C, or T_min_C/T_max_C/T_num."
        )

    # ── Pressure ─────────────────────────────────────────────────────────────
    if P_path_bar is not None:
        P_vec = np.asarray(P_path_bar, dtype=np.float32).ravel()
    elif P_start_bar is not None and P_end_bar is not None and dp_bar is not None:
        if dp_bar == 0:
            raise ValueError("dp_bar must be non-zero.")
        step = abs(float(dp_bar)) if P_end_bar >= P_start_bar else -abs(float(dp_bar))
        P_vec = np.arange(P_start_bar, P_end_bar + step * 0.5, step, dtype=np.float32)
    elif P_min_bar is not None and P_max_bar is not None and P_num is not None:
        P_vec = np.linspace(P_min_bar, P_max_bar, int(P_num), dtype=np.float32)
    elif P_bar is not None:
        P_vec = np.atleast_1d(np.asarray(P_bar, dtype=np.float32)).ravel()
    else:
        raise ValueError(
            "No pressure specified. Provide P_bar, P_path_bar, P_start_bar/P_end_bar/dp_bar, or P_min_bar/P_max_bar/P_num."
        )

    return P_vec, T_vec, is_path


def condition_inputs_as_lists(inputs):
    if isinstance(inputs, pd.Series):
        return list(inputs)
    if isinstance(inputs, np.ndarray):
        if inputs.size == 1:
            return [inputs.item()]
        else:
            return inputs.tolist()
    if isinstance(inputs, list):
        return inputs
    if isinstance(inputs, int) or isinstance(inputs, float):
        return [inputs]
    raise ValueError(f"Unsupported input type: {type(inputs)}. Expected pd.Series, np.ndarray, integer/float scalar, or list.")


def _is_liquid_only_mask(mass_g_df):
    """Boolean array: True where only liquid1 (and optionally fluid1) carry non-zero mass."""
    solid_cols = [c for c in mass_g_df.columns if c not in ('liquid1', 'fluid1')]
    if not solid_cols:
        return np.ones(len(mass_g_df), dtype=bool)
    return (mass_g_df[solid_cols].values <= 1e-6).all(axis=1)


def trim_superliquidus(nGibbs_output, grouping_column=None, tableNo=None):
    """
    Remove superliquidus rows, keeping the liquidus boundary and all sub-liquidus rows.

    For each group, finds the lowest-T all-liquid row (the liquidus) and drops all
    rows with higher T. The liquidus row itself and every row below it (where crystals
    are present) are retained. If no all-liquid row exists in a group, all rows are kept.

    grouping_column : None
        Treat the entire output as one group.
    grouping_column : 'Pressure(System_main)'
        Group by the P_bar column of the Conditions table.
    grouping_column : anything else
        Group by the provided tableNo integer array (one int per row).
    tableNo : ndarray of int, shape (n_rows,)
        Required when grouping_column is not None or 'Pressure(System_main)'.
    """
    T_vals = nGibbs_output['Conditions']['T_C'].values
    liquid_mask = _is_liquid_only_mask(nGibbs_output['mass_g'])
    n = len(T_vals)
    keep = np.zeros(n, dtype=bool)

    if grouping_column is None:
        liq_idx = np.where(liquid_mask)[0]
        if len(liq_idx):
            liq_T = T_vals[liq_idx[np.argmin(T_vals[liq_idx])]]
            keep = T_vals <= liq_T
        else:
            keep[:] = True  # no liquidus found — keep everything
    elif grouping_column == 'Pressure(System_main)':
        group_arr = nGibbs_output['Conditions']['P_bar'].values
        for g in np.unique(group_arr):
            g_bool = group_arr == g
            g_idx = np.where(g_bool)[0]
            liq_in_g = g_idx[liquid_mask[g_idx]]
            if len(liq_in_g):
                liq_T = T_vals[liq_in_g[np.argmin(T_vals[liq_in_g])]]
                keep |= g_bool & (T_vals <= liq_T)
            else:
                keep |= g_bool  # no liquidus found in this P group — keep all its rows
    else:
        if tableNo is None:
            raise ValueError(
                f"trim_superliquidus: tableNo is required when grouping_column='{grouping_column}'"
            )
        for g in np.unique(tableNo):
            g_bool = tableNo == g
            g_idx = np.where(g_bool)[0]
            liq_in_g = g_idx[liquid_mask[g_idx]]
            if len(liq_in_g):
                liq_T = T_vals[liq_in_g[np.argmin(T_vals[liq_in_g])]]
                keep |= g_bool & (T_vals <= liq_T)
            else:
                keep |= g_bool  # no liquidus found — keep all rows of this group

    for key, val in nGibbs_output.items():
        if isinstance(val, pd.DataFrame):
            nGibbs_output[key] = val.loc[keep].reset_index(drop=True)

    return nGibbs_output


def _concat_ptt_tables(tables):
    """Row-concatenate a list of PTT output dicts (each a dict of DataFrames)."""
    all_keys = set().union(*(t.keys() for t in tables))
    return {
        key: pd.concat([t[key] for t in tables if key in t], ignore_index=True)
        for key in all_keys
    }


def run_nGibbs(P, T, comp, overrides={}, Model="nMELTSv1.2.0NoProp", is_path=False, label=None, find_liquidus=False, verbose=False):
    """
    Run nGibbs for a 1D-2D PTX grid.

    is_path=True  → P and T are correlated step-by-step (PT path / single point).
                    Scalar P or T is broadcast to match the length of the other.
    is_path=False → P and T are independent axes (grid mode).
                    Total rows = len(P) × len(T) [× len(fO2)].

    Only supported for isothermal calculations for now. Isentropic calculations will require more engineering

    """

    if Model.endswith('NoProp'):
        modelname = Model[:-6] 
        run_props = False
    else:
        modelname = Model
        run_props = True

    if modelname == 'nMELTSv1.2.0':
        CO2_OK = True
    else:        
        CO2_OK = False

    grouping_column = None # used to organize pandas tables outputs
    fO2 = overrides.pop('fO2', None)
    grid_trigger = False
    grid_name = None
    for key, val in overrides.items():
        if len(val) > 1:
            if grid_trigger:
                raise ValueError(f"Multiple multi-valued overrides found.\n" 
                                    "Only one multi-valued override arg (CO2, H2O, Fe3FeT liq or init) " 
                                    "is allowed to trigger a grid calculation")
            grid_name = key
            grid_trigger = True
        

    if modelname not in _nGibbs_models:
        raise ValueError(f"Model '{Model}' not recognised. Available models: {list(_nGibbs_models.keys())}")

    if isinstance(comp, pd.DataFrame): # 1D (parallel) mode, runs = rows in comp table, with explicitly specified compositions for each PT point.
        if verbose: print(f"Running nGibbs in parallel mode for {comp.shape[0]} parallel compositions")
        if grid_trigger:
            raise ValueError(f"Multi-valued override '{grid_name}' found in arguments.\n"
                             "Multi-valued compositional overrides trigger grid calculations, which are supported only for a single composition.\n"
                             "For parallel, custom calculations, place these free-to-vary values in the composition table instead")
        headers, comp_array = extract_composition_from_table(comp, overrides=overrides, CO2_OK=CO2_OK, verbose=verbose) # Will error if any override is multi-valued!!
        P_vec = np.array(P, dtype=np.float32).reshape(-1, 1)
        T_vec = np.array(T, dtype=np.float32).reshape(-1, 1)
        lengths = [len(P), len(T), comp_array.shape[0]]

        if fO2 is not None:
            fO2_vec = np.array(fO2, dtype=np.float32).reshape(-1, 1)
            lengths.append(len(fO2))
            headers = ["Pressure(System_main)", "Temperature(System_main)", 'logfO2-QFM(System_main)'] + headers

            if not np.all(np.isin(np.unique(lengths), [1, max(lengths)])):      
                raise ValueError(f"Inconsistent input lengths: P={len(P)}, T={len(T)}, X={comp_array.shape[0]}, fO2={len(fO2)}.\n" 
                                "For parallel runs of variable composition all inputs must have length 1 or the same length as the composition matrix,\n" 
                                "or X must be a series/dictionary for constant-composition grid calculations")

            P_vec, T_vec, fO2_vec, comp_array= [
                np.repeat(obj, max(lengths), axis=0) if obj.shape[0] == 1 else obj
                for obj in [P_vec, T_vec, fO2_vec, comp_array] ]
        else:
            headers = ["Pressure(System_main)", "Temperature(System_main)"] + headers

            if not np.all(np.isin(np.unique(lengths), [1, max(lengths)])):      
                raise ValueError(f"Inconsistent input lengths: P={len(P)}, T={len(T)}, X={comp_array.shape[0]}.\n" 
                                "For parallel runs of variable composition all inputs must have length 1 or the same length as the composition matrix,\n" 
                                "or X must be a series/dictionary for constant-compositiongrid calculations")
            
            P_vec, T_vec, comp_array = [
                np.repeat(obj, max(lengths), axis=0) if obj.shape[0] == 1 else obj
                for obj in [P_vec, T_vec, comp_array] ]

        
        
        input_array = np.hstack([P_vec, T_vec, comp_array]) if fO2 is None else np.hstack([P_vec, T_vec, fO2_vec, comp_array])

    elif isinstance(comp, pd.Series) or isinstance(comp, dict): # PT Path OR 2D/3D grid mode (or 1D if only one of PTfO2 varies)
        #1D outputs: single ptt out. 
        #2D outputs: sort tables by pressure, or 
        #3D outputs:
        headers, comp_array = extract_composition_1D(comp, overrides=overrides, CO2_OK=CO2_OK, verbose=verbose)

        P_vec = condition_inputs_as_lists(P)
        T_vec = condition_inputs_as_lists(T)

        grid_warning_text = f" at {len(overrides[grid_name])} values for {grid_name}" if grid_trigger else ""

        if fO2 is not None:
            fO2_vec = condition_inputs_as_lists(fO2) 

            headers = ["Pressure(System_main)", "Temperature(System_main)", 'logfO2-QFM(System_main)'] + headers
        else:
            headers = ["Pressure(System_main)", "Temperature(System_main)"] + headers

        if is_path: # PT path / single-point mode
            # Broadcast scalar P or T to match the length of the other axis.
            if len(P_vec) == 1 and len(T_vec) > 1:
                P_vec = np.repeat(P_vec, len(T_vec))
            elif len(T_vec) == 1 and len(P_vec) > 1:
                T_vec = np.repeat(T_vec, len(P_vec))
            PT_arr = np.column_stack([P_vec, T_vec])
            if fO2 is None: # Simple case: one PT path, no fO2 variation
                if verbose: print(f"Running nGibbs in PT path mode: {len(P_vec)} points{grid_warning_text}")
                input_array = np.hstack([PT_arr, np.tile(comp_array, (len(P_vec), 1))])
            else: # More complex: One PT path at variable fO2
                if verbose: print(f"Running nGibbs in PT path mode for: {len(P_vec)} points for {len(fO2_vec)} fO2 values{grid_warning_text}")
                if len(fO2_vec) == 1:
                    # Single fO2: broadcast across path (simple hstack)
                    fO2_col = np.full((len(P_vec), 1), fO2_vec[0], dtype=np.float32)
                    input_array = np.hstack([PT_arr, fO2_col, np.tile(comp_array, (len(P_vec), 1))])
                else:
                    # Multiple fO2: cross PT path with fO2 values
                    # grid_sample_explicit requires 2-D column-vector params
                    fO2_arr = np.array(fO2_vec, dtype=np.float32).reshape(-1, 1)
                    input_array = grid_sample_explicit([PT_arr, fO2_arr, comp_array])
                    grouping_column = 'logfO2-QFM(System_main)'


        else: # Grid mode: P and T are independent axes. fO2 or a composition override adds a third axis.
            # grid_sample_explicit requires 2-D column-vector params
            P_arr = np.array(P_vec, dtype=np.float32).reshape(-1, 1)
            T_arr = np.array(T_vec, dtype=np.float32).reshape(-1, 1)
            if fO2 is not None:
                if (len(fO2_vec) != 1):
                    grouping_column = 'logfO2-QFM(System_main)'
                    if grid_trigger:
                        raise ValueError(f"Attempted to do grid calculation with 3 or more dimensions: P, T, fO2, and {grid_name}. Not yet supported.")
                else:
                    grouping_column = 'Pressure(System_main)'
                if verbose: print(f"Running nGibbs in grid mode of {len(P_vec)} P points, {len(T_vec)} T points, and {len(fO2_vec)} fO2 values")
                fO2_arr = np.array(fO2_vec, dtype=np.float32).reshape(-1, 1)
                input_array = grid_sample_explicit([P_arr, T_arr, fO2_arr], comp_array)
            else:
                if verbose: print(f"Running nGibbs in grid mode of {len(P_vec)} P points, {len(T_vec)} T points{grid_warning_text}")
                grouping_column = 'Pressure(System_main)'
                input_array = grid_sample_explicit([P_arr, T_arr], comp_array)


        if grid_trigger: # Prepare multi-runs at different specified conditions
            if grid_name == 'Fe3FeT': # The ugliest op. First is already done and Fe2O3 columns already exist.  
                ferricIDX = headers.index('Fe2O3') # Get total FeOT
                ferrousIDX = headers.index('FeO')
                FeOT = input_array[:, ferrousIDX] + input_array[:, ferricIDX ] * _FE2O3_TO_FEO # Should be constant

            base_input_array = input_array.copy() 
            for i, override_val in enumerate(overrides[grid_name]):
                if i == 0:
                    continue # First run is already done in the base input array, so skip to next override value
                new_arr = base_input_array.copy()

                if grid_name == 'Fe3FeT':
                    FeR = override_val
                    Fe2O3 = FeOT * FeR / _FE2O3_TO_FEO
                    new_arr[:, ferricIDX] = Fe2O3
                    new_arr[:, ferrousIDX] = FeOT * (1.0 - FeR)
                else:
                    col_idx = headers.index(grid_name)
                    new_arr[:, col_idx] = override_val

                input_array = np.vstack([input_array, new_arr])

            grouping_column = grid_name # Add grouping column to organize outputs by override variable in the final table. This is used in the ptt output processing step to split tables by override variable and add a column with the override variable value for each run.
            

    else:
        raise ValueError(f"Unsupported composition input type: {type(comp)}. Expected pd.DataFrame, pd.Series, or dict.")


    if input_array.shape[0] > 500_000:
        raise ValueError(f"Input array has {input_array.shape[0]} rows. This errors because we think you did this on accident.")


    # Prepare input array for nGibbs!
    if fO2 is not None and 'Fe2O3' in headers: 
        if verbose:
            print("Warning: fO2 is specified in overrides, so iron will be recalculated as total iron from Fe2O3 and FeO.\n"
                  "Table Fe3FeT values are ignored.")
        ferricIDX = headers.index('Fe2O3')
        ferrousIDX = headers.index('FeO')
        input_array[:, ferrousIDX] += input_array[:, ferricIDX] * _FE2O3_TO_FEO
        input_array = np.delete(input_array, ferricIDX, axis=1)
        headers.remove('Fe2O3')


    colsums = np.sum(input_array, axis=0)
    zerocols = np.where(colsums == 0)[0]
    delIDX = []
    if np.any(zerocols):
        #### Delete if not pressure or fO2.
        for idx in zerocols:
            if headers[idx] not in ["Pressure(System_main)", "Temperature(System_main)", 'logfO2-QFM(System_main)']:
                delIDX.append(idx)

    if delIDX: # Need to delete absent stuff (esp. Cr to trigger NoCr model)
        input_array = np.delete(input_array, delIDX, axis=1)
        headers = [h for i, h in enumerate(headers) if i not in delIDX]

    nGibbs_output = _nGibbs_models[modelname].ForwardMB( # Support ForwardNN for phase diagrams? May be cleaner. Speed is no issue. 
            input_array, headers=headers, outputs=["ptt_out"]
        )["ptt_out"]

    _is_parallel = isinstance(comp, pd.DataFrame)

    if grouping_column is not None:
        # Compute row-to-group assignment. P_bar is readable from Conditions; fO2 and
        # composition override values are not stored there so are inferred from input shape.
        total_rows = nGibbs_output['Conditions'].shape[0]
        if grouping_column == 'Pressure(System_main)':
            unique_vals, tableNo = np.unique(
                nGibbs_output['Conditions']['P_bar'], return_inverse=True
            )
        elif grouping_column == 'logfO2-QFM(System_main)':
            k = len(fO2_vec)
            tableNo = np.tile(np.arange(k), total_rows // k)
            unique_vals = fO2_vec
        else:
            n_overrides = len(overrides[grid_name])
            n_base = total_rows // n_overrides
            tableNo = np.repeat(np.arange(n_overrides), n_base)
            unique_vals = overrides[grid_name]

        n_groups = len(unique_vals)

        if find_liquidus and not _is_parallel:
            if grouping_column == 'Pressure(System_main)':
                # P_bar is in Conditions so we can trim directly without a full divide.
                nGibbs_output = trim_superliquidus(nGibbs_output, grouping_column='Pressure(System_main)')
                # Recompute tableNo — rows were removed so the old one is stale.
                unique_vals, tableNo = np.unique(
                    nGibbs_output['Conditions']['P_bar'], return_inverse=True
                )
            else:
                # fO2 / composition grids: divide → trim each sub-table by P → re-merge.
                # This keeps add_phase_props to a single call on the already-trimmed data.
                divided = _nGibbs_models[modelname].divide_ptt_tables(nGibbs_output, tableNo)
                trimmed_list = []
                for i in range(n_groups):
                    sub = divided[f'Run {i}']
                    p_unique = sub['Conditions']['P_bar'].unique()
                    sub_gc = 'Pressure(System_main)' if len(p_unique) > 1 else None
                    trimmed_list.append(trim_superliquidus(sub, grouping_column=sub_gc))
                # Rebuild tableNo from the (now unequal) trimmed group sizes.
                sizes = [len(s['Conditions']) for s in trimmed_list]
                tableNo = np.repeat(np.arange(n_groups), sizes)
                nGibbs_output = _concat_ptt_tables(trimmed_list)

        if run_props:
            nGibbs_output = add_phase_props(nGibbs_output, Model=modelname)

        
        divide_name = grouping_column.split('(')[0]
        divided_tables = _nGibbs_models[modelname].divide_ptt_tables(nGibbs_output, tableNo)
        if divide_name == 'Pressure':
            return {
                f"P = {uniq:.1f} bars": divided_tables[f'Run {i}']
                for i, uniq in enumerate(unique_vals)
            }
        
        return {
            f"{divide_name} = {uniq:.1f}": divided_tables[f'Run {i}']
            for i, uniq in enumerate(unique_vals)
        }

    else:
        if find_liquidus and not _is_parallel:
            nGibbs_output = trim_superliquidus(nGibbs_output, grouping_column=None)

        if run_props:
            nGibbs_output = add_phase_props(nGibbs_output, Model=modelname)

        return nGibbs_output

@contextlib.contextmanager
def _suppress_stdout():
    """Redirect stdout/stderr at the OS fd level to silence C-extension output."""
    sys.stdout.flush()
    sys.stderr.flush()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            saved_out = os.dup(1)
            saved_err = os.dup(2)
            with open(os.devnull, 'w') as devnull:
                null_fd = devnull.fileno()
                os.dup2(null_fd, 1)
                os.dup2(null_fd, 2)
                try:
                    yield
                finally:
                    os.dup2(saved_out, 1)
                    os.dup2(saved_err, 2)
                    os.close(saved_out)
                    os.close(saved_err)
        except OSError:
            yield


def add_phase_props(
    Results,
    Model="MELTSv1.2.0",
    #fO2_buffer=None,
    #fO2_offset=None,
    melts=None,
    parallel=None,
    n_workers=None,
):
    """
    Augment a PTT/nGibbs output dict with per-phase MELTS thermodynamic properties.

    Wraps ``calc_phase_props_MELTS`` (serial) or ``calc_phase_props_MELTS_parallel``
    and returns the same ``Results`` dict with ``'<phase>_prop'`` DataFrames and an
    updated ``'Conditions'`` DataFrame added in-place.

    Parameters
    ----------
    Results : dict
        PTT-format dict from ``emulator.ForwardMB(..., outputs=["ptt_out"])["ptt_out"]``.
        Must contain ``'Conditions'`` (with ``'T_C'`` and ``'P_bar'``), per-phase
        composition DataFrames, and ``'mass_g'``.
    Model : str
        MELTS model string.  Accepts both PTT-style (``"MELTSv1.0.2"``) and
        nGibbs-style (``"nMELTSv1.0.2"``); the leading ``'n'`` is stripped
        automatically.  Recognised values: ``"MELTSv1.0.2"``, ``"MELTSv1.1.0"``,
        ``"MELTSv1.2.0"``, ``"pMELTS"``.
    fO2_buffer : str, optional
        fO2 buffer name (``'FMQ'``, ``'IW'``, ``'NNO'``).  Not supported in
        parallel mode.
    fO2_offset : float, optional
        Log-unit offset from the buffer.
    melts : MELTSdynamic, optional
        Pre-initialised MELTS instance for serial mode.  Created automatically
        if not supplied.
    parallel : bool
        ``True`` → ``calc_phase_props_MELTS_parallel`` (ProcessPoolExecutor).
        ``False`` (default) → ``calc_phase_props_MELTS`` (single process).
    n_workers : int or None
        Number of worker processes for parallel mode.
        ``None`` uses ``core_config.MAX_WORKERS``.

    Returns
    -------
    Results : dict
        The input dict augmented with ``'<phase>_prop'`` DataFrames and an
        updated ``'Conditions'`` DataFrame containing bulk thermodynamic columns.
    """
    # Accept both "nMELTSv1.0.2" and "MELTSv1.0.2" style strings
    melts_model = Model[1:] if Model.startswith('n') else Model

    if parallel is None:
        # Auto-enable parallel only for large batches. Callers must guard their
        # script with `if __name__ == '__main__':` + multiprocessing.freeze_support()
        # to prevent Windows spawn from re-running top-level code in each worker.
        parallel = Results['Conditions'].shape[0] > 1000

    with _suppress_stdout():
        if parallel:
            calc_phase_props_MELTS_parallel(
                Results,
                Model=melts_model,
                #fO2_buffer=fO2_buffer,
                #fO2_offset=fO2_offset,
                n_workers=n_workers,
            )
        else:
            calc_phase_props_MELTS(
                Results,
                Model=melts_model,
                #fO2_buffer=fO2_buffer,
                #fO2_offset=fO2_offset,
                melts=melts,
            )

    return Results

def nGibbsAPI(Model="nMELTSv1.0.2NoProp", comp=None, T_C=None, T_path_C=None, T_start_C=None, T_end_C=None, dt_C=None,
              T_min_C=None, T_max_C=None, T_num=None,
              P_bar=None, P_path_bar=None, P_start_bar=None, P_end_bar=None, dp_bar=None,
              P_min_bar=None, P_max_bar=None, P_num=None,
              Fe3Fet_init=None, H2O_init=None, CO2_init=None,
              Fe3Fet_Liq=None, H2O_Liq=None, CO2_Liq=None, fO2_buffer=None, fO2_offset=None,
              Suppress=None, Suppress_except=None, label=None, find_liquidus=False, verbose=False, **kwargs):

    """nGibbsAPI that receives similar input to PTT's MELTS API"""

    if comp is None:
        raise ValueError("Composition must be provided as a pandas DataFrame, pandas Series, or dictionary.")

    for arg, contents in kwargs.items():
        if contents is not None:
            if verbose:
                print(f"Warning: nGibbs in ptt does not support keyword argument '{arg}'. Ignoring {arg}={contents}.")

    modelname = Model[:-6] if Model.endswith('NoProp') else Model

    if modelname == "nMELTSv1.2.0":
        CO2_OK = True
    else:
        CO2_OK = False

    overrides = handle_kwargs(modelname, CO2_init=CO2_init, H2O_init=H2O_init, Fe3Fet_init=Fe3Fet_init,
                      Fe3Fet_Liq=Fe3Fet_Liq, H2O_Liq=H2O_Liq, CO2_Liq=CO2_Liq, fO2_buffer=fO2_buffer, fO2_offset=fO2_offset,
                      Suppress=Suppress, Suppress_except=Suppress_except, CO2_OK=CO2_OK, verbose=verbose)

    # ── find_liquidus: default T_start_C to 1800 when only range end is given ─
    if find_liquidus and T_start_C is None and T_end_C is not None and dt_C is not None:
        T_start_C = 1800

    P_vec, T_vec, is_path = build_PT_vectors(
        T_C=T_C, T_path_C=T_path_C, T_start_C=T_start_C, T_end_C=T_end_C, dt_C=dt_C,
        T_min_C=T_min_C, T_max_C=T_max_C, T_num=T_num,
        P_bar=P_bar, P_path_bar=P_path_bar, P_start_bar=P_start_bar, P_end_bar=P_end_bar, dp_bar=dp_bar,
        P_min_bar=P_min_bar, P_max_bar=P_max_bar, P_num=P_num,
    )

    return run_nGibbs(P_vec, T_vec, comp=comp, overrides=overrides, Model=Model,
                      is_path=is_path, label=label, find_liquidus=find_liquidus, verbose=verbose)

