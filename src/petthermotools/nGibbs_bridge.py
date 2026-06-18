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

import numpy as np
import pandas as pd
from petthermotools.MELTS import calc_phase_props_MELTS_parallel, calc_phase_props_MELTS


try:
    import ngibbs as _ngibbs_pkg
    from ngibbs.config.constants import ALLOWED_MELTS_OXIDES, REQUIRED_MELTS_OXIDES
    from ngibbs.utils.math_utils import grid_sample_explicit

except ImportError:
    _ngibbs_pkg = None
    ALLOWED_MELTS_OXIDES = []
    REQUIRED_MELTS_OXIDES = []
    print("nGibbs not found. Please install nGibbs.")

#BJT Prep
_FE2O3_TO_FEO = 0.8998

_nGibbs_models = {}

try:
    _nGibbs_models['nMELTSv1.2.0'] = _ngibbs_pkg.engine.API.MELTS120EmulatorCPU
    _nGibbs_models['nMELTSv1.0.2'] = _ngibbs_pkg.engine.API.MELTS102EmulatorCPU
except:
    pass

def find_oxide_cols(label):
    """
    Identify which oxide is in the string label.

    Returns the oxide name (e.g. 'SiO2', 'FeO') if found, else None.
    """

    for oxide in ALLOWED_MELTS_OXIDES:
        if oxide in label:
            return oxide
        if 'fe3fet' in label.lower():
            return 'Fe3FeT'
    return None


### FUNCTION TO GET P, T, fO2 vectos from crazy ptt input args. Make it a **kwargs type

def handle_kwargs(modelname= "nMELTSv1.2.0", **kwargs):
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
                raise ValueError(f"fO2 buffer '{val}' not recognised/supported by nGibbs. nGibbs only supports 'FMQ' and 'QFM' buffers. Other buffers will be supported in future versions")
            overrides['fO2'] = 0 #If offset is not provided, default to FMQ (0 log units from FMQ)
            continue
        if key == 'fO2_offset':
            overrides['fO2'] = val #If buffer is provided, apply offset to get logfO2 relative to FMQ
            continue
        if val is not None:
            if  "suppress" in key.lower():
                print(f"Phase suppression not supported by nGibbs, although new models may be trained for your use case. Email: benjthyer@gmail.com")
                print(f"Phases in the {modelname} isothermal Cr-bearing model:")
                print(_nGibbs_models[modelname].cr.isothermal_emulator.ml_indexer.all_phases)
                 #This can be specified further when not all emulators share the same phase. But they probably should...
            else:
                col = find_oxide_cols(key)
                if col is not None:
                    overrides[col] = val
                else:                   
                    print(f"Warning: nGibbs is ignoring Unrecognised keyword argument '{key}'")

    if 'fO2' in overrides and 'Fe3FeT' in overrides:
        raise ValueError("Simultaneous fO2 and Fe3FeT overrides are mutually exclusive!")

    # Format overrides as lists, check that no more than one is multi-valued
    multi_run = False
    for key, val in overrides.items():
        listval = condition_inputs_as_lists(val)
        if len(listval) > 1 and key != 'fO2': # fO2 is handled differently because sometimes it is needed for parallel calculation mode
            if multi_run:
                raise ValueError(f"Multiple multi-valued overrides found.\n" 
                                 "Only one multi-valued override arg (CO2, H2O, Fe3FeT liq or init) " 
                                 "is allowed to trigger a grid calculation.")
            multi_run = True
        overrides[key] = listval

    return overrides

def extract_composition_from_table(table: pd.DataFrame, overrides = {}):

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
            if oxide == 'Fe3FeT':
                if FeR is not None:
                    raise ValueError("Multiple Fe3FeT columns found in composition table.")
                else:
                    if 'Fe3FeT' in overrides:
                        FeR = overrides['Fe3FeT'][0]
                    else:
                        FeR = table.to_numpy(dtype=np.float32)[:, idx]
                    continue
            keepIDX.append(idx)
            headers.append(oxide)

    for oxide, override_val in overrides.items():
        if oxide != 'Fe3FeT':
            if oxide in headers:
                idx = headers.index(oxide)
                table.iloc[:, idx] = override_val[0]
            else:
                table[oxide] = override_val[0]
                keepIDX.append(table.columns.get_loc(oxide))
                headers.append(oxide)

    out_table = table.to_numpy(dtype=np.float32)[:, keepIDX]

    if FeR is None and 'Fe3FeT' in overrides:
        FeR = overrides['Fe3FeT'][0]

    if FeR is not None:
        if 'Fe2O3' in headers:            
            raise ValueError("Fe2O3 column already exists in composition table. Cannot apply Fe3FeT override.")
        # Recast Fe3FeT into Fe2O3 and append to the table 
        FeOt = out_table[:, headers.index('FeO')]
        Fe2O3 = FeOt * FeR / _FE2O3_TO_FEO
        out_table = np.column_stack([out_table, Fe2O3])
        out_table[:, headers.index('FeO')] = FeOt * (1.0 - FeR)
        headers.append('Fe2O3')

    assert np.all(np.isin(REQUIRED_MELTS_OXIDES, headers)), f"Missing nGibbs-required oxides in composition: {set(REQUIRED_MELTS_OXIDES) - set(headers)}"
    assert len(headers) == len(set(headers)), f"Duplicate oxide columns found in composition table: {headers}"

    return headers, out_table


def extract_composition_1D(vector: pd.Series | dict, overrides = {}):

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
            if oxide == 'Fe3FeT':
                if FeR is not None:
                    raise ValueError("Multiple Fe3FeT columns found in composition table.")
                else:
                    if 'Fe3FeT' in overrides:
                        FeR = overrides['Fe3FeT'][0]
                    else:
                        FeR = val
                    continue
            values.append(val)
            headers.append(oxide)

    for oxide, override_val in overrides.items():
        if oxide != 'Fe3FeT':
            if oxide in headers:
                idx = headers.index(oxide)
                values[idx] = override_val[0]
            else:
                values.append(override_val[0])
                headers.append(oxide)

    if FeR is None and 'Fe3FeT' in overrides:
        FeR = overrides['Fe3FeT'][0]

    if FeR is not None:
        if 'Fe2O3' in headers:            
            raise ValueError("Fe2O3 column already exists in composition table. Cannot apply Fe3FeT override.")
        # Recast Fe3FeT into Fe2O3 and append to vector
        FeOt = values[headers.index('FeO')]
        Fe2O3 = FeOt * FeR / _FE2O3_TO_FEO
        values.append(Fe2O3)
        values[headers.index('FeO')] = FeOt * (1.0 - FeR)
        headers.append('Fe2O3')

    assert np.all(np.isin(REQUIRED_MELTS_OXIDES, headers)), f"Missing nGibbs-required oxides in composition: {set(REQUIRED_MELTS_OXIDES) - set(headers)}"
    assert len(headers) == len(set(headers)), f"Duplicate oxide columns found in composition table: {headers}"

    return headers, np.array(values, dtype=np.float32).reshape(1, -1)

def build_PT_vectors(
    T_C=None, T_path_C=None, T_start_C=None, T_end_C=None, dt_C=None,
    P_bar=None, P_path_bar=None, P_start_bar=None, P_end_bar=None, dp_bar=None,
):
    """
    Resolve PTT-style T and P arguments into 1D numpy arrays (or length-1 arrays
    for scalar inputs).

    Resolution priority for each axis:
      *_path   >  *_start / *_end / d*  >  scalar (*_C or *_bar)

    Parameters
    ----------
    T_C : float or array-like
        Fixed temperature(s) in °C.
    T_path_C : array-like
        Explicit temperature path in °C (used as-is).
    T_start_C, T_end_C, dt_C : float
        Generate a temperature path with ``np.arange``.  Direction is inferred
        from start vs end; the sign of dt_C is ignored.
    P_bar : float or array-like
        Fixed pressure(s) in bar.
    P_path_bar : array-like
        Explicit pressure path in bar (used as-is).
    P_start_bar, P_end_bar, dp_bar : float
        Generate a pressure path with ``np.arange``.  Direction is inferred
        from start vs end; the sign of dp_bar is ignored.

    Returns
    -------
    P_vec, T_vec : np.ndarray
        1D float32 arrays of pressures and temperatures.
    """
    # ── Temperature ──────────────────────────────────────────────────────────
    if T_path_C is not None:
        T_vec = np.asarray(T_path_C, dtype=np.float32).ravel()
    elif T_start_C is not None and T_end_C is not None and dt_C is not None:
        if dt_C == 0:
            raise ValueError("dt_C must be non-zero.")
        step = abs(float(dt_C)) if T_end_C >= T_start_C else -abs(float(dt_C))
        T_vec = np.arange(T_start_C, T_end_C + step * 0.5, step, dtype=np.float32)
    elif T_C is not None:
        T_vec = np.atleast_1d(np.asarray(T_C, dtype=np.float32)).ravel()
    else:
        raise ValueError(
            "No temperature specified. Provide T_C, T_path_C, or T_start_C/T_end_C/dt_C."
        )

    # ── Pressure ─────────────────────────────────────────────────────────────
    if P_path_bar is not None:
        P_vec = np.asarray(P_path_bar, dtype=np.float32).ravel()
    elif P_start_bar is not None and P_end_bar is not None and dp_bar is not None:
        if dp_bar == 0:
            raise ValueError("dp_bar must be non-zero.")
        step = abs(float(dp_bar)) if P_end_bar >= P_start_bar else -abs(float(dp_bar))
        P_vec = np.arange(P_start_bar, P_end_bar + step * 0.5, step, dtype=np.float32)
    elif P_bar is not None:
        P_vec = np.atleast_1d(np.asarray(P_bar, dtype=np.float32)).ravel()
    else:
        raise ValueError(
            "No pressure specified. Provide P_bar, P_path_bar, or P_start_bar/P_end_bar/dp_bar."
        )

    return P_vec, T_vec


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


def run_nGibbs(P, T, comp, overrides = {}, Model="nMELTSv1.2.0NoProp"):
    """
    Run nGibbs for a 1D-2D PTX grid.
    1D mode: len P == 1|comp.shape[0], len T == 1|comp.shape[0], len fO2  == 1|comp.shape[0],
    2D mode: comp.shape[0] == 1. len(P) != len(T) Total sims = len(P) * len(T) * len(fO2)
    PT path mode: comp.shape[0] == 1, len(P) == len(T), len(fO2) == 1. Total sims = len(P)


    Only supported for isothermal calculations for now. Isentropic calculations will require more engineering 

    """

    if Model.endswith('NoProp'):
        modelname = Model[:-6] 
        run_props = False
    else:
        modelname = Model
        run_props = True

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
        print(f"Running nGibbs in parallel mode for {comp.shape[0]} parallel compositions")
        if grid_trigger:
            raise ValueError(f"Multi-valued override '{grid_name}' found in arguments.\n"
                             "Multi-valued compositional overrides trigger grid calculations, which are supported only for a single composition.\n"
                             "For parallel, custom calculations, place these free-to-vary values in the composition table instead")
        headers, comp_array = extract_composition_from_table(comp, overrides=overrides) # Will error if any override is multi-valued!!
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
        headers, comp_array = extract_composition_1D(comp, overrides=overrides)

        P_vec = condition_inputs_as_lists(P)
        T_vec = condition_inputs_as_lists(T)

        grid_warning_text = f" at {len(overrides[grid_name])} values for {grid_name}" if grid_trigger else ""

        if fO2 is not None:
            fO2_vec = condition_inputs_as_lists(fO2) 

            headers = ["Pressure(System_main)", "Temperature(System_main)", 'logfO2-QFM(System_main)'] + headers
        else:
            headers = ["Pressure(System_main)", "Temperature(System_main)"] + headers

        if len(P_vec) == len(T_vec): # Trigger PT path mode
            PT_arr = np.column_stack([P_vec, T_vec])
            if fO2 is None: # Simple case: one PT path, no fO2 variation
                print(f"Running nGibbs in PT path mode: {len(P_vec)} points{grid_warning_text}")
                input_array = np.hstack([PT_arr, np.tile(comp_array, (len(P_vec), 1))])
            else: # More complex: One PT path at variable fO2
                print(f"Running nGibbs in PT path mode for: {len(P_vec)} points for {len(fO2_vec)} fO2 values{grid_warning_text}")
                input_array = grid_sample_explicit([PT_arr.tolist(), fO2_vec, comp_array])
       
        
        else: # Trigger grid mode, where P and T are independent axes. If fO2 or an H2O, CO2, or Fe3FeT, override is provided, it is a third independent axis.
            if fO2 is not None:
                if (len(fO2_vec) != 1):
                    grouping_column = 'logfO2-QFM(System_main)'
                    if grid_trigger: 
                        raise ValueError(f"Attempted to do grid calculation with 3 or more dimensions: P, T, fO2, and {grid_name}. Not yet supported.")
                print(f"Running nGibbs in grid mode of {len(P_vec)} P points, {len(T_vec)} T points, and {len(fO2_vec)} fO2 values")
                input_array = grid_sample_explicit([P_vec, T_vec, fO2_vec], comp_array)
            else:
                print(f"Running nGibbs in grid mode of {len(P_vec)} P points, {len(T_vec)} T points{grid_warning_text}")
                grouping_column = 'Pressure(System_main)'
                input_array = grid_sample_explicit([P_vec, T_vec], comp_array)


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

    nGibbs_output = _nGibbs_models[modelname].ForwardMB(
            input_array, headers=headers, outputs=["ptt_out"]
        )["ptt_out"]
    
    if run_props:
        nGibbs_output = add_phase_props(nGibbs_output, Model=modelname)

    if grouping_column is not None:
                
        unique_vals, tableNo = np.unique(nGibbs_output['Conditions'][grouping_column], return_inverse=True) # Get labels and positions to sort tables

        divide_name = grouping_column.split('(')[0]
        
        divided_tables = _nGibbs_models[modelname].divide_ptt_tables(nGibbs_output, tableNo) # Recursive table divider, not model-specific

        for i, uniq in enumerate(unique_vals):
            divided_tables[f"{divide_name}={uniq:.2f}"] = divided_tables.pop(f'Run {i}') # Rename tables to reflect their organization
        
        return divided_tables

    return nGibbs_output

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

    if parallel is None and Results['Conditions'].shape[0] > 1000: # Serial is faster for small calls due to async overhead
        parallel = True
    elif parallel is None:
        parallel = False

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

def nGibbsAPI(Model = "nMELTSv1.0.2NoProp", comp = None,  T_C=None, T_path_C=None, T_start_C=None, T_end_C=None, dt_C=None,
                      P_bar=None, P_path_bar=None, P_start_bar=None, P_end_bar=None, dp_bar=None,
                      Fe3Fet_init = None, H2O_init = None, CO2_init = None, 
                      Fe3Fet_Liq = None, H2O_Liq = None, CO2_Liq = None, fO2_buffer = None, fO2_offset = None,
                      Suppress = None, Suppress_except = None, **kwargs):
    
    """nGibbsAPI that recieves similar input to PTT's MELTS API"""

    if comp is None:
        raise ValueError("Composition must be provided as a pandas DataFrame, pandas Series, or dictionary.")

    for arg, contents in kwargs.items():
        if contents is not None: #No worries about defaults the user isn't thinking about
            print(f"Warning: nGibbs in ptt does not support keyword argument '{arg}'. Ignoring {arg}={contents}.")

    overrides = handle_kwargs(Model, CO2_init = CO2_init, H2O_init = H2O_init, Fe3Fet_init = Fe3Fet_init,
                      Fe3Fet_Liq = Fe3Fet_Liq, H2O_Liq = H2O_Liq, CO2_Liq = CO2_Liq, fO2_buffer = fO2_buffer, fO2_offset = fO2_offset,
                      Suppress = Suppress, Suppress_except = Suppress_except)
    
    P_vec, T_vec = build_PT_vectors(T_C=T_C, T_path_C=T_path_C, T_start_C=T_start_C, T_end_C=T_end_C, dt_C=dt_C,
                                    P_bar=P_bar, P_path_bar=P_path_bar, P_start_bar=P_start_bar, P_end_bar=P_end_bar, dp_bar=dp_bar)
    
    return run_nGibbs(P_vec, T_vec, comp=comp, overrides=overrides, Model=Model)

