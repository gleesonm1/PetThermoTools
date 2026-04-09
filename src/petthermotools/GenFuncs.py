import numpy as np
import pandas as pd
import copy
from pathlib import Path
import psutil
import multiprocessing
import os
import re
import platform
import subprocess
# from petthermotools.Barom import *
# from petthermotools.Liq import *
# from petthermotools.Crystallise import *
# from petthermotools.MELTS import *
from petthermotools.Compositions import *
from petthermotools.core_config import MAX_WORKERS
# try:
#     from petthermotools.Holland import *
# except:
#     pass

Names = {'liquid1': '_Liq',
        'olivine1': '_Ol',
        'orthopyroxene1': '_Opx',
        'clinopyroxene1': '_Cpx',
        'garnet1': '_Grt',
        'spinel1': '_Sp',
        'alkali-feldspar1': '_Kspar',
        'quartz1': '_Qtz',
        'rhm-oxide1': '_Rhm',
        'apatite1': '_Apa',
        'olivine2': '_Ol2',
        'plagioclase1': '_Plag',
        'clinopyroxene2': '_Cpx2',
        'plagioclase2': '_Plag2',
        'spinel2': '_Sp2',
        'alkali-feldspar2': '_Kspar2',
        'garnet2': '_Grt2',
        'rhm-oxide2': '_Rhm2',
        'quartz2': '_Qtz2',
        'orthopyroxene2': '_Opx2',
        'apatite2': '_Apa2',
        'fluid1': '_Fl',
        'liquid2': '_Liq2',
        'liquid3': '_Liq3',
        'liquid4': '_Liq4',}
        # 'feldspar1': '_Fsp',
        # 'feldspar2': '_Fsp2'}

Names_MM = {'liq1': '_Liq',
            'ol1': '_Ol',
            'opx1': '_Opx',
            'cpx1': '_Cpx',
            'g1': '_Grt',
            'spl1': '_Sp',
            'fsp1': '_Fsp',
            'ol2': '_Ol2',
            'cpx2': '_Cpx2',
            'opx2': '_Opx2',
            'g2': '_Grt2',
            'fsp2': '_Fsp2',
            'spl2': '_Sp2',
            'fl1': '_Fl',
            'liq2': '_Liq2',
            'liq3': '_Liq3',
            'liq4': '_Liq4',
            'pl1': '_Plag',
            'pl2': '_Plag2',
            'afs1': '_Kspar',
            'afs2': '_Kspar2'}

Names_MM_replace = {'liq1': 'liquid1',
            'ol1': 'olivine1',
            'opx1': 'orthopyroxene1',
            'cpx1': 'clinopyroxene1',
            'g1': 'garnet1',
            'spl1': 'spinel1',
            'fsp1': 'feldspar1',
            'ol2': 'olivine1',
            'cpx2': 'clinopyroxene2',
            'opx2': 'orthopyroxene2',
            'g2': 'garnet2',
            'fsp2': 'feldspar2',
            'spl2': 'spinel2',
            'ilm1': 'rhm-oxide1',
            'ilm2': 'rhm-oxide2',
            'fl1': 'fluid1',
            'liq2': 'liquid2',
            'liq3': 'liquid3',
            'liq4': 'liquid4',
            'pl1': 'plagioclase1',
            'pl2': 'plagioclase2',
            'afs1': 'alkali-feldspar1',
            'afs2': 'alkali-feldspar2'}

# Molar masses
molar_masses = {
    "SiO2_Liq": 60.0843,
    "Al2O3_Liq": 101.9613,
    "TiO2_Liq": 79.866,
    "FeOt_Liq": 71.844,
    "MgO_Liq": 40.304,
    "MnO_Liq": 70.937,
    "CaO_Liq": 56.077,
    "Na2O_Liq": 61.979,
    "K2O_Liq": 94.196
}

# Number of cations per oxide
cation_numbers = {
    "SiO2_Liq": 1,
    "Al2O3_Liq": 2,
    "TiO2_Liq": 1,
    "FeOt_Liq": 1,
    "MgO_Liq": 1,
    "MnO_Liq": 1,
    "CaO_Liq": 1,
    "Na2O_Liq": 2,
    "K2O_Liq": 2
}

def control_path_parameters(L, P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C,
                            dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset):
    if P_bar is None:
        P_bar = [None] * int(L)
    elif type(P_bar) == float or type(P_bar) == int:
        P_bar = np.zeros(int(L)) + P_bar
    if T_C is None:
        T_C = [None] * int(L)
    elif type(T_C) == float or type(T_C) == int:
        T_C = np.zeros(int(L)) + T_C

    if P_start_bar is None:
        P_start_bar = [None] * int(L)
    elif type(P_start_bar) == float or type(P_start_bar) == int:
        P_start_bar = np.zeros(int(L)) + P_start_bar
    if T_start_C is None:
        T_start_C = [None] * int(L)
    elif type(T_start_C) == float or type(T_start_C) == int:
        T_start_C = np.zeros(int(L)) + T_start_C

    if P_end_bar is None:
        P_end_bar = [None] * int(L)
    elif type(P_end_bar) == float or type(P_end_bar) == int:
        P_end_bar = np.zeros(int(L)) + P_end_bar
    if T_end_C is None:
        T_end_C = [None] * int(L)
    elif type(T_end_C) == float or type(T_end_C) == int:
        T_end_C = np.zeros(int(L)) + T_end_C

    if dp_bar is None:
        dp_bar = [None] * int(L)
    elif type(dp_bar) == float or type(dp_bar) == int:
        dp_bar = np.zeros(int(L)) + dp_bar
    if dt_C is None:
        dt_C = [None] * int(L)
    elif type(dt_C) == float or type(dt_C) == int:
        dt_C = np.zeros(int(L)) + dt_C

    if P_path_bar is None:
        P_path_bar = [None] * int(L)
    elif len(np.shape(P_path_bar))  == 1:
        P_path_bar = np.vstack([P_path_bar] * int(L))
    if T_path_C is None:
        T_path_C = [None] * int(L)
    elif len(np.shape(T_path_C))  == 1:
        T_path_C = np.vstack([T_path_C] * int(L))

    if fO2_offset is None:
        fO2_offset = [None] * int(L)
    elif type(fO2_offset) == float or type(fO2_offset) == int:
        fO2_offset = np.zeros(int(L)) + fO2_offset

    return P_bar, T_C, P_start_bar, T_start_C, P_end_bar, T_end_C, dp_bar, dt_C, P_path_bar, T_path_C, fO2_offset

def merge_sequential_phases(output_dict):
    """
    Identifies mineral phases of the same group and merges them into a single 
    dataframe if they never coexist (mass_g > 0) in the same row.
    """
    # 1. Identify groups (e.g., 'cpx', 'spl', 'ilm')
    groups = {}
    for key in output_dict.keys():
        if key in ["sys", "Conditions"] or key.endswith("_prop"):
            continue
        
        match = re.match(r"([a-zA-Z\-]+)([0-9]+)", key)
        if match:
            prefix, _ = match.groups()
            if prefix not in groups:
                groups[prefix] = []
            groups[prefix].append(key)

    new_dict = {k: v for k, v in output_dict.items() if k in ["sys", "Conditions"]}
    
    for prefix, instances in groups.items():
        # Sort instances to ensure deterministic merging (cpx1, cpx2...)
        instances.sort(key=lambda x: int(re.search(r'\d+', x).group()))
        
        # We will keep a list of "active" merged bins for this mineral
        # Each bin: {'data': df, 'prop': df, 'indices_with_mass': set}
        bins = []

        for inst in instances:
            df_data = output_dict[inst]
            df_prop = output_dict[f"{inst}_prop"]
            
            # Identify row indices where this phase exists (mass > 0)
            current_active_rows = set(df_prop.index[df_prop['mass_g'] > 0].tolist())
            
            merged = False
            for b in bins:
                # Check for overlap: Is the intersection of active rows empty?
                if not current_active_rows.intersection(b['active_rows']):
                    # NO OVERLAP: Merge this instance into the existing bin
                    # We combine data; assuming identical columns, we fill zeros/NaNs
                    b['data'] = b['data'].add(df_data, fill_value=0)
                    b['prop'] = b['prop'].add(df_prop, fill_value=0)
                    b['active_rows'].update(current_active_rows)
                    merged = True
                    break
            
            if not merged:
                # OVERLAP or FIRST INSTANCE: Create a new bin
                bins.append({
                    'data': df_data.copy(),
                    'prop': df_prop.copy(),
                    'active_rows': current_active_rows
                })

        # 2. Re-insert merged/separated bins back into the dictionary with new numbering
        for i, b in enumerate(bins, 1):
            new_label = f"{prefix}{i}"
            new_dict[new_label] = b['data']
            new_dict[f"{new_label}_prop"] = b['prop']

    return new_dict

def estimate_t(comp=None, P_bar=None):
    # Use eq22 from Putirka 2008 to estimate the liquidus temperate of a melt - provides a starting point for the liquidus search routine
    if isinstance(comp, dict):
        df = pd.DataFrame([comp])
        comp_is_single = True
    elif isinstance(comp, pd.DataFrame):
        df = comp.copy()
        comp_is_single = False
    else:
        raise TypeError("comp must be dict or pandas DataFrame")

    df = df.fillna(0.0)

    n_rows = len(df)

    if np.isscalar(P_bar):
        P = pd.Series(P_bar, index=df.index)
    else:
        P = pd.Series(P_bar)

        if comp_is_single:
            # Single composition, multiple pressures
            df = pd.concat([df]*len(P), ignore_index=True)
            P = pd.Series(P.values, index=df.index)
            n_rows = len(df)
        else:
            if len(P) != n_rows:
                raise ValueError("Length of P_bar must match number of compositions")

            P.index = df.index

    P_GPa = P / 10000.0

    oxide_moles = pd.DataFrame(index=df.index)

    for oxide in molar_masses:
        oxide_moles[oxide] = df.get(oxide, 0.0) / molar_masses[oxide]

    cation_moles = pd.DataFrame(index=df.index)

    for oxide in oxide_moles.columns:
        moles = oxide_moles[oxide] * cation_numbers[oxide]
        cation_moles[oxide] = moles

    cation_moles = cation_moles.div(cation_moles.sum(axis=1), axis=0)

    # Extract components
    Fe = cation_moles.get('FeOt_Liq', 0)
    Mg = cation_moles.get('MgO_Liq', 0)
    Mn = cation_moles.get('MnO_Liq', 0)
    Ca = cation_moles.get('CaO_Liq', 0)
    Si = cation_moles.get('SiO2_Liq', 1e-12)
    Al = cation_moles.get('Al2O3_Liq', 0)
    Ti = cation_moles.get('TiO2_Liq', 0)

    CNM = Fe + Mg + Mn + Ca
    CSi = Si

    # Numerical safety
    CNM = np.clip(CNM, 1e-12, None)
    CSi = np.clip(CSi, 1e-12, None)
    Al  = np.clip(Al, 0, 0.999999)
    Ti  = np.clip(Ti, 0, 0.999999)
    Mn  = np.clip(Mn, 0, 0.999999)
    Fe  = np.clip(Fe, 0, 0.999999)
    Mg  = np.clip(Mg, 1e-12, None)

    NF = (7/2)*np.log(1-Al) + 7*np.log(1-Ti)

    H2O = df.get("H2O_Liq", 0)

    lnD = np.log(
        (0.666 - (-0.049*Mn + 0.027*Fe))
        /
        (Mg + 0.259*Mn + 0.299*Fe)
    )

    T = (
        (15294.6 + 1318.8*P_GPa + 2.4834*(P_GPa**2))
        /
        (8.048
         + 2.8352*lnD
         + 2.097*np.log(1.5*CNM)
         + 2.575*np.log(3*CSi)
         - 1.41*NF
         + 0.222*H2O
         + 0.5*P_GPa)
    )

    T[T < 750] = 750

    if comp_is_single and np.isscalar(P_bar):
        return T.iloc[0]

    return T

def memory_limit(cores = multiprocessing.cpu_count()):
    if psutil.virtual_memory().available/(1034**3) < 2:
        if cores > 2:
            print("Limiting calculations to 2 processes due to limitations in available memory")
            return 2
        else:
            return cores
    elif psutil.virtual_memory().available/(1034**3) < 4:
        if cores > 4:
            print("Limiting calculations to 4 processes due to limitations in available memory")
            return 4
        else:
            return cores
    elif psutil.virtual_memory().available/(1024**3)<6:
        if cores > 6:
            print("Limiting calculations to 6 processes due to limitations in available memory")
            return 6
        else:
            return cores
    else:
        return cores

def check_array(var):
    '''
    Checks the provided data is a numpy array. If a list is provided instead it is converted.
    '''
    if type(var) == list:
        var_new = np.array(var)
        return var_new
    else:
        return var

def rename_keys_with_prefix(d, mapping):
    '''
    Renames keys in a dictionary based on an exact match or a prefix match defined 
    in a mapping dictionary.

    This is primarily used to convert short-form phase names (e.g., 'ol1') from 
    MAGEMinCalc into longer names used by alphaMELTS for Python (e.g., 'olivine1').

    Parameters:
    ----------
    d : dict
        The input dictionary whose keys need renaming.
    mapping : dict
        A mapping dictionary where keys are the old names/prefixes (e.g., 'ol1') 
        and values are the new names/prefixes (e.g., 'olivine1').

    Returns:
    ---------
    new_dict : dict
        A new dictionary with the keys renamed according to the mapping logic.
    '''
    new_dict = {}
    for k, v in d.items():
        # 1. Exact match
        if k in mapping:
            new_key = mapping[k]

        # 2. Prefix match (e.g. ol1_ → olivine1_)
        else:
            new_key = k
            for old, new in mapping.items():
                if k.startswith(old + "_"):
                    new_key = k.replace(old, new, 1)
                    break

        new_dict[new_key] = v
    return new_dict

# Create a global flag to ensure we only setup once per session
_JULIA_INITIALIZED = False

def _ensure_julia_ready():
    """
    Private helper to configure Julia only when actually needed.
    This prevents 'import PetThermoTools' from being slow.
    """
    global _JULIA_INITIALIZED
    
    if not _JULIA_INITIALIZED:
        import os
        import juliapkg
        from pathlib import Path
        
        env_path = Path.home() / ".petthermotools_julia_env"
        
        if env_path.exists():
            # Register the project and version requirement
            os.environ["JULIA_PROJECT"] = str(env_path)
            juliapkg.require_julia("~1.11")
            _JULIA_INITIALIZED = True
            print("Julia environment detected.")
        else:
            print("Julia environment not available. Please run the `install_MAGEMinCalc()` function first.")

def get_performance_core_count():
    system = platform.system()
    arch = platform.machine() # 'arm64' or 'x86_64'
    
    # --- macOS Logic ---
    if system == "Darwin":
        if arch == "arm64":
            try:
                # Apple Silicon: perflevel0 is Performance
                cmd = ["sysctl", "-n", "hw.perflevel0.physicalcpu"]
                return int(subprocess.check_output(cmd).decode().strip())
            except Exception:
                return psutil.cpu_count(logical=False) or 1
        else:
            # Intel Mac: No E-cores, so all physical cores are P-cores
            return psutil.cpu_count(logical=False) or 1

    # --- Windows & Linux Logic (Intel Hybrid/AMD) ---
    lt = psutil.cpu_count(logical=True) or 1
    lf = psutil.cpu_count(logical=False) or 1

    if lt > lf:
        # Heuristic: P-cores have Hyper-threading, E-cores don't.
        # Math: P = logical - physical
        p_cores = lt - lf
        # If calculation seems wrong (e.g., E-cores > P-cores), 
        # fallback to physical core count.
        return p_cores if p_cores > 0 else lf
    
    # Fallback for non-hyperthreaded or non-hybrid CPUs
    return lf

def activate_petthermotools_env():
    '''
    Activates the custom Julia environment (.petthermotools_julia_env) required 
    for running MAGEMinCalc calculations via the JuliaCall interface.

    Parameters:
    ----------
    None

    Returns:
    ----------
    None. Activates the environment using `Pkg.activate`.
    '''
    _ensure_julia_ready()
    from juliacall import Main as jl
    env_dir = Path.home() / ".petthermotools_julia_env"
    jl_env_path = env_dir.as_posix()
    jl.seval(f"""import Pkg
             Pkg.activate(expanduser("{jl_env_path}"))
             Pkg.precompile()""")
    jl.seval("using Distributed")
    jl.Distributed.addprocs(memory_limit(cores = get_performance_core_count()))
    jl.seval("@everywhere using MAGEMinCalc")
    jl.seval(f"""
        if pkgversion(MAGEMinCalc) != v"0.6.1"
            error("Incorrect version! Expected v0.6.1. Please run ptt.update_MAGEMinCalc()")
        else
             println("Julia Environment Ready. MAGEMinCalc v0.6.1 detected")
        end"""
    )
    

def to_float(x):
    '''
    Converts input data (scalar, list, tuple, or numpy array) to float type.

    Parameters:
    ----------
    x : scalar (int, float) or array-like (list, tuple, np.ndarray)
        The input data to be converted. Can also be None.

    Returns:
    ----------
    float or list of float or numpy.ndarray of float or None
        The input data converted to float type. Returns None if input is None.
    '''
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, (list, tuple)):
        return [float(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.astype(float)
    return x  # leave unchanged if unexpected type

def label_results(Result,label):
    '''
    Renames the keys of a results dictionary based on a specified input parameter 
    (e.g., pressure, H2O content) to create descriptive, sortable labels.

    Example: Renames 'Run 1' to 'P = 1000 bars' if label='P_bar'.

    Parameters:
    ----------
    Result : dict
        A dictionary where keys are generic (e.g., 'Run 0') and values are the 
        standard simulation output dictionaries (which contain an 'Input' key).
    label : str
        The input parameter to use for labeling. Common options include:
        'pressure', 'P_bar', 'CO2', 'CO2_init', 'fO2', 'fO2_offset', 'H2O', 'H2O_init'.

    Returns:
    ----------
    new_out : dict
        The result dictionary with keys replaced by descriptive, sorted strings.
        If the label is not found, a copy of the original dictionary is returned.
    '''
    Results = Result.copy()
    new_out = {}
    if  label == "CO2" or label == "CO2_init":
        for r in Results:
            new_out['CO2 = ' + str(Results[r]['Input']['comp']['CO2_Liq']) + ' wt%'] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[1])))
    elif label == "pressure" or label == "P" or label == "P_bar":
        for r in Results:
            new_out['P = ' + str(Results[r]['Input']['P_bar']) + ' bars'] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[1])))
    elif label == "P_start_bar":
        for r in Results:
            new_out['P_start = ' + str(Results[r]['Input']['P_start_bar']) + ' bars'] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[1])))
    elif label == "fO2" or label == "fO2_offset":
        for r in Results:
            new_out['fO2 = ' + Results[r]['Input']['fO2_buffer'] + ' ' + str(round(Results[r]['Input']['fO2_offset'],2))] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[2])))
    elif label == 'H2O' or label == "H2O_init":
        for r in Results:
            new_out['H2O = ' + str(Results[r]['Input']['comp']['H2O_Liq']) + ' wt%'] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[1])))
    elif label == 'Fe3Fet' or label == "Fe3Fet_init":
        for r in Results:
            new_out['Fe3Fet = ' + str(Results[r]['Input']['comp']['Fe3Fet_Liq'])] = Results[r].copy()
        new_out = dict(sorted(new_out.items(), key=lambda x: float(x[0].split('=')[1].split(' ')[1])))
    
    if len(new_out) == 0:
        new_out = Results.copy()
    
    return new_out

def comp_check(comp_lith, Model, MELTS_filter, Fe3Fet):
    '''
    Checks and preprocesses the input bulk composition before full standardization.

    Handles string inputs (e.g., sample names) by looking them up in a global 
    Compositions dictionary (if not pyMelt), converts Series to dicts, and ensures 
    minimum required volatile/trace components are present for MELTS.

    Parameters:
    ----------
    comp_lith : str, dict, or pd.Series
        The input composition, oxide concentration in wt% or a string name corresponding to one of the saved compositions (G2, KLB-1, KG1).
    Model : str
        Thermodynamic model to be used (e.g., "MELTSv1.0.2", "Green2025").
    MELTS_filter : bool
        If True and using a MELTS model, ensures K2O, P2O5, H2O, and CO2 are 
        initialized to zero if missing.
    Fe3Fet : float, optional
        Initial Fe³⁺/Fe_total ratio override.

    Returns:
    ----------
    comp : dict or pd.DataFrame
        The preprocessed composition dictionary or DataFrame.
    '''
    if type(comp_lith) == str:
        if Model != "pyMelt":
            comp = Compositions[comp_lith]
        else:
            comp = comp_lith
    else:
        comp = comp_lith.copy()

    # if comp is entered as a pandas series, it must first be converted to a dict
    if Model != "pyMelt":
        if type(comp) == pd.core.series.Series:
            comp = comp.to_dict()

        comp = comp_fix(Model = Model, comp = comp, Fe3Fet_Liq = Fe3Fet)

    if "MELTS" in Model and MELTS_filter == True:
        if type(comp) == pd.core.frame.DataFrame:
            comp['K2O_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['P2O5_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['H2O_Liq'] = np.zeros(len(comp['SiO2_Liq']))
            comp['CO2_Liq'] = np.zeros(len(comp['SiO2_Liq']))
        else:
            comp['K2O_Liq'] = 0
            comp['P2O5_Liq'] = 0
            comp['H2O_Liq'] = 0
            comp['CO2_Liq'] = 0
    
    return comp

def comp_fix(Model=None, comp=None, Fe3Fet_Liq=None, H2O_Liq=None, CO2_Liq=None, keep_columns=False):
    '''
    Standardizes and ensures correct formatting of input compositions for thermodynamic calculations.

    This function performs several critical operations:
    1. **Fe Speciation:** Calculates FeOt and Fe3Fet columns if FeO/Fe2O3 are present.
    2. **Header Standardization:** Renames common headers (e.g., 'SiO2') to 
       standard "ideal" headers (e.g., 'SiO2_Liq') that are used in subsequent functions.
    3. **Initialization:** Initializes any missing required oxide/volatile columns to 0.0.
    4. **Redox/Volatile Overrides:** Applies optional external `Fe3Fet_Liq`, `H2O_Liq`, 
       and `CO2_Liq` overrides, handling scalar and array-like inputs.
    5. **Normalization:** Normalizes the non-volatile oxides to 100% *minus* the 
       specified volatile contents (H2O_Liq, CO2_Liq).
    6. **Type Enforcement:** Ensures all required component values are float type.

    Parameters:
    ----------
    Model : str, default "MELTSv1.0.2"
        Thermodynamic model (e.g., "MELTSvx.x.x" or "Green2025").
    comp : dict or pd.DataFrame
        Input composition, oxide values in wt%.
    Fe3Fet_Liq, H2O_Liq, CO2_Liq : float or np.ndarray, optional
        Optional overrides for redox or volatiles. If arrays are used, the output
        will be a DataFrame with the corresponding number of rows.
    keep_columns : bool, default False
        If False, returns only the standardized columns required for calculations 
        (Columns_ideal). If True, preserves all original columns.

    Returns:
    ----------
    comp : dict or pd.DataFrame
        New composition with correct headers, types, and normalization.
    '''
    if Model is None:
        Model = "MELTSv1.0.2"

    if isinstance(comp, pd.DataFrame):
        comp = comp.copy(deep=True)
    elif isinstance(comp, dict):
        comp = copy.deepcopy(comp)

    # Create a copy to avoid SettingWithCopy warnings
    Comp_start = comp.copy()

    # --- 1. PRE-CALCULATE Fe SPECIATION IF AVAILABLE ---
    if isinstance(comp, pd.DataFrame):
        # Helper to ensure existing Fe cols are float before math
        fe_cols = [c for c in ['FeO_Liq', 'Fe2O3_Liq', 'FeO', 'Fe2O3'] if c in comp.columns]
        comp[fe_cols] = comp[fe_cols].astype(float)

        if "FeO_Liq" in comp.columns and "Fe2O3_Liq" in comp.columns:
            if "FeOt_Liq" not in comp.columns:
                comp['FeOt_Liq'] = comp['FeO_Liq'] + 0.8998 * comp['Fe2O3_Liq']
            if "Fe3Fet_Liq" not in comp.columns:
                comp['Fe3Fet_Liq'] = (1 - comp['FeO_Liq'] / (comp['FeO_Liq'] + 0.8998 * comp['Fe2O3_Liq']))
        elif "FeOt_Liq" not in comp.columns and "FeO_Liq" in comp.columns and "Fe2O3_Liq" not in comp.columns:
            comp['FeOt_Liq'] = comp["FeO_Liq"]
        elif "FeOt_Liq" not in comp.columns and "FeO_Liq" not in comp.columns and "Fe2O3_Liq" in comp.columns:
            comp['FeOt_Liq'] = 0.8998*comp['Fe2O3_Liq']
        
        if "FeO" in comp.columns and "Fe2O3" in comp.columns:
            if "FeOt" not in comp.columns:
                comp['FeOt'] = comp['FeO'] + 0.8998 * comp['Fe2O3']
            if "Fe3Fet" not in comp.columns:
                comp['Fe3Fet'] = 1 - comp['FeO'] / (comp['FeO'] + 0.8998 * comp['Fe2O3'])
        elif "FeOt" not in comp.columns and "FeOt_Liq" not in comp.columns and "FeO" in comp.columns and "Fe2O3" not in comp.columns:
            comp['FeOt'] = comp["FeO"]
        elif "FeOt" not in comp.columns and "FeOt_Liq" not in comp.columns and "FeO" not in comp.columns and "Fe2O3" in comp.columns:
            comp['FeOt'] = 0.8998*comp['Fe2O3']

        Comp_start = comp.copy()
    
    elif isinstance(comp, dict):
        for k in ['FeO_Liq', 'Fe2O3_Liq', 'FeO', 'Fe2O3']:
            if k in comp: comp[k] = float(comp[k])

        if "FeO_Liq" in comp and "Fe2O3_Liq" in comp:
            if "FeOt_Liq" not in comp:
                comp['FeOt_Liq'] = comp['FeO_Liq'] + 0.8998 * comp['Fe2O3_Liq']
            if "Fe3Fet_Liq" not in comp:
                comp['Fe3Fet_Liq'] = (1 - comp['FeO_Liq'] / (comp['FeO_Liq'] + 0.8998 * comp['Fe2O3_Liq']))
        elif "FeOt_Liq" not in comp and "FeO_Liq" in comp and "Fe2O3_Liq" not in comp:
            comp['FeOt_Liq'] = comp["FeO_Liq"]
        elif "FeOt_Liq" not in comp and "FeO_Liq" not in comp and "Fe2O3_Liq" in comp:
            comp['FeOt_Liq'] = 0.8998*comp['Fe2O3_Liq']
                
        if "FeO" in comp and "Fe2O3" in comp:
            if "FeOt" not in comp:
                comp['FeOt'] = comp['FeO'] + 0.8998 * comp['Fe2O3']
            if "Fe3Fet" not in comp:
                comp['Fe3Fet'] = 1 - comp['FeO'] / (comp['FeO'] + 0.8998 * comp['Fe2O3'])
        elif "FeOt_Liq" not in comp and "FeOt" not in comp and "FeO" in comp and "Fe2O3" not in comp:
            comp['FeOt'] = comp["FeO"]
        elif "FeOt_Liq" not in comp and "FeOt" not in comp and "FeO" not in comp and "Fe2O3" in comp:
            comp['FeOt'] = 0.8998*comp['Fe2O3']
        
        Comp_start = comp.copy()

    # --- 2. DEFINE COLUMNS BASED ON MODEL ---
    if "MELTS" in Model:
        Columns_bad = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeOt', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'CO2', 'Fe3Fet']
        Columns_ideal = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'Cr2O3_Liq', 'FeOt_Liq', 'MnO_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'P2O5_Liq', 'H2O_Liq', 'CO2_Liq', 'Fe3Fet_Liq']
    else:
        Columns_bad = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O', 'H2O', 'Cr2O3', 'Fe3Fet']
        Columns_ideal = ['SiO2_Liq', 'TiO2_Liq', 'Al2O3_Liq', 'FeOt_Liq', 'MgO_Liq', 'CaO_Liq', 'Na2O_Liq', 'K2O_Liq', 'Cr2O3_Liq', 'H2O_Liq', 'Fe3Fet_Liq']

    # --- 3. STANDARDIZE HEADERS AND ENFORCE FLOAT TYPE ---
    if isinstance(comp, pd.DataFrame):
        # Rename bad columns
        rename_dict = {el: el + '_Liq' for el in Comp_start.columns if el in Columns_bad}
        comp = comp.rename(columns=rename_dict)

        # Initialize missing IDEAL columns with 0.0
        for el in Columns_ideal:
            if el not in comp.columns:
                comp[el] = 0.0
        
        # Enforce float ONLY on the ideal columns (preserve strings in other columns)
        cols_to_numeric = [col for col in comp.columns if col in Columns_ideal]
        comp[cols_to_numeric] = comp[cols_to_numeric].astype(float)

    elif isinstance(comp, dict):
        # Move bad keys to good keys
        for el in list(Comp_start.keys()):
            if el in Columns_bad:
                comp[el + '_Liq'] = float(comp[el])
                del comp[el]

        # Initialize missing keys
        for el in Columns_ideal:
            if el not in comp:
                comp[el] = 0.0
        
        # Enforce float ONLY on ideal keys
        for key in comp:
            if key in Columns_ideal:
                comp[key] = float(comp[key])

    # --- 4. HANDLE Fe3Fet_Liq ---
    if Fe3Fet_Liq is not None:
        if isinstance(comp, dict):
            if not isinstance(Fe3Fet_Liq, np.ndarray):
                comp['Fe3Fet_Liq'] = float(Fe3Fet_Liq)
            else:
                Comp = pd.DataFrame([comp] * len(Fe3Fet_Liq))
                Comp['Fe3Fet_Liq'] = Fe3Fet_Liq.astype(float)
                comp = Comp.copy()
        else:
            comp['Fe3Fet_Liq'] = np.full(len(comp), Fe3Fet_Liq, dtype=float)

    # --- 5. HANDLE VOLATILES AND NORMALIZATION ---
    excluded_from_sum = ["Fe3Fet_Liq"]

    if H2O_Liq is None:
        H2O = False
    else:
        H2O = True
        excluded_from_sum = excluded_from_sum + ["H2O_Liq"]

    if CO2_Liq is None:
        CO2 = False
    else:
        CO2 = True
        excluded_from_sum = excluded_from_sum + ["CO2_Liq"]

    if H2O_Liq is None and CO2_Liq is not None:
        H2O_Liq = np.zeros(len(CO2_Liq)) if isinstance(CO2_Liq, np.ndarray) else 0.0
    elif H2O_Liq is not None and CO2_Liq is None:
        CO2_Liq = np.zeros(len(H2O_Liq)) if isinstance(H2O_Liq, np.ndarray) else 0.0
    
    if isinstance(H2O_Liq, np.ndarray) and not isinstance(CO2_Liq, np.ndarray):
        CO2_Liq = np.full(len(H2O_Liq), CO2_Liq, dtype=float)
    elif isinstance(CO2_Liq, np.ndarray) and not isinstance(H2O_Liq, np.ndarray):
        H2O_Liq = np.full(len(CO2_Liq), H2O_Liq, dtype=float)

    
    
    # Case A: No external volatiles
    if H2O_Liq is None and CO2_Liq is None:
        if isinstance(comp, dict):
            total = sum(v for k, v in comp.items() if k in Columns_ideal and k != "Fe3Fet_Liq")
            for n in comp:
                if n in Columns_ideal and n != "Fe3Fet_Liq":
                    comp[n] = (comp[n] / total) * 100.0
        else:
            oxides = [c for c in comp.columns if c in Columns_ideal and c != "Fe3Fet_Liq"]
            total = comp[oxides].sum(axis=1)
            comp.loc[:, oxides] = comp[oxides].div(total, axis=0) * 100.0

    # Case B: External volatiles provided
    else:
        volatiles = H2O_Liq + CO2_Liq

        if isinstance(comp, dict):
            if isinstance(H2O_Liq, np.ndarray) or isinstance(CO2_Liq, np.ndarray):
                length = len(H2O_Liq) if isinstance(H2O_Liq, np.ndarray) else len(CO2_Liq)
                comp = pd.DataFrame([comp] * length)
                if H2O:
                    comp['H2O_Liq'] = H2O_Liq
                if CO2:
                    comp['CO2_Liq'] = CO2_Liq
            else:
                if H2O:
                    comp['H2O_Liq'] = float(H2O_Liq)
                if CO2:
                    comp['CO2_Liq'] = float(CO2_Liq)
        else:
            if H2O:
                comp['H2O_Liq'] = H2O_Liq if isinstance(H2O_Liq, np.ndarray) else float(H2O_Liq)
            if CO2:
                comp['CO2_Liq'] = CO2_Liq if isinstance(CO2_Liq, np.ndarray) else float(CO2_Liq)

        if isinstance(comp, dict):
            total = sum(v for k, v in comp.items() if k in Columns_ideal and k not in excluded_from_sum)
            for n in comp:
                if n in Columns_ideal and n not in excluded_from_sum:
                    comp[n] = (comp[n] / total) * (100.0 - volatiles)
        else:
            oxides = [c for c in comp.columns if c in Columns_ideal and c not in excluded_from_sum]
            total = comp[oxides].sum(axis=1)
            scale = 100.0 - volatiles
            comp.loc[:, oxides] = comp[oxides].div(total, axis=0).mul(scale, axis=0)

    # --- 6. FINAL FILTERING BASED ON KEEP_COLUMNS ---
    if not keep_columns:
        if isinstance(comp, pd.DataFrame):
            # Select only columns present in Columns_ideal
            cols_to_return = [c for c in comp.columns if c in Columns_ideal]
            comp = comp[cols_to_return]
        elif isinstance(comp, dict):
            # Remove keys not in Columns_ideal
            keys_to_remove = [k for k in comp.keys() if k not in Columns_ideal]
            for k in keys_to_remove:
                del comp[k]

    return comp

def compile_results(results):
    # identify all unique keys to set up output dataframes
    unique_keys = set()
    properties = {}

    for step_data in results.values():
        for key, value in step_data.items():
            if "_keys" not in key:
                unique_keys.add(key)
            else:
                base_key = key[:-5]
                if base_key not in properties:
                    properties[base_key] = value

    unique_keys = sorted(unique_keys)

    sorted_steps = sorted(results.keys())
    n_steps = len(sorted_steps)

    index_map = {step: i for i, step in enumerate(sorted_steps)}

    # create new output to populate with results
    new_results = {}

    for key in unique_keys:
        if key == "Conditions":
            columns = ['temperature', 'pressure', 'h', 's', 'v', 'mass', 'dvdp', 'logfO2']
            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

        elif "_prop" not in key:
            columns = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'Cr2O3', 'FeO',
                    'MnO', 'MgO', 'CaO', 'Na2O', 'K2O',
                    'P2O5', 'H2O', 'CO2']
            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

        elif "_key" not in key:
            columns = properties[key]
            new_results[key] = pd.DataFrame(index=range(n_steps), columns=columns, dtype=float)

    # fill output
    for step, step_data in results.items():

        row_i = index_map[step]

        for key, value in step_data.items():

            if key.endswith("_keys"):
                continue

            if key not in new_results:
                continue
            
            if len(value) == len(new_results[key].iloc[row_i]):
                new_results[key].iloc[row_i] = value
            else:
                print(f"Size mismatch, potentially errored results on step {row_i}")
    
    return new_results

def stich(Res, multi = None, Model = None, Frac_fluid = None, Frac_solid = None):
    '''
    Processes, cleans, and combines the raw output from single or multiple 
    alphaMELTS/MAGEMinCalc calculations into a standardized format.

    Key operations performed by this function and its worker (`stich_work`):
    1. **Unit Conversion/Renaming:** Converts units (e.g., P_kbar to P_bar, rho to kg/m³), 
       renames property columns (e.g., `h` to `h_J`), and adds density and volume columns.
    2. **Fe Speciation:** Calculates FeOt and Fe3Fet for all phases in MELTS outputs.
    3. **Cleanup:** Removes invalid steps (e.g., where mass/h/T is 0.0).
    4. **Key Standardization:** Renames short-form phase keys (e.g., 'liq1') to 
       long-form (e.g., 'liquid1').
    5. **Suffix Application:** Adds standard suffixes (e.g., `_Liq`, `_Ol`) to the 
       composition and property columns of each phase DataFrame.
    6. **Aggregation:** Creates three new key DataFrames: 'All' (combined conditions and 
       all phase properties), 'mass_g', 'volume_cm3', and 'rho_kg/m3' (phase mass/volume/density 
       arrays, including cumulative values if requested).

    Parameters:
    ----------
    Res : dict
        Final results from the crystallisation/decompression calculations.
    multi : bool, optional
        If True, `Res` is a dictionary of dictionaries (from a multi-process run).
        If None/False, `Res` is a single calculation output dictionary.
    Model : str
        e.g., "MELTSvx.x.x" or "Green2025" (for MAGEMinCalc) to determine the correct processing logic.
    Frac_fluid : bool, optional
        If True, calculates cumulative mass for the fluid phase(s).
    Frac_solid : bool, optional
        If True, calculates cumulative mass for the solid phase(s).

    Returns:
    ----------
    Results : dict
        A copy of the input dictionary (`Res`) with standardized column names, units, 
        and the aggregated 'All', 'mass_g', 'volume_cm3', and 'rho_kg/m3' DataFrames.
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
        if Model == "Weller2024":
            Order = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O', 'Fe3Fet']
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
    Internal worker function for `stich`. Performs the core unit conversions, 
    data cleaning, Fe speciation, and property standardization on a single 
    simulation result dictionary.

    Parameters:
    ----------
    Results : dict
        The raw output dictionary from a single MELTS/MAGEMinCalc run.
    Order : list
        List of oxide names defining the desired output column order.
    Model : str, default "MELTS"
        Thermodynamic model type ("MELTS" or "Green2025").
    Frac_fluid, Frac_solid : bool, optional
        Flags for cumulative mass calculation.

    Returns:
    ----------
    Results : dict
        The processed and cleaned result dictionary.
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
        if "MELTS" in Model:
            Results[R].loc[:,'FeOt'] = Results[R].loc[:,'FeO'].fillna(0.0) + 71.844/(159.69/2)*Results[R].loc[:,'Fe2O3'].fillna(0.0)
            Results[R].loc[:,'Fe3Fet'] = (71.844/(159.69/2)*Results[R].loc[:,'Fe2O3'].fillna(0.0))/Results[R].loc[:,'FeOt'].fillna(0.0)
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
            Results[R][Results[R + '_prop']['mass_g'] == 0.0] = np.nan
            Results[R] = Results[R][Order]

    if "MELTS" in Model:
        Remove = np.where(Results['Conditions']['h'] == 0.0)[0]
    else:
        Remove = np.where(Results['Conditions']['T_C'] == 0.0)[0]

    for R in Results:
        if len(Remove) > 0:
            Results[R] = Results[R].drop(labels = Remove)
            
    ### Fix units
    if "MELTS" in Model:
        Results['Conditions'] = Results['Conditions'].rename(columns=dict(zip(['h', 's', 'v', 'mass', 'logfO2', 'dvdp'],['H_J', 'S_J/K', 'V_cm^3', 'mass_g', 'log10(fO2)', 'dVdP_cm^3/bar'])))
        Results['Conditions']['rho_kg/m^3'] = 1000*Results['Conditions']['mass_g']/Results['Conditions']['V_cm^3']
        for n in SN:
            Results[n + '_prop']['rho'] = Results[n + '_prop']['rho']*1000 #convert to kg/m3
            Results[n + '_prop']['cp'] = Results[n + '_prop']['cp']/(Results[n + '_prop']['mass']/1000) #
            Results[n + '_prop'] = Results[n + '_prop'].rename(columns=dict(zip(['g','h','s','v','cp','dcpdt','dvdt','dpdt',
                                                                                 'd2vdt2','d2vdtdp','d2vdp2','rho','mass'],
                                                                                 ['G_J','H_J','S_J/K', 'V_cm^3','Cp_J/(kg.K^2)',
                                                                                 'dCpdT_J/(kg.K^2)','dVdT_cm^3/K','dPdT_bar/K',
                                                                                 'd2VdT2_cm^3/K^2', 'd2VdTdP_cm^3/(bar.K)', 'd2VdP2_cm^3/bar^2',
                                                                                 'rho_kg/m^3','mass_g'])))
                                                                                # ['g_J','h_J','s_J/K', 'v_cm3','cp_J/kg/K',
                                                                                #  'dcpdt_J/K','dvdt_cm3/K','dpdt_bar/K',
                                                                                #  'd2vdt2_cm3/K2', 'd2vdtdp_cm3/bar.K', 'd2vdp2_cm3/bar2',
                                                                                #  'rho_kg/m3','mass_g'])))
    else:
        Results['Conditions']['h_kJ/mol'] = 1000*Results['Conditions']['h_kJ/mol']*Results['Conditions']['mass_g']/Results['Conditions']['mass_per_mole_g/mol']
        Results['Conditions']['s_kJ/K'] = 1000*Results['Conditions']['s_kJ/K']*Results['Conditions']['mass_g']/Results['Conditions']['mass_per_mole_g/mol']
        Results['Conditions']['v_cm3'] = Results['Conditions']['mass_g']/(Results['Conditions']['rho_kg/m3']/1000)
        Results['Conditions'] = Results['Conditions'].rename(columns=dict(zip(['h_kJ/mol','s_kJ/K', 'v_cm3',
                                                                               'rho_kg/m3', 'cp_J/K/kg'],
                                                                              ['H_J','S_J/K','V_cm^3',
                                                                               'rho_kg/m^3', 'Cp_J/(kg.K^2)'])))

        for n in SN:
            Results[n + '_prop']['v_cm3'] = Results['Conditions']['V_cm^3']*Results[n+'_prop']['vol%']
            Results[n+'_prop']['h_kJ/mol'] = 1000*Results[n+'_prop']['h_kJ/mol']*Results[n+'_prop']['mol%']*(Results[n+'_prop']['mass_g']/(Results['Conditions']['mass_per_mole_g/mol']*Results[n+'_prop']['mass%']))
            Results[n+'_prop']['s_kJ/K'] = 1000*Results[n+'_prop']['s_kJ/K']*(Results[n+'_prop']['mass_g']/(Results['Conditions']['mass_per_mole_g/mol']*Results[n+'_prop']['mass%']))
            Results[n+'_prop'] = Results[n+'_prop'].rename(columns=dict(zip(['h_kJ/mol','s_kJ/K', 'v_cm3',
                                                                             'cp_J/kg/K', 'rho_kg/m3'],
                                                                             ['H_J', 'S_J/K', 'V_cm^3',
                                                                              'Cp_J/(kg.K^2)', 'rho_kg/m^3'])))

    for s in Results:
        if s == "Conditions":
            priority = ['T_C', 'P_bar', 'mass_g', 'H_J', 'S_J/K', 'V_cm^3', 'rho_kg/m^3', 'log10(fO2)']
            others = [c for c in Results['Conditions'].columns if c not in priority]
            Results['Conditions'] = Results['Conditions'][priority + others]
        elif s != "sys":
            if '_prop' in s:
                if "MELTS" in Model:
                    priority = ['mass_g', 'rho_kg/m^3', 'V_cm^3', 'G_J','H_J','S_J/K', 'Cp_J/(kg.K^2)',
                                'dCpdT_J/(kg.K^2)','dVdT_cm^3/K','dPdT_bar/K',
                                'd2VdT2_cm^3/K^2', 'd2VdTdP_cm^3/(bar.K)', 'd2VdP2_cm^3/bar^2']
                    others = [c for c in Results[s].columns if c not in priority]
                    Results[s] = Results[s][priority + others]
                else:
                    priority = ['mass_g', 'rho_kg/m^3','V_cm^3', 'H_J', 'S_J/K', 'Cp_J/(kg.K^2)']
                    others = [c for c in Results[s].columns if c not in priority]
                    Results[s] = Results[s][priority + others]

    Results = rename_keys_with_prefix(Results, Names_MM_replace)

    if "MELTS" not in Model:
        SN = []
        for R in Results:
            if '_prop' not in R and R != 'Conditions' and R != "sys":
                SN += [R]

    Results_Mass = pd.DataFrame(data = np.zeros((len(Results['Conditions']['T_C']), len(SN))), columns = SN)
    Results_Volume = Results_Mass.copy()
    Results_rho = Results_Mass.copy()
    for n in SN:
        Results_Mass[n] = Results[n + '_prop']['mass_g'].fillna(0.0)
        Results_Volume[n] = Results[n + '_prop']['V_cm^3'].fillna(0.0)
        Results_rho[n] = Results[n + '_prop']['rho_kg/m^3'].fillna(0.0)

    if Frac_solid is True or Frac_fluid is True:
        if Frac_solid is None:
            Results_Mass['fluid1_cumsum'] = Results_Mass['fluid1'].fillna(0.0).cumsum()
        elif Frac_fluid is None:
            for n in SN:
                if n != 'liquid1' and n!= 'fluid1' and n != 'liq1' and n != 'fl1':
                    Results_Mass[n + '_cumsum'] = Results_Mass[n].fillna(0.0).cumsum()
        else:
            for n in SN:
                if n != 'liquid1' and n != 'liq1':
                    Results_Mass[n + '_cumsum'] = Results_Mass[n].fillna(0.0).cumsum()

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
                elif any(n in R for n in Names):
                    for n in Names:
                        if n in R:
                            Results[R] = Results[R].add_suffix(Names[n])
                else:       
                    if '_prop' in R:
                        Results[R] = Results[R].add_suffix('_' + R[:-5])
                    else:
                        Results[R] = Results[R].add_suffix('_' + R)

            Results_All = pd.concat([Results_All, Results[R]], axis = 1)

    Results['All'] = Results_All
    def combine_headers(row):
        return ','.join([col[7:] for col in Results['All'].loc[:, Results['All'].columns.str.contains('mass_g_')].columns if row[col] > 0.0 and not pd.isna(row[col])])

    Results['PhaseList'] = Results['All'].apply(combine_headers, axis=1)

    Results['mass_g'] = Results_Mass.fillna(0.0)
    Results['volume_cm3'] = Results_Volume.fillna(0.0)
    Results['rho_kg/m3'] = Results_rho.fillna(0.0)

    if Results['mass_g'].sum(axis = 1).iloc[-1] == 0.0:
        for R in Results:
            Results[R].drop(Results[R].tail(1).index,inplace=True)

    return Results

def standardize_mineral_labels(output_dict):
    # Phase list definitions
    # Note: Using tuples for the abbreviations matrix for Pythonic structure
    dict_ss = {
        "sp": ("spinel", [("sp", "spinel"), ("mt", "magnetite")]),
        "spl": ("spinel", [("spl", "spinel"), ("cm", "chromite"), ("usp", "uvospinel"), ("mgt", "magnetite")]),
        "mu": ("muscovite", [("pat", "paragonite"), ("mu", "muscovite")]),
        "amp": ("amphibole", [("gl", "glaucophane"), ("act", "actinolite"), ("cumm", "cummingtonite"), ("tr", "tremolite"), ("amp", "amphibole")]),
        "ilm": ("ilmenite", [("hem", "hematite"), ("ilm", "ilmenite")]),
        "ilmm": ("ilmenite", [("hemm", "hematite"), ("ilmm", "ilmenite")]),
        "nph": ("nepheline", [("K-nph", "K-rich nepheline"), ("nph", "nepheline")]),
        "cpx": ("clinopyroxene", [("pig", "pigeonite"), ("Na-cpx", "Na-rich clinopyroxene"), ("cpx", "clinopyroxene")]),
        "dio": ("diopside", [("dio", "diopside"), ("omph", "omphacite"), ("jd", "jadeite")]),
        "occm": ("carbonate", [("sid", "siderite"), ("ank", "ankerite"), ("mag", "magnesite"), ("cc", "calcite")]),
        "oamp": ("orthoamhibole", [("anth", "anthophyllite"), ("ged", "gedrite")])
    }

    # 1. Define the Group Merges
    group_merges = {
        "sp": "spl",
        "spl": "spl",
        "cpx": "cpx",
        "dio": "cpx",
        "ilm": "ilm",
        "ilmm": "ilm",
        "amp": "amp",
        "oamp": "amp"
    }

    # 2. Build a Reverse Lookup Map
    reverse_lookup = {}
    for group_key, (name, abbr_list) in dict_ss.items():
        target_label = group_merges.get(group_key, group_key)
        for abbr, full_name in abbr_list:
            reverse_lookup[abbr] = target_label

    # 3. Process the output_dict keys
    new_dict = {}
    rename_map = {}
    group_counters = {}

    # Identify unique base keys (e.g., "pl1", "mgt1")
    all_keys = output_dict.keys()
    base_instances = sorted(list(set(
        k.replace("_prop", "") for k in all_keys 
        if k not in ["sys", "Conditions"] and "_prop" in k or not k.endswith("_prop")
    )))

    # Narrow down to specific phase instances
    for orig_instance in base_instances:
        if orig_instance in ["sys", "Conditions"]:
            continue
            
        # Extract prefix (letters) and suffix (numbers)
        match = re.match(r"([a-zA-Z\-]+)([0-9]*)", orig_instance)
        if match:
            prefix, suffix = match.groups()
            
            if prefix in reverse_lookup:
                target_group = reverse_lookup[prefix]
                
                # Increment counter
                group_counters[target_group] = group_counters.get(target_group, 0) + 1
                new_instance = f"{target_group}{group_counters[target_group]}"
                rename_map[orig_instance] = new_instance
            else:
                rename_map[orig_instance] = orig_instance

    # 4. Construct the final dictionary
    for old_key, value in output_dict.items():
        if old_key in ["sys", "Conditions"]:
            new_dict[old_key] = value
            continue

        is_prop = "_prop" in old_key
        base = old_key.replace("_prop", "")
        
        if base in rename_map:
            new_key = f"{rename_map[base]}_prop" if is_prop else rename_map[base]
            new_dict[new_key] = value
        else:
            new_dict[old_key] = value

    return new_dict