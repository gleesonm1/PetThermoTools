import numpy as np 
import pandas as pd 
from petthermotools.trace_engine.partition_standards import *
from petthermotools.trace_engine.trace_core import *

def calc_trace_conc(results=None, elements=None, partition_models=None, D_coeffs=coeffs["D_Gibson_2010"],
                    initial_conc = None, fractional = False):

    # check element list
    if elements is None:
        if isinstance(initial_conc, dict):
            elements = list(initial_conc.keys())
        else:
            elements = common_elements

    # check and format initial conc
    if initial_conc is None:
        initial_conc = np.ones(len(elements))
    
    if isinstance(initial_conc, (np.ndarray,list)):
        initial_conc = dict(zip(elements, initial_conc))
    
    initial_conc = {
        el: val for el, val in initial_conc.items() if el in elements
    }

    df_initial = pd.DataFrame([initial_conc])
    df_initial = df_initial.reindex(columns = elements, fill_value = 1.0).fillna(0)

    # check and format partition coeffs
    if isinstance(D_coeffs, str):
        D_coeffs = coeffs[D_coeffs]

    if D_coeffs is None: ## setup background/standard partitioning behaviour
        D_coeffs = coeffs["D_Gibson_2010"]

    check_partition_keys(D_coeffs)
    if type(D_coeffs) == dict:
        d_coeffs = pd.DataFrame.from_dict(D_coeffs, orient="index")

    d_coeffs = d_coeffs.reindex(columns=elements, fill_value=0).fillna(0)
    

    def process_trace_calcs(res):
        trace_d = get_trace_d(res, d_coeffs, partition_models)
        trace_split = trace_distribution(res, trace_d)
        conc_dict = get_trace_conc(
            res, trace_d, trace_split, df_initial, fractional
        )

        mass_df = res["mass_g"]

        # Mask rows where mass is 0.0 or NaN to 0.0 concentration
        for phase_name in list(conc_dict.keys()):
            if phase_name in mass_df.columns:
                phase_mass = mass_df[phase_name].values
                # True where mass is 0, NaN, or negative
                zero_mask = np.isnan(phase_mass) | (phase_mass <= 0.0)
                # Set all element concentrations in those rows to 0.0
                conc_dict[phase_name][zero_mask, :] = 0.0

        # Calculate total mass-weighted bulk trace sum
        first = next(iter(conc_dict))
        sum_trace = np.zeros_like(conc_dict[first], dtype=float)
        for phase_name, tr in conc_dict.items():
            if phase_name in mass_df.columns:
                sum_trace += tr * mass_df[phase_name].values[:, np.newaxis]

        conc_dict["sum"] = sum_trace

        # Convert matrices to DataFrames
        for d in conc_dict:
            conc_dict[d] = pd.DataFrame(
                data=conc_dict[d], columns=d_coeffs.columns
            )

        # Convert matrices to DataFrames
        for d in trace_d:
            trace_d[d] = pd.DataFrame(
                data=trace_d[d], columns=d_coeffs.columns
            )

        return conc_dict, trace_d

    # calculate how trace elements are split between the relevant phases
    if isinstance(results, dict) and "mass_g" in results.keys():
        trace_conc, trace_ds = process_trace_calcs(results)

    elif isinstance(results, dict) and "mass_g" not in results.keys():
        trace_conc = {}
        trace_ds = {}
        for i in results:
            if isinstance(results[i], dict) and "mass_g" in results[i].keys():
                trace_conc[i], trace_ds[i] = process_trace_calcs(results[i])
    else:
        print("Unable to find results dictionary")
        return None

    return trace_conc, trace_ds

