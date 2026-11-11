import numpy as np 
import pandas as pd 
import re
import warnings
from petthermotools.trace_engine.partition_standards import phases
from petthermotools.trace_engine.px_partitioning_models import PxTraceElementCalculator
from petthermotools.trace_engine.ol_partitioning_models import OlTraceElementCalculator
from petthermotools.trace_engine.grt_partitioning_models import GrtTraceElementCalculator


def check_partition_keys(partition):
    if type(partition) == dict:
        for p in partition.keys():
            if p not in phases:
                warnings.warn(f"Supplied dictionary of partition coefficients includes the phase {p} which is not in the alphaMELTS list of phases. These values will not be used in the following calculations.")

def get_trace_d(results, d_coeffs, partition_models):
    mass = results['mass_g'].copy()
    mass = mass.loc[:, ~mass.columns.str.contains('_cumsum')]

    trace_d = {}
    for phase in mass.columns:
        if "liquid" not in phase and "fluid" not in phase:
            # Initialize as DataFrame matching structure of d_coeffs and mass
            trace_d[phase] = pd.DataFrame(
                -1.0, index=mass.index, columns=d_coeffs.columns
            )

    if partition_models is not None:
        for key in partition_models:
            key_lower = str(key).lower()
            if key_lower in ["cpx", "clinopyroxene", "Cpx", "clinopyroxene1"]:
                cpx_calc = PxTraceElementCalculator()
                cpx_dict = {k: v for k, v in results.items() if "clinopyroxene" in k and "prop" not in k}
                if len(cpx_dict)>0.0:
                    for cpx_num in cpx_dict.keys():
                        if 'liquid1' in results.keys():
                            kd_calc = cpx_calc.calculate_partitioning(
                                    data = results[cpx_num],
                                    liquid = results['liquid1'],
                                    model_name=partition_models[key],
                                    temperature_c=results['All']['T_C'],
                                    pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        else:
                            kd_calc = cpx_calc.calculate_partitioning(
                                    data = results[cpx_num],
                                    liquid = None,
                                    model_name=partition_models[key],
                                    temperature_c=results['All']['T_C'],
                                    pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        
                        valid_mask = kd_calc.notna() & (kd_calc > 0.0)
                        trace_d[cpx_num]= trace_d[cpx_num].mask(valid_mask, kd_calc)
            elif key_lower in ["opx", "orthopyroxene", "Opx", "orthopyroxene1"]:
                opx_calc = PxTraceElementCalculator()
                opx_dict = {k: v for k, v in results.items() if "orthopyroxene" in k and "prop" not in k}
                if len(opx_dict)>0.0:
                    for opx_num in opx_dict.keys():
                        if 'liquid1' in results.keys():
                            kd_calc = opx_calc.calculate_partitioning(
                                    data = results[opx_num],
                                    liquid = results['liquid1'],
                                    model_name=partition_models[key],
                                    temperature_c=results['All']['T_C'],
                                    pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        else:
                            kd_calc = opx_calc.calculate_partitioning(
                                    data = results[opx_num],
                                    liquid = None,
                                    model_name=partition_models[key],
                                    temperature_c=results['All']['T_C'],
                                    pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        
                        valid_mask = kd_calc.notna() & (kd_calc > 0.0)
                        trace_d[opx_num]= trace_d[opx_num].mask(valid_mask, kd_calc)
            elif key_lower in ["ol", "olivine", "olivine1", "Ol", "Olivine"]:
                ol_calc = OlTraceElementCalculator()
                ol_dict = {k: v for k, v in results.items() if "olivine" in k and "prop" not in k}
                if len(ol_dict) > 0.0:
                    for ol_num in ol_dict.keys():
                        if 'liquid1' in results.keys():
                            kd_calc = ol_calc.calculate_partitioning(
                                data = results[ol_num],
                                liquid = results['liquid1'],
                                model_name=partition_models[key],
                                temperature_c=results['All']['T_C'],
                                pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        else:
                            kd_calc = ol_calc.calculate_partitioning(
                                data = results[ol_num],
                                liquid = None,
                                model_name=partition_models[key],
                                temperature_c=results['All']['T_C'],
                                pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                    
                        valid_mask = kd_calc.notna() & (kd_calc > 0.0)
                        trace_d[ol_num] = trace_d[ol_num].mask(valid_mask, kd_calc)
            elif key_lower in ["grt", "garnet", "g"]:
                grt_calc = GrtTraceElementCalculator()
                grt_dict = {k: v for k, v in results.items() if "garnet" in k and "prop" not in k}
                if len(grt_dict) > 0.0:
                    for grt_num in grt_dict.keys():
                        if 'liquid1' in results.keys():
                            kd_calc = grt_calc.calculate_partitioning(
                                data = results[grt_num],
                                liquid = results['liquid1'],
                                model_name=partition_models[key],
                                temperature_c=results['All']['T_C'],
                                pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                        else:
                            kd_calc = grt_calc.calculate_partitioning(
                                data = results[grt_num],
                                liquid = None,
                                model_name=partition_models[key],
                                temperature_c=results['All']['T_C'],
                                pressure_bar=results['All']['P_bar']).reindex(columns=d_coeffs.columns, fill_value=0.0)
                    
                        valid_mask = kd_calc.notna() & (kd_calc > 0.0)
                        trace_d[grt_num] = trace_d[grt_num].mask(valid_mask, kd_calc)

            # elif key == "opx" or key == "orthopyroxene" or key == "Opx":
            #     opx_calc = PxTraceElementCalculator()
            
            else:
                print(f"{key} not identified, using default partition coefficients instead.")

    if d_coeffs is not None:
        for phase, df in trace_d.items():
            mineral_name = re.sub(r"\d+$", "", phase)
            
            if mineral_name in d_coeffs.index:
                # Get the row of element Kd values as a 1D numpy array
                kd_values = d_coeffs.loc[mineral_name]

                trace_d[phase] = df.mask(df == -1.0, kd_values, axis = 1)

    final_trace_d = {}
    for phase, df in trace_d.items():
        cleaned_df = df.replace(-1.0, np.nan).fillna(0.0)
        final_trace_d[phase] = cleaned_df.values

    return final_trace_d

def trace_distribution(results, trace_d):
    mass = results['mass_g'].copy()
    mass = mass.loc[:, ~mass.columns.str.contains('_cumsum')]
    mass_sum = mass.div(mass.sum(axis = 1), axis = 0)

    trace_split = {}
    for d in trace_d:
        trace_split[d] = np.zeros_like(trace_d[d], dtype = float)

    first_phase = next(iter(trace_d))
    DX_sum = np.zeros_like(trace_d[first_phase], dtype = float)

    for phase_name, D_i in trace_d.items():
        if phase_name != "liquid1" and phase_name in mass.columns:
            X_i = mass_sum[phase_name].values[:, np.newaxis]
            DX_sum += D_i*X_i

    if "liquid1" in mass.columns:
        X_L = mass_sum["liquid1"].values[:, np.newaxis]
        denominator = X_L + DX_sum
        trace_split["liquid1"] = np.where(denominator > 0.0, 1.0/denominator, 0.0)
        for phase_name, D_i in trace_d.items():
            if phase_name != "liquid1":
                trace_split[phase_name] = trace_split['liquid1']*D_i
    else:
        trace_split["liquid1"] = np.zeros_like(trace_d[first_phase], dtype = float)
        trace_split["liquid1"] = np.where(DX_sum > 0.0, 1.0/DX_sum, 0.0)
        for phase_name, D_i in trace_d.items():
            if phase_name != "liquid1":
                trace_split[phase_name] = trace_split['liquid1']*D_i

    return trace_split

def get_trace_conc(results, trace_d, trace_split, df_initial, fractional):
    # use initial concentration to update the trace split to give a true concentration in each phase
    if any(results['mass_g'].columns.str.contains('_cumsum')):
        fractional = True

    trace_conc = {}
    for d in trace_split:
        trace_conc[d] = np.zeros_like(trace_split[d], dtype = float)

    if fractional:
        first = next(iter(trace_split))
        for idx, row in enumerate(trace_split[first]):
            for phase_name, tr in trace_split.items():
                trace_conc[phase_name][idx] = df_initial.iloc[0,:].values*tr[idx]
            if "residual_mass_g" not in results.keys():
                if results['mass_g'].loc[idx,'liquid1'] > 0.0:
                    df_initial.iloc[0,:] = trace_conc["liquid1"][idx]
            else:
                if results['mass_g'].loc[idx,'liquid1'] > 0.0:
                    if idx < len(results['residual_mass_g']):
                        new_initial = np.zeros_like(df_initial.iloc[0,:].values)
                        for phase_residual in results['residual_mass_g'].keys():
                            if phase_residual in trace_conc.keys():
                                new_initial += trace_conc[phase_residual][idx]*results['residual_mass_g'].loc[idx,phase_residual]

                        new_initial = new_initial/results['residual_mass_g'].sum(axis = 1)[idx]
                        df_initial.iloc[0,:] = new_initial
    else:
        for phase_name, tr in trace_split.items():
            trace_conc[phase_name] = df_initial.iloc[0,:].values*tr

    return trace_conc

def calc_accumulated_melt(results, trace_conc):
    liq_mass = results['mass_g']['liquid1']
    mass_trace_sum = trace_conc['liquid1'].mul(liq_mass, axis=0).cumsum()
    liq_sum = liq_mass.cumsum()
    return mass_trace_sum.div(liq_sum, axis=0)