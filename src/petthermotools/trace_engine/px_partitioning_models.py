import numpy as np 
import pandas as pd
from abc import ABC, abstractmethod
from petthermotools.trace_engine.core_trace_calcs import *


# 1. Registry to store available models
MODEL_REGISTRY = {}


def register_model(name: str):
    """Decorator to register a partitioning model class under a string name."""

    def decorator(cls):
        MODEL_REGISTRY[name] = cls
        return cls

    return decorator


# 2. Abstract Base Class for all models
class PartitionModel(ABC):

    @abstractmethod
    def calculate(
        self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs
    ):
        """Calculate partition coefficient (D). Must be implemented by subclasses."""
        pass


# 3. Individual Model Implementations
@register_model("Wood1997")
class Wood1997(PartitionModel):

    def calculate(
        self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs
    ):
        cpx = data.copy()
        cpx = ensure_oxide_columns(cpx)
        liquid = ensure_oxide_columns(liquid.copy())

        p = pressure_bar.values/10000
        temp = temperature_c.values + 273.15

        # calc apfu composition
        cpx_apfu = px_apfu(cpx)

        # get parameters needed for trace model
        XCa_M2 = cpx_apfu['CaO'].values
        XSi = cpx_apfu['SiO2'].values
        XAl_tet = 2 - XSi
        XAl_vi = cpx_apfu['Al2O3'].values - XAl_tet

        M2 = 1 - cpx_apfu['CaO'].values - cpx_apfu['Na2O'].values - cpx_apfu['K2O'].values
        MgFe = cpx_apfu['MgO']/(cpx_apfu['MgO'] + cpx_apfu['FeO'] + cpx_apfu['Fe2O3']) # assume they're treating Fe2+ and Fe3+ the same?

        XMg_M2 = M2*MgFe.values
        XMg_M1 = cpx_apfu['MgO'].values - XMg_M2

        # calculate parameters used in lattice strain model
        ro = 0.974+0.067*XCa_M2 - 0.051*XAl_vi
        E = (318.6+6.9*p - 0.036*temp)*10**9
        C = 88750-65.644*temp + 7050*p-770*p**2

        MgNo_cpx = (cpx_apfu['MgO'])/(cpx_apfu['MgO'] + cpx_apfu['FeO'] + cpx_apfu['Fe2O3']) # assume they're treating Fe2+ and Fe3+ the same?
        K_femg = 0.109+0.186*MgNo_cpx.values
        MgNo_melt = (K_femg*MgNo_cpx.values)/(1-MgNo_cpx.values+MgNo_cpx.values*K_femg)
        if liquid is not None:
            FeO = liquid['FeOt_Liq'].values/71.844
            MgO = liquid['MgO_Liq'].values/40.304

            MgNo_melt_true = MgO/(MgO+FeO) # assume that they only consider the total Mg#, i.e., ignoring ferric Fe...
        
            valid_mask = np.nan_to_num(MgNo_melt_true) > 0.0
            MgNo_melt[valid_mask] = MgNo_melt_true[valid_mask]

        Do = np.exp(np.log(MgNo_melt) - np.log(XMg_M1)+C/(8.3145*temp))

        ro_e = ro*10**(-10)
        rj = np.repeat(trivalent_radii_df.values, len(data), axis=0)*10**(-10)

        N = 6.0221408*10**23
        trivalent_radii_names = pd.DataFrame.from_dict(trivalent_radii, orient='index').T.columns
        Kd = Do.reshape(-1, 1)*np.exp(((-4*np.pi*E.reshape(-1, 1)*N)/(8.3145*(temp.reshape(-1, 1))))*((ro_e.reshape(-1, 1)/2)*(rj - ro_e.reshape(-1, 1))**2 + (1/3)*(rj-ro_e.reshape(-1, 1))**3))

        return pd.DataFrame(data = Kd, columns = trivalent_radii_names)

R_GAS = 8.3145
N_A = 6.0221408e23   # macro uses 6.02e23 ("602"); difference in D <= 0.09 %
 

 
# Oxide molar masses as in the Sun 2012 macro (O = 15.9994)
_O = 15.9994
_OXIDES = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeOt", "MnO", "MgO", "CaO", "Na2O", "K2O"]
_MW = np.array([28.0855 + 2*_O, 47.867 + 2*_O, 26.9815*2 + 3*_O, 51.9961*2 + 3*_O,
                55.845 + _O, 54.938 + _O, 24.305 + _O, 40.078 + _O,
                22.9898*2 + _O, 39.0983*2 + _O])
_NCAT = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2])
_NOX = np.array([2, 2, 3, 3, 1, 1, 1, 1, 1, 1])
 
 
def sun2012_xh2o(h2o_wt):
    """Melt H2O wt% -> X_H2O as hard-coded in the macro (their fit to the
    Wood & Blundy 2002 mole fraction). NOT a simple oxide mole fraction."""
    w = np.asarray(h2o_wt, dtype=float)
    return (0.000170685747186617 * w**3 - 0.006831778563526 * w**2
            + 0.109999984719062 * w + 0.000971902247238525)
 
 
def sun2012_sites(cpx_wt, suffix="_Cpx", macro_si_quirk=True):
    """Cpx site occupancies on a 6-oxygen basis, all Fe as Fe2+, Fe-Mg random
    over M1/M2 (Wood & Banno 1973), Mn on M2 -- as in the macro's PyxForm()."""
    wt = np.column_stack([np.asarray(cpx_wt[f"{ox}{suffix}"], dtype=float)
                          for ox in _OXIDES])
    oxy = wt / _MW * _NOX
    cat = oxy * (6.0 / oxy.sum(axis=1, keepdims=True)) / _NOX * _NCAT
    Si, Ti, Al, Cr, Fe, Mn, Mg, Ca, Na, K = cat.T
 
    AlIV = np.clip(np.minimum(2.0 - Si, Al), 0.0, None)
    AlVI = Al - AlIV
    if macro_si_quirk:
        # Macro sets AlIV = AlVI = 0 when Si >= 2 (only possible for Si-rich analyses)
        AlVI = np.where(Si >= 2.0, 0.0, AlVI)
 
    M2_other = Ca + Na + K + Mn
    MgM2 = np.where(M2_other <= 1.0, (1.0 - M2_other) * Mg / (Mg + Fe), 0.0)
    return AlIV, AlVI, MgM2
 
 
def _suffix(df, default):
    for c in df.columns:
        if isinstance(c, str) and c.startswith("SiO2_"):
            return c[c.find("_"):]
    return default

# Shannon (1976) VIII-fold radii (Angstrom) -- exactly as used in the Sun macro
SUN2012_RADII = {"La": 1.160, "Ce": 1.143, "Pr": 1.126, "Nd": 1.109, "Sm": 1.079,
                 "Eu": 1.066, "Gd": 1.053, "Tb": 1.040, "Dy": 1.027, "Y": 1.019,
                 "Ho": 1.015, "Er": 1.004, "Tm": 0.994, "Yb": 0.985, "Lu": 0.977}
 
 
@register_model("Sun2012")
class Sun2012Model(PartitionModel):
 
    def calculate(self, data, liquid=None, temperature_c=None, pressure_bar=None, **kwargs):
        n = len(data)
        T_K = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15
 
        # melt water (wt%) -> macro's X_H2O; no liquid or no H2O column = anhydrous
        h2o = np.zeros(n)
        if liquid is not None:
            liq = strip(liquid, LIQ_SUFFIX)
            h2o = np.broadcast_to(num(liq, "H2O").to_numpy(), (n,))
        XH2O = sun2012_xh2o(h2o)
 
        # 6-O formula, all Fe as Fe2+, Fe-Mg random over M1/M2, Mn on M2 (as macro)
        a = px_apfu(data, fe_mode="total").reindex(
            columns=["SiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O"], fill_value=0.0)
        Si, Al, Fe, Mg = (a[k].to_numpy() for k in ("SiO2", "Al2O3", "FeO", "MgO"))
        M2_other = a[["CaO", "Na2O", "K2O", "MnO"]].sum(axis=1).to_numpy()
 
        AlIV = np.clip(np.minimum(2.0 - Si, Al), 0.0, None)
        AlVI = np.where(Si >= 2.0, 0.0, Al - AlIV)        # macro: no Al assigned if Si >= 2
        with np.errstate(divide="ignore", invalid="ignore"):
            MgM2 = np.where(M2_other <= 1.0, (1.0 - M2_other) * Mg / (Mg + Fe), 0.0)
 
        # Eqs. 8-10 with the macro's unrounded coefficients
        D0 = np.exp(-7.13615731281538 + 71864.8736892434 / (R_GAS * T_K)
                    + 4.37096420622614 * AlIV + 1.9813056415739 * MgM2
                    - 0.908067851839806 * XH2O)
        r0 = 1.06596147404798 - 0.103654121912384 * AlVI - 0.211803756387744 * MgM2
        E = (2.27190977329864 * r0 - 1.99606952151084) * 1e12                 # Pa
 
        ri = np.array(list(SUN2012_RADII.values()))[None, :] * 1e-10
        r0m = r0[:, None] * 1e-10
        strain = r0m / 2 * (ri - r0m)**2 + (ri - r0m)**3 / 3
        D = D0[:, None] * np.exp(-4 * np.pi * E[:, None] * N_A * strain / (R_GAS * T_K[:, None]))
        return pd.DataFrame(D, columns=list(SUN2012_RADII), index=data.index)

# @register_model("Sun2012")
# class Sun2012Model(PartitionModel):

#     def calculate(self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs):
#         cpx = data.copy()
#         cpx = ensure_oxide_columns(cpx)

#         suffix = None
#         for col in cpx.columns:
#             if '_' in col:
#                 suffix = col[col.find('_'):]
#                 break
#         suffix = suffix or '_Cpx'

#         liquid = ensure_oxide_columns(liquid.copy())

#         p = pressure_bar.values/10000
#         temp = temperature_c.values + 273.15

#         ## fix all FeO to FeO and remove Fe2O3 from consideration
#         cpx[f"FeO{suffix}"] = cpx[f"FeOt{suffix}"]
#         cpx[f"Fe2O3{suffix}"] = np.zeros_like(cpx[f"FeO{suffix}"])

#         cpx_apfu = px_apfu(cpx)

#         # get parameters

#         XSi = cpx_apfu['SiO2'].values
#         XAl_tet = 2 - XSi
#         XAl_vi = cpx_apfu['Al2O3'].values - XAl_tet 

#         M2 = 1 - cpx_apfu['CaO'].values - cpx_apfu['Na2O'].values - cpx_apfu['K2O'].values
#         MgFe = cpx_apfu['MgO']/(cpx_apfu['MgO'] + cpx_apfu['FeO'])

#         XMg_M2 = M2*MgFe.values

#         if liquid is not None:
#             XH2O = 1.7069*10**(-4)*liquid['H2O_Liq'].values**3 - 6.8318*10**(-3)*liquid['H2O_Liq'].values**2 + 0.11*liquid['H2O_Liq'] + 9.719*10**(-4)
#         else:
#             XH2O = np.zeros(len(XSi))

#         lnDo = -7.14+((7.19*10**4)/(8.3145*(temp))) + 4.37*XAl_tet + 1.98*XMg_M2 -0.91*XH2O
#         Do = np.exp(lnDo)
#         ro = (1.066 - 0.104*XAl_vi-0.212*XMg_M2)
#         E = ((2.27191*ro - 1.99607)*10**3)*1e9
#         N = 6.0221408*10**23

#         ro_e = ro*10**(-10)
#         rj = np.repeat(trivalent_radii_df.values, len(data), axis=0)*10**(-10)

#         trivalent_radii_names = pd.DataFrame.from_dict(trivalent_radii, orient='index').T.columns
#         Kd = Do.reshape(-1, 1)*np.exp(((-4*np.pi*E.reshape(-1, 1)*N)/(8.3145*(temp.reshape(-1, 1))))*((ro_e.reshape(-1, 1)/2)*(rj - ro_e.reshape(-1, 1))**2 + (1/3)*(rj-ro_e.reshape(-1, 1))**3))

#         return pd.DataFrame(data = Kd, columns = trivalent_radii_names)

# @register_model("Sun_2012_cpx")
# class Sun2012Model(PartitionModel):

#     def calculate(
#         self, element: str, temperature_k: float, pressure_gpa: float, **kwargs
#     ) -> float:
#         # Example calculation placeholder logic for Sun & Liang (2012)
#         return 0.42 * (temperature_k / 1200.0)


# 4. Main Framework Class for the user
class PxTraceElementCalculator:

    def calculate_partitioning(
        self,
        *,
        data: pd.DataFrame,
        temperature_c: np.ndarray | pd.Series,
        pressure_bar: np.ndarray | pd.Series,
        liquid: pd.DataFrame | None = None,
        model_name: str = "Wood1997",
        **kwargs,
    ):
        """Main user-facing function."""
        if model_name not in MODEL_REGISTRY:
            available = list(MODEL_REGISTRY.keys())
            raise ValueError(
                f"Model '{model_name}' not found. Available models: {available}"
            )

        # Instantiate the class dynamically from the registry
        model_cls = MODEL_REGISTRY[model_name]
        model_instance = model_cls()

        # Run calculation
        return model_instance.calculate(
            data, liquid, temperature_c, pressure_bar, **kwargs
        )