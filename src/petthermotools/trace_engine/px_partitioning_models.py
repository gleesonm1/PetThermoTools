import numpy as np 
import pandas as pd
from abc import ABC, abstractmethod
from petthermotools.trace_engine.core_trace_calcs import *



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

# 1. Registry to store available models
MODEL_REGISTRY = {}


def register_model(name: str, notes: str = "", citations: str = ""):
    """Decorator to register a partitioning model class under a string name;
    `notes` is stored on the class and printed by trace_engine.show_trace_models()."""

    def decorator(cls):
        cls.notes = notes
        cls.citations = citations
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
@register_model("Wood1997", 
                notes = "Clinopyroxene partitioning model for 3+ cations on the M2 site from Wood and Blundy 1997.",
                citations = "Wood, B.J. and Blundy, J.D., 1997. A predictive model for rare earth element partitioning between clinopyroxene and anhydrous silicate melt. Contributions to Mineralogy and Petrology, 129(2), pp.166-181.")
class Wood1997(PartitionModel):

    def calculate(
        self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs
    ):
        # cpx = data.copy()
        # cpx = ensure_oxide_columns(cpx)
        if liquid is not None:
            liq = oxides(liquid, "liq", "total")

        n = len(data)
        temp = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15
        p = np.broadcast_to(np.asarray(pressure_bar, dtype=float).ravel(), (n,)) / 1e4

        # calc apfu composition
        a = px_apfu(data, fe_mode="as_given").reindex(
            columns=["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "Fe2O3", "MgO", "MnO", "CaO", "Na2O", "K2O"], fill_value=0.0)
        Si, Al, Fe2, Fe3, Mg, Mn, Ca = (a[k].to_numpy() for k in ("SiO2", "Al2O3", "FeO", "Fe2O3", "MgO", "MnO", "CaO"))
        M2_other = a[["CaO", "Na2O", "K2O"]].sum(axis=1).to_numpy()

        AlIV = np.clip(np.minimum(2.0 - Si, Al), 0.0, None)
        AlVI = Al - AlIV


        # get parameters needed for trace model
        XCa_M2 = Ca
        XSi = Si
        XAl_tet = AlIV
        XAl_vi = AlVI

        # absent-phase (NaN) or zero rows give NaN without divide/invalid warnings
        with np.errstate(divide="ignore", invalid="ignore"):
            MgFe = Mg/(Mg + Fe2 + Fe3) # assume they're treating Fe2+ and Fe3+ the same?
            MgFeMn = Mg/(Mg + Fe2 + Fe3 + Mn)

            XMg_M2 = (1-M2_other)*MgFeMn
            XMg_M1 = Mg - XMg_M2

            # calculate parameters used in lattice strain model
            ro = 0.974+0.067*XCa_M2 - 0.051*XAl_vi
            E = (318.6+6.9*p - 0.036*temp)*10**9
            C = 88750-65.644*temp + 7050*p-770*p**2

            MgNo_cpx = MgFe # assume they're treating Fe2+ and Fe3+ the same?
            K_femg = 0.109+0.186*MgNo_cpx
            MgNo_melt = (K_femg*MgNo_cpx)/(1-MgNo_cpx+MgNo_cpx*K_femg)
            if liquid is not None:
                FeO = liq["FeO"].to_numpy() / OXIDES.at["FeO", "mw"]
                MgO = liq["MgO"].to_numpy() / OXIDES.at["MgO", "mw"]

                MgNo_melt_true = np.broadcast_to(MgO/(MgO+FeO), (n,)) # assume that they only consider the total Mg#, i.e., ignoring ferric Fe...

                valid_mask = np.nan_to_num(MgNo_melt_true) > 0.0
                MgNo_melt[valid_mask] = MgNo_melt_true[valid_mask]

            Do = np.exp(np.log(MgNo_melt) - np.log(XMg_M1)+C/(8.3145*temp))

            ri = np.array(list(SUN2012_RADII.values()))[None, :] * 1e-10
            r0m = ro[:, None] * 1e-10
            strain = r0m / 2 * (ri - r0m)**2 + (ri - r0m)**3 / 3
            D = Do[:, None] * np.exp(-4 * np.pi * E[:, None] * N_A * strain / (R_GAS * temp[:, None]))
        return pd.DataFrame(D, columns=list(SUN2012_RADII), index=data.index)

## Extra functions required for Sun 2012
 
def sun2012_xh2o(h2o_wt):
    """Melt H2O wt% -> X_H2O as hard-coded in the macro (their fit to the
    Wood & Blundy 2002 mole fraction). NOT a simple oxide mole fraction."""
    w = np.asarray(h2o_wt, dtype=float)
    return (0.000170685747186617 * w**3 - 0.006831778563526 * w**2
            + 0.109999984719062 * w + 0.000971902247238525)
 
 
def _suffix(df, default):
    for c in df.columns:
        if isinstance(c, str) and c.startswith("SiO2_"):
            return c[c.find("_"):]
    return default

# Shannon (1976) VIII-fold radii (Angstrom) -- exactly as used in the Sun macro
SUN2012_RADII = {"La": 1.160, "Ce": 1.143, "Pr": 1.126, "Nd": 1.109, "Sm": 1.079,
                 "Eu": 1.066, "Gd": 1.053, "Tb": 1.040, "Dy": 1.027, "Y": 1.019,
                 "Ho": 1.015, "Er": 1.004, "Tm": 0.994, "Yb": 0.985, "Lu": 0.977}
 
 
@register_model("Sun2012", notes = """
                    Clinopyroxene partitioning model for 3+ cations 
                    (REE + Y) on the M2 site. Benchmarked against 
                    the excel macro file from Prof. Chenguang Sun. 
                    Original studies - Sun and Liang (2012, 2014). 
                    """,
                    citations = """
                    Sun, C. and Liang, Y., 2012. Distribution of REE between clinopyroxene and basaltic melt along a mantle adiabat: effects of major element composition, water, and temperature. Contributions to Mineralogy and Petrology, 163(5), pp.807-823.
                    Sun, C. and Liang, Y., 2014. An assessment of subsolidus re-equilibration on REE distribution among mantle minerals olivine, orthopyroxene, clinopyroxene, and garnet in peridotites. Chemical Geology, 372, pp.80-91.
                    """)
class Sun2012Model(PartitionModel):
    N_A = 6.02e23  
    def calculate(self, data, liquid=None, temperature_c=None, pressure_bar=None, **kwargs):
        n = len(data)
        T_K = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15
 
        # melt water (wt%) -> macro's X_H2O; no liquid or no H2O column = anhydrous
        h2o = np.zeros(n)
        if liquid is not None:
            # liq = strip(liquid, LIQ_SUFFIX)
            # h2o = np.broadcast_to(num(liq, "H2O").to_numpy(), (n,))
            h2o = np.broadcast_to(oxides(liquid, "liq")["H2O"].to_numpy(), (n,))
        XH2O = sun2012_xh2o(h2o)
 
        # 6-O formula, all Fe as Fe2+, Fe-Mg random over M1/M2, Mn on M2 (as macro)
        a = px_apfu(data, fe_mode="total").reindex(
            columns=["SiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O"], fill_value=0.0)
        Si, Al, Fe, Mg = (a[k].to_numpy() for k in ("SiO2", "Al2O3", "FeO", "MgO"))
        M2_other = a[["CaO", "Na2O", "K2O", "MnO"]].sum(axis=1).to_numpy() ## Mn pooled in M2 following the Sun macro
 
        AlIV = np.clip(np.minimum(2.0 - Si, Al), 0.0, None)
        AlVI = Al - AlIV
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
        D = D0[:, None] * np.exp(-4 * np.pi * E[:, None] * self.N_A * strain / (R_GAS * T_K[:, None]))
        return pd.DataFrame(D, columns=list(SUN2012_RADII), index=data.index)

@register_model("Sun2013", notes = """
                    Orthopyroxene (low-Ca pyroxene) partitioning model. Benchmarked against 
                    the excel macro file from Prof. Chenguang Sun. 
                    Original studies - Sun and Liang (2013 2014); Liang et al. (2012). 
                    """,
                    citations = """
                    Yao, L., Sun, C. and Liang, Y., 2012. A parameterized model for REE distribution between low-Ca pyroxene and basaltic melts with applications to REE partitioning in low-Ca pyroxene along a mantle adiabat and during pyroxenite-derived melt and peridotite interaction. Contributions to Mineralogy and Petrology, 164(2), pp.261-280. \n
                    Sun, C. and Liang, Y., 2013. Distribution of REE and HFSE between low-Ca pyroxene and lunar picritic melts around multiple saturation points. Geochimica et Cosmochimica Acta, 119, pp.340-358.
                    Sun, C. and Liang, Y., 2014. An assessment of subsolidus re-equilibration on REE distribution among mantle minerals olivine, orthopyroxene, clinopyroxene, and garnet in peridotites. Chemical Geology, 372, pp.80-91.
                    """)
class Sun2013Model(PartitionModel):
    N_A = 6.02e23  
    def calculate(self, data, liquid=None, temperature_c = None, pressure_bar = None, **kwargs):
        n = len(data)
        T_K = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15

        # 6-O formula, all Fe as Fe2+, Fe-Mg random over M1/M2, Mn on M2 (as macro)
        a = px_apfu(data, fe_mode="total").reindex(
            columns=["SiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O"], fill_value=0.0)
        Si, Al, Fe, Mg, Ca = (a[k].to_numpy() for k in ("SiO2", "Al2O3", "FeO", "MgO", "CaO"))
        M2_other = a[["CaO", "Na2O", "K2O", "MnO"]].sum(axis=1).to_numpy() ## Mn pooled in M2 following the Sun macro

        AlIV = np.clip(np.minimum(2.0 - Si, Al), 0.0, None)
        AlVI = Al - AlIV        
        with np.errstate(divide="ignore", invalid="ignore"):
            MgM2 = np.where(M2_other <= 1.0, (1.0 - M2_other) * Mg / (Mg + Fe), 0.0)

        Ti = np.zeros(n)
        if liquid is not None:
            Ti = 0.0304*np.broadcast_to(oxides(liquid, "liq")["TiO2"].to_numpy(), (n,))

        Do = np.exp(-5.37068370408382+(3.87339838585151*10**4)/(R_GAS*T_K) + 3.54081780772732*AlIV + 3.56148426371104*Ca - 0.844452956724455*Ti)
        ro = 0.692505524356551 + 0.431788725148506*Ca + 0.227605456716806*MgM2 
        E = (1.85482137214735*ro - 1.37247337291936 - 0.530577358996128*Ca)*10**12  # Pa

        ri = np.array(list(SUN2012_RADII.values()))[None, :] * 1e-10
        r0m = ro[:, None] * 1e-10
        strain = r0m / 2 * (ri - r0m)**2 + (ri - r0m)**3 / 3
        D = Do[:, None] * np.exp(-4 * np.pi * E[:, None] * self.N_A * strain / (R_GAS * T_K[:, None]))
        return pd.DataFrame(D, columns=list(SUN2012_RADII), index=data.index)

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