import numpy as np 
import pandas as pd
from abc import ABC, abstractmethod
from petthermotools.trace_engine.core_trace_calcs import *


# 1. Registry to store available models
MODEL_REGISTRY = {}


def register_model(name: str, notes: str = "", citations: str = ""):
    """Decorator to register a partitioning model class under a string name;
    `notes` is stored on the class and printed by trace_engine.show_models()."""

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
@register_model("Beattie1991", 
                notes = "Olivine partitioning models from Beattie et al. (1991), including Ca, Sc, Mn, Ni, Co, and Yb",
                citations = "Beattie, P., Ford, C. and Russell, D., 1991. Partition coefficients for olivine-melt and orthopyroxene-melt systems. Contributions to Mineralogy and Petrology, 109(2), pp.212-224.")
class Beattie1991(PartitionModel):

    ol_data = {
        'Element': ['Ca', 'Sc', 'Mn', 'Co', 'Ni', 'Yb'],
        'A_i_OlL': [0.0056, 0.063, 0.259, 0.786, 3.346, 0.0072],
        'B_i_OlL': [0.0135, -0.034, -0.049, -0.385, -3.665, -0.0005],
    }
    ol_AB = pd.DataFrame(ol_data).set_index('Element')

    def calculate(
        self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs
    ):
        ol = data.copy()
        ol, ol_tot = any_cation_fractions(ol, True, phase="ol")
        liquid, liq_tot = any_cation_fractions(liquid, True, phase="liq")

        # p = pressure_bar.values/10000
        # temp = temperature_c.values + 273.15

        # # calc cation fractions
        # ol_cat, ol_total_cat = any_cation_fractions(ol, return_total_cations=True)
        # liq_cat, liq_total_cat = any_cation_fractions(liquid, return_total_cations=True)

        D_Mg = ol['MgO'] / liquid['MgO']

        A = self.ol_AB['A_i_OlL']
        B = self.ol_AB['B_i_OlL']

        # D on cation-fraction basis: (N samples) x (6 elements)
        D_cation = pd.DataFrame(
            np.outer(D_Mg.values, A.values) + B.values,
            index=D_Mg.index,
            columns=A.index,
        )

        # Convert to mass (ppm) basis
        mbar_ratio = ol_tot / liq_tot  # Series, length N samples
        D_mass = D_cation.mul(mbar_ratio, axis=0)

        return D_mass

@register_model("Sun2014", notes = """
                    Olivine partitioning model for 3+ cations 
                    (REE + Y). Benchmarked against 
                    the excel macro file from Prof. Chenguang Sun. 
                    Original study - Sun and Liang (2014). 
                    """,
                    citations = """
                    Sun, C. and Liang, Y., 2014. An assessment of subsolidus re-equilibration on REE distribution among mantle minerals olivine, orthopyroxene, clinopyroxene, and garnet in peridotites. Chemical Geology, 372, pp.80-91.                    
                    """)
class Sun2014(PartitionModel):
    # Shannon (1976) VIII-fold radii (Angstrom) -- exactly as used in the Sun macro
    SUN2012_RADII = {"La": 1.160, "Ce": 1.143, "Pr": 1.126, "Nd": 1.109, "Sm": 1.079,
                    "Eu": 1.066, "Gd": 1.053, "Tb": 1.040, "Dy": 1.027, "Y": 1.019,
                    "Ho": 1.015, "Er": 1.004, "Tm": 0.994, "Yb": 0.985, "Lu": 0.977}

    N_A = 6.02e23  
    def calculate(self, data, liquid=None, temperature_c = None, pressure_bar = None, **kwargs):
        n = len(data)
        T_K = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15
        P_GPa = np.broadcast_to(np.asarray(pressure_bar, dtype=float).ravel(), (n,))*10**(-4)

        ol = any_cation_fractions(data, False, phase="ol", fe_mode="total")
        # ol_ox = oxides(data, phase = "ol")

        Mg, Fe = (ol[k].to_numpy() for k in ("MgO", "FeO"))
        Fo = (100*Mg/(Mg+Fe))

        Al = normalise(oxides(data, "ol", "total", volatiles=False), "oxygen", 4)['Al2O3'].values

        Do = np.exp(-0.435048903684216 - 0.183783970943972*P_GPa + 123.748633183663*Al - 1.48890441979127*10**(-2)*Fo)
        ro = 0.808800831576529
        E = 297.625641481357*10**9

        ri = np.array(list(self.SUN2012_RADII.values()))[None, :] * 1e-10
        r0m = ro * 1e-10
        strain = r0m / 2 * (ri - r0m)**2 + (ri - r0m)**3 / 3
        D = Do[:, None] * np.exp(-4 * np.pi * E * self.N_A * strain / (8.3145 * T_K[:, None]))
        return pd.DataFrame(D, columns=list(self.SUN2012_RADII), index=data.index)

# 4. Main Framework Class
class OlTraceElementCalculator:

    def calculate_partitioning(
        self,
        *,
        data: pd.DataFrame,
        temperature_c: np.ndarray | pd.Series,
        pressure_bar: np.ndarray | pd.Series,
        liquid: pd.DataFrame | None = None,
        model_name: str = "Beattie1991",
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