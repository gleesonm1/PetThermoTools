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

@register_model("Sun2013", notes = """
                    Garnet partitioning model for 3+ cations 
                    (REE + Y). Benchmarked against 
                    the excel macro file from Prof. Chenguang Sun. 
                    Original studies - Sun and Liang (2013; 2014). 
                    """,
                    citations = """
                    Sun, C. and Liang, Y., 2013. The importance of crystal chemistry on REE partitioning between mantle minerals (garnet, clinopyroxene, orthopyroxene, and olivine) and basaltic melts. Chemical Geology, 358, pp.23-36.
                    Sun, C. and Liang, Y., 2014. An assessment of subsolidus re-equilibration on REE distribution among mantle minerals olivine, orthopyroxene, clinopyroxene, and garnet in peridotites. Chemical Geology, 372, pp.80-91.                    
                    """)
class Sun2013(PartitionModel):
    # Shannon (1976) VIII-fold radii (Angstrom) -- exactly as used in the Sun macro
    SUN2012_RADII = {"La": 1.160, "Ce": 1.143, "Pr": 1.126, "Nd": 1.109, "Sm": 1.079,
                    "Eu": 1.066, "Gd": 1.053, "Tb": 1.040, "Dy": 1.027, "Y": 1.019,
                    "Ho": 1.015, "Er": 1.004, "Tm": 0.994, "Yb": 0.985, "Lu": 0.977}

    N_A = 6.02e23  
    def calculate(self, data, liquid=None, temperature_c = None, pressure_bar = None, **kwargs):
        n = len(data)
        T_K = np.broadcast_to(np.asarray(temperature_c, dtype=float).ravel(), (n,)) + 273.15
        P_GPa = np.broadcast_to(np.asarray(pressure_bar, dtype=float).ravel(), (n,))*10**(-4)

        grt = normalise(oxides(data, "grt", "as_given", volatiles=False), "oxygen", 12)
        
        Do = np.exp(-2.00867101293334 + (9.03205750936295*10**4 - 93.0176694902273*P_GPa*(37.78-P_GPa))/(8.3145*T_K) - 1.04210813433143*grt['CaO'].values)
        ro = 0.784816766995275 + 0.152525704709021*grt['CaO'].values
        E = (-1.66719609963883+2.3445036098238*ro)*10**12

        ri = np.array(list(self.SUN2012_RADII.values()))[None, :] * 1e-10
        r0m = ro[:, None] * 1e-10
        strain = r0m / 2 * (ri - r0m)**2 + (ri - r0m)**3 / 3
        D = Do[:, None] * np.exp(-4 * np.pi * E[:, None] * self.N_A * strain / (8.3145 * T_K[:, None]))
        return pd.DataFrame(D, columns=list(self.SUN2012_RADII), index=data.index)


# 4. Main Framework Class
class GrtTraceElementCalculator:

    def calculate_partitioning(
        self,
        *,
        data: pd.DataFrame,
        temperature_c: np.ndarray | pd.Series,
        pressure_bar: np.ndarray | pd.Series,
        liquid: pd.DataFrame | None = None,
        model_name: str = "Sun2013",
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