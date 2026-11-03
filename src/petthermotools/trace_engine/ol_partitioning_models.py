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
@register_model("Beattie1993")
class Beattie1993(PartitionModel):

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
        ol = ensure_oxide_columns(ol)
        liquid = ensure_oxide_columns(liquid.copy())

        p = pressure_bar.values/10000
        temp = temperature_c.values + 273.15

        # calc cation fractions
        ol_cat, ol_total_cat = any_cation_fractions(ol, return_total_cations=True)
        liq_cat, liq_total_cat = any_cation_fractions(liquid, return_total_cations=True)

        D_Mg = ol_cat['MgO'] / liq_cat['MgO']

        A = self.ol_AB['A_i_OlL']
        B = self.ol_AB['B_i_OlL']

        # D on cation-fraction basis: (N samples) x (6 elements)
        D_cation = pd.DataFrame(
            np.outer(D_Mg.values, A.values) + B.values,
            index=D_Mg.index,
            columns=A.index,
        )

        # Convert to mass (ppm) basis
        mbar_ratio = ol_total_cat / liq_total_cat   # Series, length N samples
        D_mass = D_cation.mul(mbar_ratio, axis=0)

        return D_mass

# 4. Main Framework Class for the user
class OlTraceElementCalculator:

    def calculate_partitioning(
        self,
        *,
        data: pd.DataFrame,
        temperature_c: np.ndarray | pd.Series,
        pressure_bar: np.ndarray | pd.Series,
        liquid: pd.DataFrame | None = None,
        model_name: str = "Beattie1993",
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