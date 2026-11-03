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
@register_model("Sun_2012_cpx")
class SunLiang2012_cpx(PartitionModel):

    def calculate(
        self, data, liquid, temperature_c: float, pressure_bar: float, **kwargs
    ):
        XH2O = liquid['H2O_Liq']/18.01

        
        # Example calculation placeholder logic for Wood & Blundy (1997) lattice strain
        r_0 = kwargs.get("r_0", 1.0)  # Strain-free ionic radius
        # Insert model-specific equation here...
        return 0.15 * (temperature_k / 1000.0) / pressure_gpa


# @register_model("Sun_2012_cpx")
# class Sun2012Model(PartitionModel):

#     def calculate(
#         self, element: str, temperature_k: float, pressure_gpa: float, **kwargs
#     ) -> float:
#         # Example calculation placeholder logic for Sun & Liang (2012)
#         return 0.42 * (temperature_k / 1200.0)


# 4. Main Framework Class for the user
class TraceElementCalculator:

    def calculate_partitioning(
        self,
        data,
        model_name: str,
        temperature_k: float,
        pressure_gpa: float,
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
            data, temperature_k, pressure_gpa, **kwargs
        )