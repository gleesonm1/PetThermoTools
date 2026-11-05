"""
Tier A tests for the Sun & Liang (2012) cpx-melt REE+Y partitioning model
(trace_engine.px_partitioning_models.Sun2012Model).

Reference: Unit_tests/golden/trace_engine/sun2012_macro_examples.csv, ten
cpx-melt examples (T, P, melt H2O wt%, cpx majors) with two sets of D:
  - `macro_<El>`: output of Sun & Liang's own spreadsheet macro (external truth)
  - `<El>`:       PTT's output when the reference was recorded (regression snapshot)

PTT uses the macro's Avogadro constant (6.02e23) and unrounded coefficients, so
it matches the macro to ~1e-5 (the rounding of the macro output).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import petthermotools.trace_engine.px_partitioning_models as px

_REF = Path(__file__).resolve().parents[1] / "golden" / "trace_engine" / "sun2012_macro_examples.csv"
REE = list(px.SUN2012_RADII)


@pytest.fixture(scope="module")
def ref():
    # round_trip: pandas' default fast parser can be 1 ulp off the stored value
    return pd.read_csv(_REF, float_precision="round_trip")


def _sun2012(ref):
    return px.PxTraceElementCalculator().calculate_partitioning(
        data=ref[[c for c in ref.columns if c.endswith("_Cpx")]],
        liquid=ref[["H2O_Liq"]],
        model_name="Sun2012",
        temperature_c=ref["T_C"],
        pressure_bar=ref["P_bar"],
    )


def _assert_close(D, ref, prefix, rtol):
    expected = ref[[prefix + el for el in REE]].to_numpy()
    rel = np.abs(D[REE].to_numpy() / expected - 1)
    worst = np.unravel_index(rel.argmax(), rel.shape)
    assert rel.max() <= rtol, (
        f"max relative difference {rel.max():.2e} > {rtol:.0e} "
        f"({ref['Sample'].iloc[worst[0]]}, {REE[worst[1]]})"
    )


def test_sun2012_reproduces_macro(ref):
    _assert_close(_sun2012(ref), ref, "macro_", rtol=5e-5)


def test_sun2012_matches_recorded_ptt_output(ref):
    _assert_close(_sun2012(ref), ref, "", rtol=1e-10)


def test_sun2012_output_shape_and_index(ref):
    D = _sun2012(ref)
    assert list(D.columns) == REE
    assert D.index.equals(ref.index)
    assert np.isfinite(D.to_numpy()).all() and (D.to_numpy() > 0).all()


def test_sun2012_water_lowers_D(ref):
    """D0 carries -0.91*X_H2O, so removing the melt (anhydrous) must raise every D."""
    wet = _sun2012(ref)
    dry = px.PxTraceElementCalculator().calculate_partitioning(
        data=ref[[c for c in ref.columns if c.endswith("_Cpx")]],
        model_name="Sun2012",
        temperature_c=ref["T_C"],
        pressure_bar=ref["P_bar"],
    )
    assert (dry[REE].to_numpy() > wet[REE].to_numpy()).all()
