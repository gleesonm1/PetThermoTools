"""
Tier A tests for the Sun & Liang (2013) / Yao et al. (2012) low-Ca pyroxene-melt
REE+Y partitioning model (trace_engine.px_partitioning_models.Sun2013Model).

Reference: Unit_tests/golden/trace_engine/sun2013_lowca_macro_examples.csv, the
eleven worked examples in Sun's 2024 workbook
"Low-Ca_pyroxene-melt_REE-Y_partitioning_model.xlsm" (T, melt TiO2 wt%, low-Ca
pyroxene majors). P_bar is not used by the model and is set to 1e4. The blank K2O
cells in the workbook are stored as 0, which is how the macro reads them.
Two sets of D:
  - `macro_<El>`: output of Sun's spreadsheet macro (external truth)
  - `<El>`:       PTT's output when the reference was recorded (regression snapshot)

PTT uses the macro's Avogadro constant (6.02e23) and unrounded coefficients, so it
matches the macro to ~1e-5; the remainder is PTT's OXIDE_MW molar masses, which differ
from the macro's in the fourth decimal place.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import petthermotools.trace_engine.px_partitioning_models as px

_REF = Path(__file__).resolve().parents[1] / "golden" / "trace_engine" / "sun2013_lowca_macro_examples.csv"
REE = list(px.SUN2012_RADII)


@pytest.fixture(scope="module")
def ref():
    # round_trip: pandas' default fast parser can be 1 ulp off the stored value
    return pd.read_csv(_REF, float_precision="round_trip")


def _opx(ref):
    return ref[[c for c in ref.columns if c.endswith("_Opx")]]


def _sun2013(ref, liquid="ref"):
    return px.PxTraceElementCalculator().calculate_partitioning(
        data=_opx(ref),
        liquid=ref[["TiO2_Liq"]] if isinstance(liquid, str) else liquid,
        model_name="Sun2013",
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


def test_sun2013_reproduces_macro(ref):
    _assert_close(_sun2013(ref), ref, "macro_", rtol=5e-5)


def test_sun2013_matches_recorded_ptt_output(ref):
    _assert_close(_sun2013(ref), ref, "", rtol=1e-10)


def test_sun2013_output_shape_and_index(ref):
    D = _sun2013(ref)
    assert list(D.columns) == REE
    assert D.index.equals(ref.index)
    assert np.isfinite(D.to_numpy()).all() and (D.to_numpy() > 0).all()


def test_sun2013_melt_tio2_lowers_D(ref):
    """D0 carries -0.8445*X_Ti (X_Ti = 0.0304*TiO2 wt%), so a Ti-free melt must raise every D,
    by exactly exp(0.8445*0.0304*TiO2)."""
    wet = _sun2013(ref)
    no_ti = _sun2013(ref, liquid=None)
    ratio = no_ti[REE].to_numpy() / wet[REE].to_numpy()
    expected = np.exp(0.844452956724455 * 0.0304 * ref["TiO2_Liq"].to_numpy())[:, None]
    assert (ratio > 1).all()
    np.testing.assert_allclose(ratio, np.broadcast_to(expected, ratio.shape), rtol=1e-12)


def test_sun2013_iron_and_column_layout_do_not_matter(ref):
    """Same answer whether Fe is given as FeOt or FeO + Fe2O3, with a lowercase suffix and no K2O column."""
    opx = _opx(ref)
    alt = (opx.drop(columns=["FeOt_Opx", "K2O_Opx"])
              .assign(FeO_Opx=opx["FeOt_Opx"] * 0.9,
                      Fe2O3_Opx=opx["FeOt_Opx"] * 0.1 / 0.89981)
              .rename(columns=lambda c: c.replace("_Opx", "_opx")))
    D = px.PxTraceElementCalculator().calculate_partitioning(
        data=alt, liquid=ref[["TiO2_Liq"]], model_name="Sun2013",
        temperature_c=ref["T_C"], pressure_bar=ref["P_bar"])
    _assert_close(D, ref, "", rtol=1e-4)
