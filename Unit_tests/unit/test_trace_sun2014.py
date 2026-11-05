"""
Tier A tests for the Sun & Liang (2014) olivine-melt REE+Y partitioning model
(trace_engine.ol_partitioning_models.Sun2014).

Reference: Unit_tests/golden/trace_engine/sun2014_olivine_macro_examples.csv, six
olivine-melt experiments from Beattie (1994) (T, P, olivine majors). Two sets of D:
  - `macro_<El>`: output of Sun's spreadsheet macro (external truth)
  - `<El>`:       PTT's output when the reference was recorded (regression snapshot)
The file also holds `macro_Sc`; PTT's Sun2014 returns REE+Y only, so Sc is not tested.

PTT uses the macro's Avogadro constant (6.02e23) and unrounded coefficients, so it
matches the macro to ~5e-6. The model needs no melt composition.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import petthermotools.trace_engine.ol_partitioning_models as ol
from petthermotools.trace_engine.core_trace_calcs import _FE2O3_TO_FEO

_REF = Path(__file__).resolve().parents[1] / "golden" / "trace_engine" / "sun2014_olivine_macro_examples.csv"
REE = list(ol.Sun2014.SUN2012_RADII)


@pytest.fixture(scope="module")
def ref():
    # round_trip: pandas' default fast parser can be 1 ulp off the stored value
    return pd.read_csv(_REF, float_precision="round_trip")


def _ol(ref):
    return ref[[c for c in ref.columns if c.endswith("_Ol")]]


def _sun2014(ref, data=None, **kwargs):
    return ol.OlTraceElementCalculator().calculate_partitioning(
        data=_ol(ref) if data is None else data,
        model_name="Sun2014",
        **{"temperature_c": ref["T_C"], "pressure_bar": ref["P_bar"], **kwargs},
    )


def _assert_close(D, ref, prefix, rtol):
    expected = ref[[prefix + el for el in REE]].to_numpy()
    rel = np.abs(D[REE].to_numpy() / expected - 1)
    worst = np.unravel_index(rel.argmax(), rel.shape)
    assert rel.max() <= rtol, (
        f"max relative difference {rel.max():.2e} > {rtol:.0e} "
        f"({ref['Sample'].iloc[worst[0]]}, {REE[worst[1]]})"
    )


def test_sun2014_reproduces_macro(ref):
    _assert_close(_sun2014(ref), ref, "macro_", rtol=2e-5)


def test_sun2014_matches_recorded_ptt_output(ref):
    _assert_close(_sun2014(ref), ref, "", rtol=1e-10)


def test_sun2014_output_shape_and_index(ref):
    D = _sun2014(ref)
    assert list(D.columns) == REE
    assert D.index.equals(ref.index)
    assert np.isfinite(D.to_numpy()).all() and (D.to_numpy() > 0).all()


def test_sun2014_pressure_lowers_D(ref):
    """D0 carries -0.1838*P (GPa), so +1 GPa must scale every D by exactly exp(-0.1838)."""
    D = _sun2014(ref)
    D_hi = _sun2014(ref, pressure_bar=ref["P_bar"] + 1e4)
    np.testing.assert_allclose(D_hi.to_numpy() / D.to_numpy(), np.exp(-0.183783970943972), rtol=1e-12)


def test_sun2014_ignores_the_liquid(ref):
    D = _sun2014(ref)
    with_liq = _sun2014(ref, liquid=pd.DataFrame({"SiO2_Liq": 50.0, "MgO_Liq": 8.0}, index=ref.index))
    pd.testing.assert_frame_equal(with_liq, D)


def test_sun2014_iron_and_column_layout_do_not_matter(ref):
    """Same answer whether Fe is given as FeOt or FeO + Fe2O3, with a lowercase suffix and no K2O column."""
    olv = _ol(ref)
    alt = (olv.drop(columns=["FeOt_Ol", "K2O_Ol"])
              .assign(FeO_Ol=olv["FeOt_Ol"] * 0.9,
                      Fe2O3_Ol=olv["FeOt_Ol"] * 0.1 / _FE2O3_TO_FEO)
              .rename(columns=lambda c: c.replace("_Ol", "_ol")))
    _assert_close(_sun2014(ref, data=alt), ref, "", rtol=1e-10)


def test_sun2014_empty_row_gives_nan(ref):
    """An absent phase (all-NaN row) gives a NaN row and leaves the others unchanged."""
    olv = _ol(ref).copy()
    olv.iloc[2] = np.nan
    D = _sun2014(ref, data=olv)
    assert D.iloc[2].isna().all()
    pd.testing.assert_frame_equal(D.drop(index=2), _sun2014(ref).drop(index=2))
