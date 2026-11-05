"""
Tier A tests for the Sun & Liang (2013) garnet-melt REE+Y partitioning model
(trace_engine.grt_partitioning_models.Sun2013).

Reference: Unit_tests/golden/trace_engine/sun2013_garnet_macro_examples.csv, nine
garnet examples (T, P, garnet majors with FeO and Fe2O3 given separately): the
macro's two worked examples (Adam & Green 2006, 1955 and 1956) and seven synthetic
variations of 1955 run through a Python port of the macro (CalGrtKd_V1), covering
1.5-6 GPa, 1100-1550 C, and one garnet with 1.5 wt% Fe2O3. Two sets of D:
  - `macro_<El>`: output of Sun's macro, or its port (external truth)
  - `<El>`:       PTT's output when the reference was recorded (regression snapshot)

PTT uses the macro's Avogadro constant (6.02e23) and unrounded coefficients, so it
matches the macro to ~2.5e-5. The model needs no melt composition.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import petthermotools.trace_engine.grt_partitioning_models as grt
from petthermotools.trace_engine.core_trace_calcs import _FE2O3_TO_FEO

_REF = Path(__file__).resolve().parents[1] / "golden" / "trace_engine" / "sun2013_garnet_macro_examples.csv"
REE = list(grt.Sun2013.SUN2012_RADII)


@pytest.fixture(scope="module")
def ref():
    # round_trip: pandas' default fast parser can be 1 ulp off the stored value
    return pd.read_csv(_REF, float_precision="round_trip")


def _grt(ref):
    return ref[[c for c in ref.columns if c.endswith("_Grt")]]


def _sun2013_grt(ref, data=None, **kwargs):
    return grt.GrtTraceElementCalculator().calculate_partitioning(
        data=_grt(ref) if data is None else data,
        model_name="Sun2013",
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


def test_sun2013_garnet_reproduces_macro(ref):
    _assert_close(_sun2013_grt(ref), ref, "macro_", rtol=5e-5)


def test_sun2013_garnet_matches_recorded_ptt_output(ref):
    _assert_close(_sun2013_grt(ref), ref, "", rtol=1e-10)


def test_sun2013_garnet_output_shape_and_index(ref):
    D = _sun2013_grt(ref)
    assert list(D.columns) == REE
    assert D.index.equals(ref.index)
    assert np.isfinite(D.to_numpy()).all() and (D.to_numpy() > 0).all()


def test_sun2013_garnet_pressure_term(ref):
    """D0 carries -93.02*P*(37.78 - P)/RT (P in GPa), so +1 GPa must scale every D by exactly
    exp(f(P + 1) - f(P)), the same factor for every element."""
    D = _sun2013_grt(ref)
    D_hi = _sun2013_grt(ref, pressure_bar=ref["P_bar"] + 1e4)
    P, T_K = ref["P_bar"].to_numpy() / 1e4, ref["T_C"].to_numpy() + 273.15
    f = lambda P: -93.0176694902273 * P * (37.78 - P) / (8.3145 * T_K)
    expected = np.exp(f(P + 1) - f(P))[:, None]
    np.testing.assert_allclose(D_hi.to_numpy() / D.to_numpy(), np.broadcast_to(expected, D.shape), rtol=1e-12)


def test_sun2013_garnet_ignores_the_liquid(ref):
    D = _sun2013_grt(ref)
    with_liq = _sun2013_grt(ref, liquid=pd.DataFrame({"SiO2_Liq": 50.0, "MgO_Liq": 8.0}, index=ref.index))
    pd.testing.assert_frame_equal(with_liq, D)


def test_sun2013_garnet_iron_and_column_layout_do_not_matter(ref):
    """Same answer with FeOt + Fe3Fet instead of FeO + Fe2O3, a lowercase suffix and no K2O column."""
    g = _grt(ref)
    feot = g["FeO_Grt"] + _FE2O3_TO_FEO * g["Fe2O3_Grt"]
    alt = (g.drop(columns=["FeO_Grt", "Fe2O3_Grt", "K2O_Grt"])
            .assign(FeOt_Grt=feot, Fe3Fet_Grt=_FE2O3_TO_FEO * g["Fe2O3_Grt"] / feot)
            .rename(columns=lambda c: c.replace("_Grt", "_grt")))
    _assert_close(_sun2013_grt(ref, data=alt), ref, "", rtol=1e-10)


def test_sun2013_garnet_uses_ferric_iron(ref):
    """Fe3+ is kept on the 12-oxygen basis, so recasting the ferric garnet's Fe2O3 as FeO
    must change its D, and only its D (every other row has no Fe2O3). The change in XCa
    moves D0 and r0 in opposite directions for the light REE (Pr barely moves), so the
    check is on the largest change in that row, ~2.4e-3."""
    g = _grt(ref)
    ferric = (g["Fe2O3_Grt"] > 0).to_numpy()
    assert ferric.sum() == 1, "reference should contain exactly one ferric garnet"
    all_feo = g.assign(FeO_Grt=g["FeO_Grt"] + _FE2O3_TO_FEO * g["Fe2O3_Grt"], Fe2O3_Grt=0.0)
    rel = np.abs(_sun2013_grt(ref, data=all_feo).to_numpy() / _sun2013_grt(ref).to_numpy() - 1)
    assert rel[ferric].max() > 1e-3
    assert (rel[~ferric] == 0).all()


def test_sun2013_garnet_empty_row_gives_nan(ref):
    """An absent phase (all-NaN row) gives a NaN row and leaves the others unchanged."""
    g = _grt(ref).copy()
    g.iloc[2] = np.nan
    D = _sun2013_grt(ref, data=g)
    assert D.iloc[2].isna().all()
    pd.testing.assert_frame_equal(D.drop(index=2), _sun2013_grt(ref).drop(index=2))
