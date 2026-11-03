"""
Tier A tests for trace_engine.core_trace_calcs.liq_mole_fractions.
"""
import numpy as np
import pandas as pd
import pytest

from petthermotools.GenFuncs import OXIDE_MW
from petthermotools.trace_engine.core_trace_calcs import liq_mole_fractions

# Column order as PTT's stich_work writes a MELTS liquid (Cr2O3 before Fe2O3,
# unlike OXIDE_MW), with FeOt/Fe3Fet alongside FeO/Fe2O3.
MELTS_LIQ = {
    "SiO2_Liq": 50.0, "TiO2_Liq": 1.0, "Al2O3_Liq": 15.0, "Cr2O3_Liq": 2.0,
    "Fe2O3_Liq": 1.5, "FeO_Liq": 7.0, "FeOt_Liq": 8.35, "MnO_Liq": 0.1,
    "MgO_Liq": 8.0, "CaO_Liq": 10.0, "Na2O_Liq": 2.5, "K2O_Liq": 0.4,
    "P2O5_Liq": 0.2, "H2O_Liq": 3.0, "CO2_Liq": 0.0, "Fe3Fet_Liq": 0.16,
}


def _expected(comp):
    """Hand calculation from FeO/Fe2O3 (FeOt/Fe3Fet ignored when they are given)."""
    moles = {k[:-4]: v / OXIDE_MW[k[:-4]] for k, v in comp.items()
             if k[:-4] in OXIDE_MW}
    total = sum(moles.values())
    return {k: m / total for k, m in moles.items()}


def test_matches_hand_calculation():
    X = liq_mole_fractions(pd.DataFrame([MELTS_LIQ])).iloc[0]
    for ox, x in _expected(MELTS_LIQ).items():
        assert X[ox] == pytest.approx(x, rel=1e-12), ox
    assert X.sum() == pytest.approx(1.0, rel=1e-12)


def test_column_order_does_not_matter():
    df = pd.DataFrame([MELTS_LIQ])
    a = liq_mole_fractions(df)
    b = liq_mole_fractions(df[list(reversed(df.columns))])
    pd.testing.assert_frame_equal(a[sorted(a.columns)], b[sorted(b.columns)])


def test_non_oxide_columns_are_ignored():
    df = pd.DataFrame([MELTS_LIQ]).assign(Sample="x", T_C=1200.0)
    pd.testing.assert_frame_equal(liq_mole_fractions(df), liq_mole_fractions(pd.DataFrame([MELTS_LIQ])))


def test_feot_only_liquid_keeps_its_iron():
    comp = {"SiO2_Liq": 50.0, "FeOt_Liq": 10.0, "MgO_Liq": 10.0, "H2O_Liq": 2.0}
    X = liq_mole_fractions(pd.DataFrame([comp])).iloc[0]
    m = {"SiO2": 50.0 / OXIDE_MW["SiO2"], "FeO": 10.0 / OXIDE_MW["FeO"],
         "MgO": 10.0 / OXIDE_MW["MgO"], "H2O": 2.0 / OXIDE_MW["H2O"]}
    total = sum(m.values())
    for ox, v in m.items():
        assert X[ox] == pytest.approx(v / total, rel=1e-12), ox
    assert X["Fe2O3"] == 0.0


def test_empty_rows_give_nan_not_error():
    df = pd.DataFrame([MELTS_LIQ, {k: np.nan for k in MELTS_LIQ}])
    X = liq_mole_fractions(df)
    assert np.isfinite(X.iloc[0]).all() and X.iloc[1].isna().all()
