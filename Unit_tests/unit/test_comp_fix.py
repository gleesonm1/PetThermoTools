"""
Tier A tests for GenFuncs.comp_fix.

These are pure Python: no alphaMELTS or Julia is imported by comp_fix, so
these run everywhere (including CI on every OS/Python version).
"""
import pandas as pd
import pytest

import petthermotools as ptt

LIQUID_COMP = {
    'SiO2_Liq': 52.72,
    'TiO2_Liq': 2.08,
    'Al2O3_Liq': 13.26,
    'FeOt_Liq': 9.21,
    'MgO_Liq': 9.32,
    'CaO_Liq': 10.26,
    'Na2O_Liq': 2.24,
    'K2O_Liq': 0.43,
    'P2O5_Liq': 0.22,
    'Fe3Fet_Liq': 0.095,
    'H2O_Liq': 0.0,
}


def test_compfix_single_df_co2_override():
    """Original regression case: a scalar CO2_Liq override is applied as-is
    (excluded from the renormalization of the other oxides)."""
    df = pd.DataFrame(LIQUID_COMP, index=[0])
    out = ptt.comp_fix(Model="pMELTS", comp=df, CO2_Liq=0.2)
    assert out["CO2_Liq"][0] == pytest.approx(0.2, abs=1e-6)


def test_compfix_accepts_dict_input():
    """The dict code path should behave the same as the DataFrame path for a
    scalar CO2_Liq override."""
    out = ptt.comp_fix(Model="MELTSv1.0.2", comp=LIQUID_COMP, CO2_Liq=0.2)
    assert out["CO2_Liq"] == pytest.approx(0.2, abs=1e-6)


def test_compfix_dict_and_df_paths_agree():
    """The dict and DataFrame code paths in comp_fix are separate
    implementations (GenFuncs.py:1050-1150); this pins them together so a
    future refactor that unifies them can be checked against today's
    behaviour."""
    df_out = ptt.comp_fix(Model="MELTSv1.0.2", comp=pd.DataFrame(LIQUID_COMP, index=[0]), Fe3Fet_Liq=0.2)
    dict_out = ptt.comp_fix(Model="MELTSv1.0.2", comp=LIQUID_COMP, Fe3Fet_Liq=0.2)

    for key in dict_out:
        assert df_out[key][0] == pytest.approx(dict_out[key], rel=1e-9), key


def test_compfix_magemin_columns_exclude_mn_p_but_keep_h2o():
    """MAGEMin's oxide set (GenFuncs.py:1117-1118) drops MnO/P2O5/CO2 relative
    to the MELTS set (:1114-1115), but *does* keep H2O_Liq -- this pins down
    a discrepancy noticed during the audit (a prior code-reading claimed
    MAGEMin drops H2O; the source does not support that for comp_fix)."""
    out = ptt.comp_fix(Model="Green2025", comp=LIQUID_COMP, H2O_Liq=1.5)
    assert "H2O_Liq" in out
    assert "MnO_Liq" not in out
    assert "P2O5_Liq" not in out
    assert "CO2_Liq" not in out


def test_compfix_defaults_missing_oxides_to_zero():
    minimal = {'SiO2_Liq': 60.0, 'MgO_Liq': 5.0}
    out = ptt.comp_fix(Model="MELTSv1.0.2", comp=minimal)
    # Untouched oxides should be present and defaulted to 0.0 before
    # normalization, not silently dropped.
    assert out["TiO2_Liq"] == pytest.approx(0.0)
    assert out["CaO_Liq"] == pytest.approx(0.0)


def test_compfix_model_none_defaults_to_meltsv1_0_2():
    """Model=None is documented as defaulting to MELTSv1.0.2
    (GenFuncs.py:1047-1048). Several *other* functions in the package test
    ``"MELTS" in Model`` before applying this same default and would raise a
    TypeError on Model=None instead (e.g. Liq.py:658-661, 849); this test
    only pins down comp_fix's own (correct) handling."""
    out = ptt.comp_fix(Model=None, comp=LIQUID_COMP)
    assert "SiO2_Liq" in out
