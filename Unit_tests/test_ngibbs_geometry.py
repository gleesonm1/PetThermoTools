"""
Smoke / geometry tests for the nGibbs neural-network emulator API.

All tests are accessed exclusively through ptt.multi_path(), matching real
end-user usage.  They verify:
  - No unhandled exceptions
  - Output is a dict (or dict-of-dicts for grid / multi-override modes)
  - Conditions DataFrame has the expected number of rows

Accuracy of predictions is NOT tested here.

Run with:
    pytest Unit_tests/test_ngibbs_geometry.py -v
    pytest Unit_tests/test_ngibbs_geometry.py -v -m "not slow"
"""

import pytest
import numpy as np
import pandas as pd

ngibbs = pytest.importorskip("ngibbs", reason="ngibbs package not installed")

from petthermotools.Path import multi_path

# ── Reference basalt composition (PTT-style _Liq headers) ────────────────────
# Uses Fe3Fet_Liq so that iron speciation is always defined.
BASALT = {
    "SiO2_Liq":   49.5,
    "TiO2_Liq":    2.7,
    "Al2O3_Liq":  13.1,
    "Cr2O3_Liq":   0.03,
    "FeO_Liq":    10.5,
    "Fe3Fet_Liq":  0.10,
    "MnO_Liq":     0.18,
    "MgO_Liq":     7.2,
    "CaO_Liq":    11.1,
    "Na2O_Liq":    2.5,
    "K2O_Liq":     0.5,
    "P2O5_Liq":    0.25,
    "H2O_Liq":     0.5,
}

MODEL          = "nMELTSv1.0.2NoProp"   # default: MELTS 1.0.2 (only version with trained models)
MODEL_102      = "nMELTSv1.0.2NoProp"
MODEL_WITH_PROPS = "nMELTSv1.0.2"       # requires meltsdynamic + MELTS 1.2.0 models


# ── Shared helpers ────────────────────────────────────────────────────────────

def _conditions(result, expected_rows=None):
    """Assert result is a PTT dict, return Conditions DataFrame."""
    assert isinstance(result, dict), f"Expected dict, got {type(result)}"
    assert "Conditions" in result, (
        f"'Conditions' missing from output; keys: {list(result)}"
    )
    cond = result["Conditions"]
    assert isinstance(cond, pd.DataFrame)
    if expected_rows is not None:
        assert cond.shape[0] == expected_rows, (
            f"Conditions has {cond.shape[0]} rows, expected {expected_rows}"
        )
    return cond


def _divided(result, expected_n=None):
    """
    Assert result is a dict-of-PTT-dicts (grid / multi-override output).
    Each sub-table must contain a 'Conditions' key.
    """
    assert isinstance(result, dict), f"Expected dict, got {type(result)}"
    if expected_n is not None:
        assert len(result) == expected_n, (
            f"Expected {expected_n} sub-tables, got {len(result)}: {list(result)}"
        )
    for key, sub in result.items():
        assert isinstance(sub, dict), f"Sub-table '{key}' is not a dict"
        assert "Conditions" in sub, f"Sub-table '{key}' missing 'Conditions'"
    return result


# ═════════════════════════════════════════════════════════════════════════════
# 1. Single-point (scalar P, T)
# ═════════════════════════════════════════════════════════════════════════════

class TestSinglePoint:
    """One P, one T → Conditions must have exactly 1 row."""

    def test_dict_comp(self):
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1200.0, P_bar=5000.0)
        _conditions(res, expected_rows=1)

    def test_series_comp(self):
        res = multi_path(Model=MODEL, comp=pd.Series(BASALT),
                         T_C=1200.0, P_bar=5000.0)
        _conditions(res, expected_rows=1)

    def test_single_element_array_PT(self):
        """1-element numpy arrays treated as scalars."""
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=np.array([1200.0]), P_bar=np.array([5000.0]))
        _conditions(res, expected_rows=1)

    def test_low_T(self):
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1050.0, P_bar=2000.0)
        _conditions(res, expected_rows=1)

    def test_high_P(self):
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1300.0, P_bar=20000.0)
        _conditions(res, expected_rows=1)


# ═════════════════════════════════════════════════════════════════════════════
# 2. PT path mode  (len(P) == len(T) > 1)
# ═════════════════════════════════════════════════════════════════════════════

class TestPTPath:
    """Same-length P and T arrays → PT path; Conditions row count == path length."""

    N = 6

    def test_isobaric_cooling(self):
        T = np.linspace(1300.0, 1050.0, self.N)
        P = np.full(self.N, 5000.0)
        res = multi_path(Model=MODEL, comp=BASALT, T_path_C=T, P_path_bar=P)
        _conditions(res, expected_rows=self.N)

    def test_decompression_path(self):
        P = np.linspace(10000.0, 2000.0, self.N)
        T = np.full(self.N, 1200.0)
        res = multi_path(Model=MODEL, comp=BASALT, T_path_C=T, P_path_bar=P)
        _conditions(res, expected_rows=self.N)

    def test_series_comp_on_path(self):
        T = np.linspace(1300.0, 1100.0, self.N)
        P = np.full(self.N, 3000.0)
        res = multi_path(Model=MODEL, comp=pd.Series(BASALT),
                         T_path_C=T, P_path_bar=P)
        _conditions(res, expected_rows=self.N)

    def test_T_start_end_dt_syntax(self):
        """1300 → 1050 °C in 50-step increments = 6 points."""
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_start_C=1300.0, T_end_C=1050.0, dt_C=50.0,
                         P_bar=5000.0)
        _conditions(res, expected_rows=6)

    def test_P_start_end_dp_syntax(self):
        """10000 → 2000 bar in 2000-bar steps = 5 points."""
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0,
                         P_start_bar=10000.0, P_end_bar=2000.0, dp_bar=2000.0)
        _conditions(res, expected_rows=5)

    def test_two_point_path(self):
        """Minimal path (2 points) should not error."""
        T = np.array([1300.0, 1100.0])
        P = np.array([5000.0, 5000.0])
        res = multi_path(Model=MODEL, comp=BASALT, T_path_C=T, P_path_bar=P)
        _conditions(res, expected_rows=2)


# ═════════════════════════════════════════════════════════════════════════════
# 3. P×T grid mode  (len(P) != len(T))
# ═════════════════════════════════════════════════════════════════════════════

class TestPTGrid:
    """
    Different-length P and T arrays → output divided by pressure.
    len(result) == len(P_vals); each sub-table has len(T_vals) rows.
    """

    P_VALS = np.array([3000.0, 5000.0, 8000.0])
    T_VALS = np.linspace(1100.0, 1300.0, 4)

    def test_grid_is_divided(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         P_bar=self.P_VALS, T_C=self.T_VALS)
        _divided(res, expected_n=len(self.P_VALS))

    def test_grid_sub_table_row_count(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         P_bar=self.P_VALS, T_C=self.T_VALS)
        for key, sub in res.items():
            assert sub["Conditions"].shape[0] == len(self.T_VALS), (
                f"Sub-table '{key}': expected {len(self.T_VALS)} rows, "
                f"got {sub['Conditions'].shape[0]}"
            )

    def test_2x2_grid(self):
        """Minimal 2×2 grid."""
        P = np.array([5000.0, 10000.0])
        T = np.array([1100.0, 1300.0])
        res = multi_path(Model=MODEL, comp=BASALT, P_bar=P, T_C=T)
        _divided(res, expected_n=2)
        for sub in res.values():
            _conditions(sub, expected_rows=2)

    def test_series_comp_on_grid(self):
        res = multi_path(Model=MODEL, comp=pd.Series(BASALT),
                         P_bar=self.P_VALS, T_C=self.T_VALS)
        _divided(res, expected_n=len(self.P_VALS))


# ═════════════════════════════════════════════════════════════════════════════
# 4. Parallel (DataFrame comp) mode
# ═════════════════════════════════════════════════════════════════════════════

class TestParallelMode:
    """DataFrame comp → one output row per input row (flat PTT dict)."""

    N = 4

    def test_df_matching_PT_arrays(self):
        """N-row DataFrame + N-element P and T → N rows in Conditions."""
        comp = pd.DataFrame([BASALT] * self.N)
        P = np.full(self.N, 5000.0)
        T = np.linspace(1100.0, 1300.0, self.N)
        res = multi_path(Model=MODEL, comp=comp, P_bar=P, T_C=T)
        _conditions(res, expected_rows=self.N)

    def test_df_broadcast_scalar_PT(self):
        """N-row DataFrame + scalar P and T → N rows (P/T broadcast)."""
        comp = pd.DataFrame([BASALT] * self.N)
        res = multi_path(Model=MODEL, comp=comp, P_bar=5000.0, T_C=1200.0)
        _conditions(res, expected_rows=self.N)

    def test_df_single_row(self):
        """1-row DataFrame → 1 row in Conditions."""
        comp = pd.DataFrame([BASALT])
        res = multi_path(Model=MODEL, comp=comp, P_bar=5000.0, T_C=1200.0)
        _conditions(res, expected_rows=1)

    def test_df_variable_compositions(self):
        """DataFrame with differing compositions (not duplicates) → N rows."""
        rows = [dict(BASALT) for _ in range(3)]
        rows[0]["H2O_Liq"] = 0.5
        rows[1]["H2O_Liq"] = 1.0
        rows[2]["H2O_Liq"] = 2.0
        comp = pd.DataFrame(rows)
        res = multi_path(Model=MODEL, comp=comp, P_bar=5000.0, T_C=1200.0)
        _conditions(res, expected_rows=3)


# ═════════════════════════════════════════════════════════════════════════════
# 5. Volatile / composition override kwargs
# ═════════════════════════════════════════════════════════════════════════════

class TestOverrides:
    """Scalar override kwargs replace the matching oxide in the composition."""

    def test_H2O_init_scalar(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, H2O_init=1.5)
        _conditions(res, expected_rows=1)

    def test_CO2_init_scalar(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, CO2_init=0.1)
        _conditions(res, expected_rows=1)

    def test_Fe3Fet_init_scalar(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, Fe3Fet_init=0.2)
        _conditions(res, expected_rows=1)

    def test_H2O_init_on_path(self):
        """Scalar H2O override applied uniformly along a PT path."""
        T = np.linspace(1300.0, 1100.0, 5)
        P = np.full(5, 5000.0)
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_path_C=T, P_path_bar=P, H2O_init=2.0)
        _conditions(res, expected_rows=5)

    def test_H2O_and_CO2_together(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         H2O_init=1.0, CO2_init=0.05)
        _conditions(res, expected_rows=1)


# ═════════════════════════════════════════════════════════════════════════════
# 6. fO2 constraints
# ═════════════════════════════════════════════════════════════════════════════

class TestFO2:

    def test_fO2_buffer_FMQ_zero_offset(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         fO2_buffer="FMQ", fO2_offset=0.0)
        _conditions(res, expected_rows=1)

    def test_fO2_buffer_QFM_zero_offset(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         fO2_buffer="QFM", fO2_offset=0.0)
        _conditions(res, expected_rows=1)

    def test_fO2_positive_offset(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         fO2_buffer="FMQ", fO2_offset=1.5)
        _conditions(res, expected_rows=1)

    def test_fO2_negative_offset(self):
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         fO2_buffer="FMQ", fO2_offset=-1.0)
        _conditions(res, expected_rows=1)

    def test_fO2_on_PT_path(self):
        T = np.linspace(1300.0, 1100.0, 6)
        P = np.full(6, 5000.0)
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_path_C=T, P_path_bar=P,
                         fO2_buffer="FMQ", fO2_offset=0.0)
        _conditions(res, expected_rows=6)

    def test_fO2_on_PT_grid(self):
        """Scalar fO2 on a P×T grid → flat table (no fO2-axis grouping)."""
        P = np.array([5000.0, 10000.0])
        T = np.array([1100.0, 1200.0, 1300.0])
        res = multi_path(Model=MODEL, comp=BASALT,
                         P_bar=P, T_C=T,
                         fO2_buffer="FMQ", fO2_offset=0.5)
        # scalar fO2 does not trigger fO2-axis grouping → flat dict, P*T rows
        assert isinstance(res, dict)
        assert "Conditions" in res
        assert res["Conditions"].shape[0] == len(P) * len(T)


# ═════════════════════════════════════════════════════════════════════════════
# 7. Multi-valued overrides → composition grid
# ═════════════════════════════════════════════════════════════════════════════

class TestCompositionGrid:
    """
    A list passed to H2O_init, CO2_init, or Fe3Fet_init triggers a grid
    calculation.  Output is a dict-of-PTT-dicts (one sub-table per value).
    """

    def test_H2O_multi_single_PT(self):
        """3 H2O values × 1 PT point → 3 sub-tables, each with 1 row."""
        vals = [0.0, 1.0, 2.0]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, H2O_init=vals)
        _divided(res, expected_n=len(vals))
        for sub in res.values():
            _conditions(sub, expected_rows=1)

    def test_H2O_multi_four_values(self):
        vals = [0.0, 0.5, 1.0, 2.0]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, H2O_init=vals)
        _divided(res, expected_n=len(vals))

    def test_CO2_multi(self):
        vals = [0.0, 0.05, 0.1]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, CO2_init=vals)
        _divided(res, expected_n=len(vals))

    def test_Fe3Fet_multi(self):
        vals = [0.05, 0.15, 0.25]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0, Fe3Fet_init=vals)
        _divided(res, expected_n=len(vals))

    def test_H2O_multi_on_PT_path(self):
        """
        3 H2O values × 3-step path → 3 sub-tables, each with 3 rows.
        Total ForwardMB rows = 3 × 3 = 9.
        """
        T = np.linspace(1300.0, 1150.0, 3)
        P = np.full(3, 5000.0)
        vals = [0.5, 1.0, 2.0]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_path_C=T, P_path_bar=P, H2O_init=vals)
        _divided(res, expected_n=len(vals))
        for sub in res.values():
            _conditions(sub, expected_rows=len(T))


# ═════════════════════════════════════════════════════════════════════════════
# 8. fO2 as a grid axis
# ═════════════════════════════════════════════════════════════════════════════

class TestFO2Grid:

    def test_multi_fO2_on_single_PT(self):
        """3 fO2 offsets × 1 PT point → divided by fO2 (3 sub-tables, 1 row each)."""
        offsets = [-1.0, 0.0, 1.0]
        res = multi_path(Model=MODEL, comp=BASALT,
                         T_C=1200.0, P_bar=5000.0,
                         fO2_buffer="FMQ", fO2_offset=offsets)
        # fO2 as a path axis: len(P)==len(T)==1, len(fO2)==3 → PT path at 3 fO2 values
        # output structure depends on run_nGibbs branch; at minimum: check it runs
        assert isinstance(res, dict)

    def test_multi_fO2_on_PT_grid(self):
        """P×T grid with 3 fO2 values → divided by fO2 (3 sub-tables)."""
        P = np.array([5000.0, 10000.0])
        T = np.array([1100.0, 1200.0, 1300.0])
        offsets = [-1.0, 0.0, 1.0]
        res = multi_path(Model=MODEL, comp=BASALT,
                         P_bar=P, T_C=T,
                         fO2_buffer="FMQ", fO2_offset=offsets)
        _divided(res, expected_n=len(offsets))
        for sub in res.values():
            _conditions(sub, expected_rows=len(P) * len(T))


# ═════════════════════════════════════════════════════════════════════════════
# 9. Model variants
# ═════════════════════════════════════════════════════════════════════════════

class TestModelVariants:

    @pytest.mark.skip(reason="MELTS 1.2.0 trained models not yet available")
    def test_nMELTS120_NoProp(self):
        res = multi_path(Model="nMELTSv1.2.0NoProp", comp=BASALT,
                         T_C=1200.0, P_bar=5000.0)
        _conditions(res, expected_rows=1)

    def test_nMELTS102_NoProp(self):
        res = multi_path(Model="nMELTSv1.0.2NoProp", comp=BASALT,
                         T_C=1200.0, P_bar=5000.0)
        _conditions(res, expected_rows=1)

    @pytest.mark.slow
    def test_nMELTS120_with_props(self):
        """nMELTSv1.2.0 runs calc_phase_props_MELTS after ForwardMB."""
        pytest.importorskip("meltsdynamic", reason="meltsdynamic not installed")
        res = multi_path(Model="nMELTSv1.2.0", comp=BASALT,
                         T_C=1200.0, P_bar=5000.0)
        _conditions(res, expected_rows=1)


# ═════════════════════════════════════════════════════════════════════════════
# 10. Output structure spot-checks
# ═════════════════════════════════════════════════════════════════════════════

class TestOutputStructure:
    """Verify expected keys and DataFrame properties in nGibbs output."""

    def test_conditions_is_dataframe(self):
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1200.0, P_bar=5000.0)
        assert isinstance(res["Conditions"], pd.DataFrame)

    def test_conditions_has_columns(self):
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1200.0, P_bar=5000.0)
        assert res["Conditions"].shape[1] > 0

    def test_liquid_table_present(self):
        """At least one liquid composition table should appear in output."""
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1200.0, P_bar=5000.0)
        liq_keys = [k for k in res if "liq" in k.lower()]
        assert len(liq_keys) > 0, f"No liquid table found; keys: {list(res)}"

    def test_mass_key_present(self):
        """mass_g or similar mass key should exist."""
        res = multi_path(Model=MODEL, comp=BASALT, T_C=1200.0, P_bar=5000.0)
        mass_keys = [k for k in res if "mass" in k.lower()]
        assert len(mass_keys) > 0, f"No mass key found; keys: {list(res)}"

    def test_conditions_row_count_matches_path(self):
        """Spot-check that Conditions rows == path length across two calls."""
        for n in [3, 7]:
            T = np.linspace(1300.0, 1100.0, n)
            P = np.full(n, 5000.0)
            res = multi_path(Model=MODEL, comp=BASALT, T_path_C=T, P_path_bar=P)
            assert res["Conditions"].shape[0] == n


# ═════════════════════════════════════════════════════════════════════════════
# 11. Error / guard tests
# ═════════════════════════════════════════════════════════════════════════════

class TestGuards:
    """Verify that invalid inputs raise appropriate errors."""

    def test_unsupported_fO2_buffer_raises(self):
        """NNO buffer is not supported by nGibbs → ValueError."""
        with pytest.raises(ValueError):
            multi_path(Model=MODEL, comp=BASALT,
                       T_C=1200.0, P_bar=5000.0, fO2_buffer="NNO")

    def test_simultaneous_fO2_and_Fe3Fet_raises(self):
        """fO2_buffer + Fe3Fet_init in overrides are mutually exclusive."""
        with pytest.raises(ValueError, match="mutually exclusive"):
            multi_path(Model=MODEL, comp=BASALT,
                       T_C=1200.0, P_bar=5000.0,
                       fO2_buffer="FMQ", fO2_offset=0.0,
                       Fe3Fet_init=0.1)

    def test_two_multi_valued_comp_overrides_raise(self):
        """Two multi-valued compositional overrides must raise ValueError."""
        with pytest.raises(ValueError):
            multi_path(Model=MODEL, comp=BASALT,
                       T_C=1200.0, P_bar=5000.0,
                       H2O_init=[0.5, 1.0], CO2_init=[0.0, 0.1])

    def test_missing_composition_raises(self):
        """Passing comp=None should raise ValueError in nGibbsAPI."""
        with pytest.raises(ValueError):
            multi_path(Model=MODEL, comp=None, T_C=1200.0, P_bar=5000.0)

    def test_no_temperature_raises(self):
        """No temperature specified → ValueError from build_PT_vectors."""
        with pytest.raises(ValueError):
            multi_path(Model=MODEL, comp=BASALT, P_bar=5000.0)

    def test_no_pressure_raises(self):
        """No pressure specified → ValueError from build_PT_vectors."""
        with pytest.raises(ValueError):
            multi_path(Model=MODEL, comp=BASALT, T_C=1200.0)
