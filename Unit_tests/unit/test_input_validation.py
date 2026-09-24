"""
Tier A failure-injection tests: bad `Model` strings and bad compositions
passed to multi_path.

Safety: multi_path picks a backend from the Model string (any name without
"MELTS" is treated as MAGEMin and boots Julia; a MELTS name spawns worker
processes). These tests must never reach either, so the `tripwired_backends`
fixture replaces the Julia start-up helpers and `Process` with tripwires that
raise `BackendTouched`. A test therefore observes only what multi_path does
in the parent, before any backend work, and cannot start Julia or MELTS on
any machine (or hang on CI, where neither is installed).

Three groups:
  1. Routing that works today and must survive the backend-registry refactor
     (valid model names reach the right backend).
  2. Clear errors that already exist today and must not regress.
  3. Proposed behaviour, marked xfail: today bad input is silently routed to
     the wrong backend, silently accepted, or fails with an obscure error
     (ZeroDivisionError, AttributeError: ... 'copy', ...). These encode what
     the cleanup plan proposes; they are a proposal, not an existing spec.
     `raises=AssertionError` means only a genuine mismatch counts as the
     expected failure.
"""
import importlib

import numpy as np
import pandas as pd
import pytest

# See test_import_hygiene.py: `import petthermotools.Path` is shadowed by
# pathlib.Path via the package's own star-import.
ptt_path = importlib.import_module("petthermotools.Path")

GOOD = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


class BackendTouched(Exception):
    """multi_path got far enough to start a backend."""


@pytest.fixture(autouse=True)
def tripwired_backends(monkeypatch):
    def _julia(*args, **kwargs):
        raise BackendTouched("JULIA")

    class _NoSpawn:
        def __init__(self, *args, **kwargs):
            raise BackendTouched("SPAWN")

    monkeypatch.setattr(ptt_path, "_ensure_julia_ready", _julia)
    monkeypatch.setattr(ptt_path, "_ensure_julia_workers", _julia)
    monkeypatch.setattr(ptt_path, "Process", _NoSpawn)


def _call(**overrides):
    kwargs = {'Model': "MELTSv1.0.2", 'comp': GOOD, 'T_C': 1200.0, 'P_bar': 1000.0,
              'Print_suppress': True, 'timeout': 5}
    kwargs.update(overrides)
    return ptt_path.multi_path(**kwargs)


def _backend_reached(**overrides):
    """Returns 'JULIA' / 'SPAWN' for the backend multi_path tried to start."""
    with pytest.raises(BackendTouched) as info:
        _call(**overrides)
    return str(info.value)


def _assert_clear_error(mention, expected=(ValueError, TypeError), **overrides):
    """Passes only if multi_path rejects the input with a ValueError/TypeError
    whose message mentions `mention`, without touching a backend."""
    try:
        _call(**overrides)
    except BackendTouched as exc:
        raise AssertionError(f"reached the {exc} backend instead of rejecting the input") from None
    except expected as exc:
        if mention.lower() not in str(exc).lower():
            raise AssertionError(
                f"raised {type(exc).__name__} but its message {str(exc)!r} does not mention {mention!r}"
            ) from None
        return
    except BaseException as exc:  # noqa: BLE001
        raise AssertionError(
            f"raised {type(exc).__name__}: {exc} instead of a clear ValueError/TypeError"
        ) from None
    raise AssertionError("returned without raising")


# --------------------------------------------------------------------------
# 1. Routing that works today (pins behaviour for the backend-registry work)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("model, backend", [
    ("MELTSv1.0.2", "SPAWN"),
    ("MELTSv1.1.0", "SPAWN"),
    ("MELTSv1.2.0", "SPAWN"),
    ("pMELTS", "SPAWN"),
    ("Green2025", "JULIA"),
    ("Weller2024", "JULIA"),
])
def test_valid_models_route_to_the_expected_backend(model, backend):
    assert _backend_reached(Model=model) == backend


def test_model_none_defaults_to_a_melts_model_in_multi_path():
    """multi_path defaults Model=None to MELTSv1.0.2. (Other entry points,
    e.g. Liq.py, test `"MELTS" in Model` first and raise TypeError instead.)"""
    assert _backend_reached(Model=None) == "SPAWN"


# --------------------------------------------------------------------------
# 2. Clear errors that already exist and must not regress
# --------------------------------------------------------------------------

def test_zero_h2o_is_rejected_with_a_hint_for_melts():
    with pytest.raises(Exception, match="H2O"):
        _call(comp={**GOOD, 'H2O_Liq': 0.0})


def test_zero_ferric_iron_without_a_buffer_is_rejected_with_a_hint():
    with pytest.raises(Exception, match="ferric"):
        _call(comp={**GOOD, 'Fe3Fet_Liq': 0.0})


def test_unknown_fo2_buffer_is_rejected():
    with pytest.raises(Exception, match="fO2 buffer"):
        _call(fO2_buffer="XYZ")


def test_unknown_label_is_rejected():
    with pytest.raises(ValueError, match="label"):
        _call(label="bogus")


def test_mismatched_condition_lengths_are_rejected():
    with pytest.raises(ValueError, match="non-identical length"):
        _call(P_bar=np.array([1000.0, 2000.0, 3000.0]), comp=pd.DataFrame([GOOD, GOOD]))


def test_non_numeric_oxide_value_is_rejected():
    with pytest.raises(ValueError):
        _call(comp={**GOOD, 'SiO2_Liq': 'abc'})


# --------------------------------------------------------------------------
# 3. Proposed behaviour (xfail today)
# --------------------------------------------------------------------------

_XFAIL = pytest.mark.xfail(
    raises=AssertionError,  # only a genuine mismatch counts as "expected"
    reason=(
        "Not validated today: bad input is silently routed to a backend, "
        "silently accepted, or fails with an obscure error. See the cleanup "
        "plan (Part (i) #3 validate_inputs / #5 backend registry)."
    ),
    strict=False,
)


@_XFAIL
@pytest.mark.parametrize("model", [
    "NotAModel",   # today: silently treated as MAGEMin, boots Julia
    "melts",       # lowercase typo: same
    "",            # empty string: same
    "MELTSv9.9",   # bad MELTS version: spawns a worker that then fails
], ids=["unknown-name", "lowercase-typo", "empty", "bad-melts-version"])
def test_unknown_model_is_rejected_before_any_backend(model):
    _assert_clear_error("Model", Model=model)


@_XFAIL
def test_non_string_model_gets_a_clear_error():
    # today: TypeError "argument of type 'int' is not iterable"
    _assert_clear_error("Model", Model=123)


@_XFAIL
@pytest.mark.parametrize("comp", [
    {},                                      # today: ZeroDivisionError
    None,                                    # today: AttributeError ... 'copy'
    [1, 2, 3],                               # today: AttributeError ... 'columns'
    {**GOOD, 'SiO2_Liq': float('nan')},      # today: accepted, sent to a worker
    {**GOOD, 'SiO2_Liq': -5.0},              # today: accepted (renormalised to nonsense)
    {k: 0.0 for k in GOOD},                  # today: ZeroDivisionError
], ids=["empty", "none", "list", "nan", "negative", "all-zero"])
def test_bad_composition_gets_a_clear_error(comp):
    _assert_clear_error("comp", comp=comp)


@_XFAIL
def test_no_temperature_or_pressure_is_rejected():
    # today: spawns a worker with nothing to calculate
    _assert_clear_error("", T_C=None, P_bar=None)
