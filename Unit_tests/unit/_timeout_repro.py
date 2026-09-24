"""
Helper script for test_timeout_handling.py -- NOT a test module (the leading
underscore keeps pytest from collecting it).

Runs multi_path against a fake worker that hangs, in its own interpreter, and
reports what multi_path did. Each scenario is one (kind, worker) pair:

  kind    'batch'  -> P_bar is an array; workers are path_multi(q, index_array)
          'single' -> P_bar is a scalar; the worker is path(q, 1)
  worker  follows the real queue protocol: batch workers send a
          [index, -1, "Start"] marker before each run, single-run workers send
          step messages only.

Usage:  python _timeout_repro.py <scenario>

Prints (one per line) so the test can parse them:
    STARTED=1     just before multi_path is called (after all imports)
    OUTCOME=RETURNED | RAISED
    EXC=<exception class name, or empty>
    EMPTY=1|0     whether a returned result was empty
    ELAPSED=<seconds spent inside multi_path>
    ALIVE=<comma-separated PIDs of workers still running afterwards>
    SIGNALS=<warnings raised during the call, ' || ' separated>

Exits via os._exit so Python's exit-time join of non-daemon children can't
hang the script; the test kills whatever ALIVE reports.
"""
import importlib
import multiprocessing
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np

# Test the working tree, not a stale installed copy (mirrors conftest.py).
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}
_STEP = {"Conditions": {}, "_dummy_keys": []}

# Finite so a worker that escapes cleanup cannot live forever.
_HANG_S = 120

# multi_path's own timeout. Large enough that a healthy worker has reported
# well before it expires even on a slow machine (workers are light: they only
# import numpy), small enough to keep the scenarios quick.
MULTI_PATH_TIMEOUT = 4


def _report_then_hang(q, index, **kwargs):
    """Batch worker: announces each run, then hangs mid-calculation."""
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
    time.sleep(_HANG_S)


def _step_then_hang_batch(q, index, **kwargs):
    """Batch worker: announces each run, delivers one step, then hangs."""
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
        q.put([i, 0, _STEP])
    time.sleep(_HANG_S)


def _hang_before_report(q, index, **kwargs):
    """Worker that hangs before sending anything (e.g. the engine hangs while
    initialising). Used for both batch and single runs."""
    time.sleep(_HANG_S)


def _step_then_hang_single(q, index, **kwargs):
    """Single-run worker: delivers one step (no Start marker), then hangs."""
    q.put([index, 0, _STEP])
    time.sleep(_HANG_S)


def _ignores_sigterm(q, index, **kwargs):
    """Batch worker stuck somewhere that never returns to the interpreter's
    signal handling (a long C call) -- modelled as ignoring SIGTERM, which is
    what Process.terminate() sends on macOS/Linux."""
    import signal

    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
    time.sleep(_HANG_S)


SCENARIOS = {
    "batch_report_then_hang": ("batch", _report_then_hang),
    "batch_step_then_hang": ("batch", _step_then_hang_batch),
    "batch_hang_before_report": ("batch", _hang_before_report),
    "batch_ignores_sigterm": ("batch", _ignores_sigterm),
    "single_hang_before_report": ("single", _hang_before_report),
    "single_step_then_hang": ("single", _step_then_hang_single),
}

# Incidental library noise that says nothing about how timeouts are reported.
_NOISE = (DeprecationWarning, PendingDeprecationWarning, FutureWarning, SyntaxWarning, ResourceWarning)


def main(name):
    kind, worker = SCENARIOS[name]
    ptt_path = importlib.import_module("petthermotools.Path")
    ptt_path.path_multi = worker  # what batch runs spawn
    ptt_path.path = worker        # what single runs spawn

    kwargs = {
        'Model': "MELTSv1.0.2", 'comp': COMP, 'T_C': 1200.0, 'cores': 2,
        'timeout': MULTI_PATH_TIMEOUT, 'multi_processing': True, 'Print_suppress': True,
        'P_bar': np.array([1000.0, 2000.0, 3000.0, 4000.0, 5000.0]) if kind == "batch" else 1000.0,
    }

    print("STARTED=1", flush=True)
    outcome, exc_name, empty = "RETURNED", "", 0
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = ptt_path.multi_path(**kwargs)
            empty = int(hasattr(result, "__len__") and len(result) == 0)
        except BaseException as exc:  # noqa: BLE001
            outcome, exc_name = "RAISED", type(exc).__name__
    elapsed = time.time() - t0

    time.sleep(0.5)  # let any cleanup multi_path did take effect
    alive = [p.pid for p in multiprocessing.active_children() if p.is_alive()]
    signals = [
        f"{w.category.__name__}: {' '.join(str(w.message).split())[:80]}"
        for w in caught if not issubclass(w.category, _NOISE)
    ]
    print(f"OUTCOME={outcome}")
    print(f"EXC={exc_name}")
    print(f"EMPTY={empty}")
    print(f"ELAPSED={elapsed:.1f}")
    print("ALIVE=" + ",".join(str(pid) for pid in alive))
    print("SIGNALS=" + " || ".join(signals), flush=True)
    os._exit(0)


if __name__ == "__main__":
    main(sys.argv[1])
