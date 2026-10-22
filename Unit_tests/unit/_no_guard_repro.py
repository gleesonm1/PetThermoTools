"""
Helper script for test_missing_main_guard.py -- NOT a test module (the leading
underscore keeps pytest from collecting it).

Behaves like a typical user script that calls a petthermotools calculation at
the top level of the file, with a fake healthy worker standing in for MELTS.

    python _no_guard_repro.py            # calculation at module level: NO guard
    python _no_guard_repro.py guarded    # same calculation under
                                         #   `if __name__ == "__main__":`

Why the guard matters: on Windows and macOS (and for `spawn`/`forkserver`
anywhere) every worker is a fresh interpreter that re-imports the main script.
Without the guard, each worker re-runs the calculation while it is still
importing, and Python refuses ("An attempt has been made to start a new
process before the current process has finished its bootstrapping phase").
The worker dies before reporting anything.

`spawn` is selected explicitly so Linux (which forks by default, hiding the
problem) behaves like Windows/macOS.

Prints:
    STARTED=1    just before the calculation (only the top-level process prints
                 it as a marker the test waits for; workers may repeat it)
    RETURNED=1   after the calculation returned (never reached by a worker)
"""
import importlib
import multiprocessing
import sys
from pathlib import Path

import numpy as np

# Test the working tree, not a stale installed copy (mirrors conftest.py).
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

multiprocessing.set_start_method("spawn", force=True)

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


def _healthy_worker(q, task, **kwargs):
    """Fake _path_task_v2: reports one step for the row and exits normally.

    multi_path (v2, the default) calls _path_task_v2(q, task, **kwargs)
    once per row -- not path_multi(q, index_array) -- so patching
    path_multi here was inert and let every row fall through to the real
    _path_task_v2 (meltsdynamic import, absent in CI), making this
    script's timing depend on real-engine failure behaviour instead of a
    fast, deterministic fake. Same fix as _timeout_repro.py.

    No "-1, 'Start'" marker here (v1's path_multi protocol): v2's
    mp_utils._worker_loop already reports a worker's liveness to the
    supervisor on its own (_STARTED, sent before task_fn even runs), and
    multi_path's _capture stores every item task_fn puts straight into
    `collected` with no v1-style filtering, so a bare "Start" string ends
    up handed to compile_results, which crashes calling .items() on it.
    """
    q.put([task.row_id, 0, {"Conditions": {}, "_dummy_keys": []}])


def run_calculation():
    ptt_path = importlib.import_module("petthermotools.Path")
    ptt_path._path_task_v2 = _healthy_worker
    print("STARTED=1", flush=True)
    ptt_path.multi_path(
        Model="MELTSv1.0.2", comp=COMP, T_C=1200.0,
        P_bar=np.array([1000.0, 2000.0]),
        cores=2, timeout=30, multi_processing=True, Print_suppress=True,
    )
    print("RETURNED=1", flush=True)


if len(sys.argv) > 1 and sys.argv[1] == "guarded":
    if __name__ == "__main__":
        run_calculation()
else:
    run_calculation()  # <- no guard: this is what a user's script often looks like
