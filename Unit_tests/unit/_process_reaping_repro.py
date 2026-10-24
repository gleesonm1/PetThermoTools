"""
Helper script for test_process_reaping.py -- NOT a test module (the leading
underscore keeps pytest from collecting it).

Runs the actual check in its own interpreter, so a hang can be bounded by
killing the whole process tree from outside -- matching every other
mp_utils-based repro in this directory (_multi_path_repro.py,
_timeout_repro.py, _no_guard_repro.py). This was the one test in that family
still calling multi_path directly inside the pytest process, with nothing to
catch a hang beyond CI's own 60-minute job timeout: confirmed repeatedly
hanging on ubuntu/3.10 in CI (never reproduced locally), landing here right
after the dead-monkeypatch fix made everything else in the file fast again.
Whatever the exact remaining cause, this bounds the damage to one failed test
instead of a stuck hour-long CI job.

Prints:
    STARTED=1                     just before multi_path is called
    RESULT=OK                     every worker exited and was joined
    RESULT=NO_WORKERS             multi_path returned but created no workers
    RESULT=STILL_ALIVE:<n>        n worker(s) still alive after multi_path returned
    RESULT=UNJOINED:<n>/<total>   n of total workers were never .join()ed

Usage:  python _process_reaping_repro.py
"""
import importlib
import os
import sys
from pathlib import Path

import numpy as np

# Test the working tree, not a stale installed copy (mirrors conftest.py).
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# `import petthermotools.Path as ptt_path` resolves to pathlib.Path instead --
# petthermotools/__init__.py's own `from .Path import *` shadows the name
# (see test_import_hygiene.py). Every other repro script in this directory
# uses importlib.import_module for the same reason.
ptt_path = importlib.import_module("petthermotools.Path")
ptt_mp_utils = importlib.import_module("petthermotools.mp_utils")

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


# Under spawn, Process.start() pickles the WHOLE Process instance (not just
# target/args/kwargs) to hand to the child, so this has to be a real,
# picklable module-level class. It subclasses mp_utils._CTX.Process as it
# stands at import time (the real spawn-context Process class, e.g.
# SpawnProcess) rather than the platform-default `multiprocessing.Process`
# -- those are not always the same class (Linux's bare `multiprocessing.
# Process` defaults to ForkProcess), and swapping in a fork-flavoured class
# under a spawn-context Queue/pipe setup would be a real, separate bug in
# this script, not the thing under test.
class _JoinTrackingProcess(ptt_mp_utils._CTX.Process):
    """A drop-in replacement that records every instance created and
    whether .join() was ever called on it. Only the parent-side wrapper is
    subclassed -- the target/args/kwargs sent to the child across the spawn
    boundary are unaffected, so this does not change what actually runs in
    the worker."""

    created = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.join_called = False
        _JoinTrackingProcess.created.append(self)

    def join(self, *args, **kwargs):
        self.join_called = True
        return super().join(*args, **kwargs)


def _fake_worker_completes_instantly(q, task, **kwargs):
    """Fake _path_task_v2: reports one step for the row and exits normally
    -- i.e. it always takes the natural-exit branch, never the timeout
    branch, isolating Mechanism 1 from the terminate()/timeout behaviour
    covered separately in test_timeout_handling.py."""
    q.put([task.row_id, 0, {"Conditions": {}, "_dummy_keys": []}])


def main():
    ptt_mp_utils._CTX.Process = _JoinTrackingProcess
    ptt_path._path_task_v2 = _fake_worker_completes_instantly

    print("STARTED=1", flush=True)
    ptt_path.multi_path(
        Model="MELTSv1.0.2",
        comp=COMP,
        T_C=1200.0,
        P_bar=np.array([1000.0, 2000.0]),  # one row per core -- see
        # _multi_path_repro.py's identical comment for why n == cores
        # matters under multi_path v2's one-row-per-process model.
        cores=2,
        timeout=5,
        multi_processing=True,
        Print_suppress=True,
    )

    created = _JoinTrackingProcess.created
    if not created:
        result = "NO_WORKERS"
    else:
        still_alive = sum(1 for p in created if p.is_alive())
        if still_alive:
            result = f"STILL_ALIVE:{still_alive}"
        else:
            unjoined = sum(1 for p in created if not p.join_called)
            result = "OK" if unjoined == 0 else f"UNJOINED:{unjoined}/{len(created)}"

    print(f"RESULT={result}", flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
