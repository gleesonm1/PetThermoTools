"""
Helper script for test_multi_path_supervisor.py -- NOT a test module
(the leading underscore keeps pytest from collecting it).

Runs multi_path against a fake worker in its *own interpreter*, so that a
genuine hang can be bounded by killing the whole process tree. Running it
on a thread inside the pytest process instead leaks the still-running
respawn loop into every later test: the loop looks up `Process` and
`path_multi` as module globals on every round, so it silently picks up
whatever a later test has patched them to (or the real, un-patched
`path_multi` once the patch is undone).

Usage:  python _multi_path_repro.py never_reports | completes

The `if __name__ == "__main__"` guard is required: under the spawn start
method (Windows/macOS) each worker re-imports this file.
"""
import importlib
import sys
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


def _worker_never_reports(q, index, **kwargs):
    """Engine construction fails before anything is put on the queue -- the
    alphaMELTS 'C library will not reinitialize after a crash' failure mode
    (the real path_multi has no try/except around MELTSdynamic(...))."""
    raise RuntimeError("simulated permanent alphaMELTS engine init failure")


def _worker_completes(q, index, **kwargs):
    """A well-behaved worker following path_multi's real queue protocol: a
    'Start' marker, then step messages, then a normal exit."""
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
        q.put([i, 0, {"Conditions": {}, "_dummy_keys": []}])


_WORKERS = {"never_reports": _worker_never_reports, "completes": _worker_completes}


def main(mode):
    ptt_path = importlib.import_module("petthermotools.Path")
    ptt_path.path_multi = _WORKERS[mode]
    ptt_path.multi_path(
        Model="MELTSv1.0.2",
        comp=COMP,
        T_C=1200.0,
        P_bar=np.array([1000.0, 2000.0, 3000.0, 4000.0, 5000.0]),
        cores=2,
        # Generous on purpose: this only bounds a *hung* worker; a well-behaved
        # one exits as soon as it has reported. (Workers here only import numpy
        # -- petthermotools is imported inside main(), not at module level -- so
        # they start quickly.)
        timeout=60,
        multi_processing=True,
        Print_suppress=True,
    )


if __name__ == "__main__":
    main(sys.argv[1])
